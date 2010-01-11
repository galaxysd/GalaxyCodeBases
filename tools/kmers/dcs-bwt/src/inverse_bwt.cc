// Copyright 2007 Google Inc.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

// Two algorithms for the inverse Burrows-Wheeler transform.
// 1. SimpleInverseBWTransformer: This is the algorithm basis in
//    Seward: Space-time tradeoffs in the Inverse B-W Transform.
// 2. FastInverseBWTransformer: This is the algorithm mergeTL
//    in Sewards paper with some additional heuristics to be
//    able to handle very large blocks.

#include "inverse_bwt.h"
#include "stream.h"

#include <iostream>
#include <string>
#include <vector>
#include <numeric>  // for partial_sum

namespace dcsbwt {

namespace {

//////////////////////////////////////////////////////////////////
// SimpleInverseBWTransformer
//////////////////////////////////////////////////////////////////

class SimpleInverseBWTransformer : public InverseBWTransformer {
 public:
  class Options : public InverseBWTransformer::AlgorithmSpecificOptions {
   public:
    Options() {}
    virtual ~Options() {}

    virtual bool Set(const std::string& options) {
      return 0 == options.length();
    }
    virtual std::string Get() const { return std::string(); }
    virtual int64 MaxBlockSize(int64 memory_budget) const {
      return memory_budget / 9;
    }

    virtual InverseBWTransformer* GetTransformer(int64 memory_budget) {
      return new SimpleInverseBWTransformer;
    }
  };

  SimpleInverseBWTransformer() {}
  virtual ~SimpleInverseBWTransformer() {}

 private:
  virtual bool DoTransform(const char* bwt, size_t bwt_size,
                            size_t eob_position, OutStream* output);

  SimpleInverseBWTransformer(const SimpleInverseBWTransformer&);
  SimpleInverseBWTransformer& operator=(const SimpleInverseBWTransformer&);
};


bool SimpleInverseBWTransformer::DoTransform(
    const char* bwt, size_t bwt_size,
    size_t eob_position, OutStream* output)
{
  // rank[i] will be the number of occurrences of bwt[i] in bwt[0..i-1]
  std::vector<uint32> rank(bwt_size, 0);

  // count[] serves two purposes:
  // 1. When the scan of the BWT reaches position i,
  //    count[c+1] is the number of occurrences of c in bwt[0..i-1]
  //    (count[0] counts the EOB symbol).
  // 2. During the ouput generation, count[c] is the total number
  //    of characters smaller than c (including the single EOB symbol).
  std::vector<uint32> count(257, 0);

  // count characters before EOB
  for (uint32 position = 0; position < eob_position; ++position) {
    int ch = static_cast<unsigned char>(bwt[position]);
    rank[position] = count[ch + 1];
    ++count[ch + 1];
  }
  // count EOB
  rank[eob_position] = 0;
  count[0] = 1;
  // count characters after EOB
  for (uint32 position = eob_position + 1; position < bwt_size; ++position) {
    int ch = static_cast<unsigned char>(bwt[position]);
    rank[position] = count[ch + 1];
    ++count[ch + 1];
  }
  std::partial_sum(count.begin(), count.end(), count.begin());
  assert(count[256] == bwt_size);

  OutStreamBuffer buffer;
  buffer.Connect(output);
  uint32 position = 0;
  uint32 string_length = 0;
  while (position != eob_position) {
    int ch = static_cast<unsigned char>(bwt[position]);
    buffer.WriteByte(ch);
    ++string_length;
    position = count[ch] + rank[position];
  }
  // If the BWT or the EOB position contain errors (or are garbage),
  // it is likely that the cycle ends prematurely.
  // In this case, we return false.
  buffer.Disconnect();
  return (string_length + 1 == bwt_size);
}

//////////////////////////////////////////////////////////////////
// FastInverseBWTransformer
//////////////////////////////////////////////////////////////////

class FastInverseBWTransformer : public InverseBWTransformer {
 public:
  class Options : public InverseBWTransformer::AlgorithmSpecificOptions {
   public:
    Options() {}
    virtual ~Options() {}

    virtual bool Set(const std::string& options) {
      return 0 == options.length();
    }
    virtual std::string Get() const { return std::string(); }
    virtual int64 MaxBlockSize(int64 memory_budget) const {
      memory_budget -= kMemoryOverhead;
      return memory_budget / sizeof(uint32);
    }

    virtual InverseBWTransformer* GetTransformer(int64 memory_budget) {
      return new FastInverseBWTransformer;
    }
   private:
    static const int64 kMemoryOverhead = 1 << 20;
  };

  FastInverseBWTransformer() {}
  virtual ~FastInverseBWTransformer() {}

 private:
  virtual bool DoTransform(const char* bwt, size_t bwt_size,
                            size_t eob_position, OutStream* output);

  FastInverseBWTransformer(const FastInverseBWTransformer&);
  FastInverseBWTransformer& operator=(const FastInverseBWTransformer&);
};

bool FastInverseBWTransformer::DoTransform(
    const char* bwt, size_t bwt_size,
    size_t eob_position, OutStream* output)
{
  // rank[i] will be the number of occurrences of bwt[i] in bwt[0..i-1]
  std::vector<uint32> bwt_rank_low24(bwt_size);
  std::vector<uint32> rank_milestone_buffer[256];

  // count[] serves two purposes:
  // 1. When the scan of the BWT reaches position i,
  //    count[c+1] is the number of occurrences of c in bwt[0..i-1]
  //    (count[0] counts the EOB symbol).
  // 2. During the ouput generation, count[c] is the total number
  //    of characters smaller than c (including the single EOB symbol).
  std::vector<uint32> count(257, 0);

  // count EOB
  bwt_rank_low24[eob_position] = 0;
  count[0] = 1;
  // count other characters
  for (uint32 position = 0; position < bwt_size; ++position) {
    if (position != eob_position) {
      uint32 ch = static_cast<unsigned char>(bwt[position]);
      uint32 rank = count[ch + 1];
      uint32 rank_low24 = rank & 0x00FFFFFF;
      //uint32 rank_hi8 = rank >> 24;
      bwt_rank_low24[position] = (ch << 24) + rank_low24;
      if (0 == rank_low24) {
        rank_milestone_buffer[ch].push_back(position);
      }
      ++count[ch + 1];
    }
  }
  std::partial_sum(count.begin(), count.end(), count.begin());
  assert(count[256] == bwt_size);

  for (int ch = 0; ch < 256; ++ch) {
    rank_milestone_buffer[ch].push_back(bwt_size);
  }

  OutStreamBuffer buffer;
  buffer.Connect(output);
  uint32 position = 0;
  uint32 string_length = 0;
  while (position != eob_position) {
    uint32 bwt_and_rank = bwt_rank_low24[position];
    uint32 ch = bwt_and_rank >> 24;
    buffer.WriteByte(ch);
    ++string_length;
    uint32 rank_low24 = bwt_and_rank & 0x00FFFFFF;
    uint32 rank_hi8 = 0;
    while (rank_milestone_buffer[ch][rank_hi8 + 1] <= position) ++rank_hi8;
    position = count[ch] + (rank_hi8 << 24) + rank_low24;
  }
  // If the BWT or the EOB position contain errors (or are garbage),
  // it is likely that the cycle ends prematurely.
  // In this case, we return false.
  buffer.Disconnect();
  return (string_length + 1 == bwt_size);
}

}  // namespace

bool InverseBWTransformer::Options::Set(const std::string& options_string) {
  if (options_string.length() < 1) return false;
  char algorithm = options_string[0];
  AlgorithmSpecificOptions* options;
  switch (algorithm) {
    case 's':
      options = new SimpleInverseBWTransformer::Options;
      break;
    case 'f':
      options = new FastInverseBWTransformer::Options;
      break;
    default:
      return false;
  }
  std::string algorithm_specific_options_string(options_string, 1);
  bool options_are_valid = options->Set(algorithm_specific_options_string);
  if (options_are_valid) {
    algorithm_id_ = algorithm;
    if (algorithm_specific_options_) delete algorithm_specific_options_;
    algorithm_specific_options_ = options;
  } else {
    delete options;
  }
  return options_are_valid;
}

bool InverseBWTransformer::Transform(
    const char* bwt_block, size_t bwt_block_size,
    size_t eob_position, OutStream* output,
    Options options,
    int64 memory_budget)
{
  assert(bwt_block_size <= options.MaxBlockSize(memory_budget));
  InverseBWTransformer* algorithm = options.GetAlgorithm(memory_budget);
  bool success = algorithm->DoTransform(bwt_block, bwt_block_size,
                                        eob_position,
                                        output);
  delete algorithm;
  return success;
}

}  // namespace dcsbwt
