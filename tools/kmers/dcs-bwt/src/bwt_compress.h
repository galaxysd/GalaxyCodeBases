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

// Compressor and decompressor based on Burrows-Wheeler transform.
//
// BASIC FUNCTIONALITY
//
// Compresses an arbitrary sequences of bytes.
// The input is divided into blocks, the Burrows-Wheeler transform
// is applied to each block, and the transformed blocks are compressed
// with a stream compressor optimized for compressing BW-transforms.
// See bw_transform.h for more information on the BW-transform and
// stream_compressor.h for more information on the stream compressor.
//
// The simplest way to use the compressor and the decompressor
// are the static member functions BwtCompressor::Compress and
// BwtDecompressor::Decompress. They should be sufficient for
// compressing and decompressing files, at least.
//
// There are also APIs with a finer control over the (de)compression process.
// Here are some reasons for using the low level APIs instead of
// Compress or Decompress.
// * Compress must know the total size of the uncompressed data in advance.
//   The low level API needs to know only one block size at a time.
// * Compress makes its own copy of each block of uncompressed data.
//   For large block sizes, this may waste a significant amount of space.
//   The low level API uses the caller's copy of the uncompressed data.
// * Decompress cannot tell the size of the uncompressed data in advance.
//   The outstream servant receiving the uncompressed data from Decompress
//   must be able to accept all the data in pieces of any size.
//   In the low level API, compression is done one block at a time and
//   the size of each uncompressed block is known before decompressing it.
//
// Both compression APIs have the same output format and
// can be decompressed with either decompression API.
//
// COMPRESSION FORMAT
//
// <compressed data>    := <format identifier><block><block>...<empty block>
// <format identifier>  := BWT<compression format>\n
// <compression format> := see stream_compressor.h
// <block>              := <block size(>0)><compressed block><EOB position>
// <empty block>        := <block size(==0)>
// <block size>         := uncompressed block size encoded with WriteUint63
// <compressed block>   := a block of uncompressed data processed thus:
//                         1. Do Burrows-Wheeler transform as described
//                            in bw_transform.h
//                         2. Compress using <compression format>
// <EOB position>       := an integer encoded with WriteUint63,
//                         see bw_transform.h for an explanation.

#ifndef DCSBWT_BWT_COMPRESS_H__
#define DCSBWT_BWT_COMPRESS_H__

#include <string>

#include "bw_transform.h"
#include "inverse_bwt.h"
#include "stream_compressor.h"
#include "stream.h"
#include "inttypes.h"

namespace dcsbwt {

///////////////////////////////////////////////////////////////////////
// BwtCompressor is a data compressor based on the Burrow-Wheeler transform.
//
// Usage:
//
// BwtCompressor::Compress(...);
//
// or:
//
// BwtCompressor compressor;
// compressor.Connect(...);
// compressor.CompressBlock(...);
// compressor.CompressBlock(...);
// ...
// compressor.Disconnect();
//////////////////////////////////////////////////////////////////////
class BwtCompressor {
 public:

  // Options is a class representing the options of BwtCompressor.
  // Some BwtCompressor member functions take an object of type Options
  // as an argument.
  //
  // Three kinds of options are supported:
  // 1. Options controlling the Burrows-Wheeler transform.
  //    These options have no direct effect on the output of the compressor
  //    but they affect the time and space complexity. Under a limited
  //    memory budget, the space complexity can affect the output block sizes.
  // 2. Options controlling the compression stage after the BW-transform.
  //    These do affect the output and they are encoded in the beginning
  //    of the output in <compression format> (see above).
  // 3. Memory budget controls the amount of main memory used
  //    by the compressor. May affect the choice of the algorithm
  //    and parameters for the BW-transform.
  //
  // The default constructor sets the default values.
  // These can be changed with the SetXxx...() member functions.
  // If the argument to SetXxx...() is invalid, it returns false
  // and does not change any options. Thus intances of Options
  // should always be valid.
  class Options {
   public:
    Options() : memory_budget_(kDefaultMemoryBudget),
                use_printable_characters_(false) {}
    ~Options() {}
    bool SetTransformOptions(std::string options) {
      return transform_options_.Set(options);
    }
    const BWTransformer::Options& GetTransformOptions() const {
      return transform_options_;
    }
    bool SetCompressionOptions(std::string options) {
      bool ok = compression_options_.Set(options);
      // If compression type is 't' (trivial == no compression),
      // use only printable characters for metadata.
      if (ok) use_printable_characters_ = ('t' == options[0]);
      return ok;
    }
    const StreamCompressor::Options& GetCompressionOptions() const {
      return compression_options_;
    }
    bool UsePrintableCharacters() const { return use_printable_characters_; }
    bool SetMemoryBudget(int64 memory_budget) {
      if (memory_budget < kMinMemoryBudget
          || memory_budget > kMaxMemoryBudget) return false;
      memory_budget_ = memory_budget;
      return true;
    }
    int64 GetMemoryBudget() const { return memory_budget_; }

    // The largest block size that can be handled under the current
    // options and memory budget.
    int64 MaxBlockSize() const {
      int64 memory_budget = GetMemoryBudget();
      memory_budget -= kMemoryOverhead;
      memory_budget -= GetCompressionOptions().SizeInBytes();
      return GetTransformOptions().MaxBlockSize(memory_budget);
    }

    // The "optimal" block size under the current options and memory budget.
    // This represents some kind of a sweet spot on the tradeoff curve
    // between block size and compression speed.
    int64 SuggestedBlockSize() const {
      int64 memory_budget = GetMemoryBudget();
      memory_budget -= kMemoryOverhead;
      memory_budget -= GetCompressionOptions().SizeInBytes();
      return GetTransformOptions().SuggestedBlockSize(memory_budget);
    }
   private:
    // Memory overhead by top level code, stack, and minor data structures
    // This does not include the memory for the transformer and the compressor.
    static const int64 kMemoryOverhead = (1LL << 20);
    // kMinMemoryBudget should be high enough for compressing at least
    // a small block under any option configuration.
    static const int64 kMinMemoryBudget = (1LL << 24);
    static const int64 kDefaultMemoryBudget = (1LL << 30);
    static const int64 kMaxMemoryBudget = (3LL << 30);

    int64 memory_budget_;
    BWTransformer::Options transform_options_;
    StreamCompressor::Options compression_options_;
    bool use_printable_characters_;
  };  // class Options

  BwtCompressor() : compressor_(NULL) {}
  ~BwtCompressor() { assert(compressor_ == NULL); }

  // Perform the whole compression in one call.
  // Reads input_size bytes from input and writes the result to output.
  // The input is divided into blocks as follows:
  // * If 0 < max_block_size <= options.MaxBlockSize(), all blocks
  //   except possibly the last one have size max_block_size
  // * If max_block_size == 0, options.SuggestedBlockSize() is used
  //   as the upper bound for the block sizes.
  // * Other values of max_block_size are invalid.
  static void Compress(InStream* input, int64 input_size, OutStream* output,
                       Options options,
                       int64 max_block_size);

  // Start compression, with the compressed data written to output.
  // Once called, can be called again only after Disconnect().
  void Connect(OutStream* output, Options options);

  // End compression.
  // Must be called (exactly once) after each Connect().
  // NOTE: BwtCompressor destructor does not call Disconnect() but fails
  // if the compression was not ended properly with Disconnect().
  // Returns the OutStream that was given as an argument to Connect().
  OutStream* Disconnect();

  // Compress one block.
  // The data to be compressed is in [block, block+block_size).
  //
  // Can be called only between Connect() and Disconnect().
  //
  // Returns true on success.
  // Returns false if block_size is too large for the memory budget.
  // In the case of a too large block size, nothing is written to output
  // and the state of the compressor does not change, i.e.,
  // the compression can be attempted again (using a smaller block size).
  bool CompressBlock(char* block, size_t block_size);

 private:
  StreamMaster<OutStream> rawout_;
  StreamCompressor* compressor_;
  Options options_;

  void WriteEmptyBlock();
  void WriteUint63(int64 value);
  void WriteUint63Printable(int64 value);

  BwtCompressor(const BwtCompressor&);
  BwtCompressor& operator=(const BwtCompressor&);
};

///////////////////////////////////////////////////////////////////
// BwtDecompressor decompresses data compressed by BwtCompress
//
// Basic usage:
//
// BwtDecompressor::Decompress(...);
//
// or:
//
// BwtDecompressor decompressor;
// decompressor.Connect(...);
// while ( ! decompressor.IsFinished()) {
//   int64 size = decompressor.GetUncompressedSizeOfNextBlock();
//   ...
//   decompressor.DecompressBlock(...);
// }
// decompressor.Disconnect();
//////////////////////////////////////////////////////////////////
class BwtDecompressor {
 public:

  // Options is a class representing the options of BwtDecompressor.
  // Some BwtDecompressor member functions take an object of type Options
  // as an (optional) argument.
  //
  // Two kinds of options are supported:
  // 1. Options controlling the inverse BW transform.
  // 2. Memory budget.
  // There are no options controlling the decompression stage preceding
  // the inverse BW transform. Those options are read from the compressed data.
  class Options {
   public:
    Options() : memory_budget_(kDefaultMemoryBudget) {}
    ~Options() {}
    bool SetTransformOptions(std::string options) {
      return transform_options_.Set(options);
    }
    const InverseBWTransformer::Options& GetTransformOptions() const {
      return transform_options_;
    }
    bool SetMemoryBudget(int64 memory_budget) {
      if (memory_budget < kMinMemoryBudget
          || memory_budget > kMaxMemoryBudget) return false;
      memory_budget_ = memory_budget;
      return true;
    }
    int64 GetMemoryBudget() const { return memory_budget_; }

    // The largest block size that can be handled under the current
    // options and memory budget.
    int64 MaxBlockSize(int64 decompressor_size) const {
      int64 memory_budget = GetMemoryBudget() - decompressor_size;
      return GetTransformOptions().MaxBlockSize(memory_budget);
    }

   private:
    static const int64 kMinMemoryBudget = (1LL << 24);
    static const int64 kDefaultMemoryBudget = (1LL << 30);
    static const int64 kMaxMemoryBudget = (3LL << 30);

    int64 memory_budget_;
    InverseBWTransformer::Options transform_options_;
  };

  BwtDecompressor() : decompressor_(NULL),
                      uncompressed_block_size_(0),
                      input_error_(0),
                      low_memory_failure_(false) {}
  ~BwtDecompressor() {
    assert(decompressor_ == NULL);
  }

  // Full decompression with one call.
  // Reads compressed data from input and writes uncompressed data to output.
  //
  // Returns 0 on success.
  // Returns a negative value if the input data is invalid.
  // The negative value is the return value of GetInputError().
  // Returns a positive value if the memory budget is too low.
  // The positive value is the uncompressed size of the block
  // that could not be decompressed under the memory budget.
  static int64 Decompress(InStreamBuffer* input, OutStream* output,
                         Options options);

  // Start decompressing data coming from input.
  // Returns true on success.
  // Returns false if the input is invalid.
  bool Connect(InStreamBuffer* input);

  // End decompression.
  // Must be called (exactly once) after each Connect() except
  // must not be called if an input error has occurred at some point.
  // NOTE: BwtCompressor destructor does not call Disconnect() but fails
  // if the compression was not ended properly with Disconnect().
  //
  // Normally called when IsFinished() returns true.
  // Can be called before IsFinished() returns true, but then
  // it leaves input in the middle of the compressed data in
  // a generally unrecoverable state.
  //
  // Returns the InStreamBuffer* that was given as an argument to Connect().
  InStreamBuffer* Disconnect();

  // Returns true when there is data to be decompressed.
  bool HasData() const { return uncompressed_block_size_ > 0; }

  // Returns the uncompressed size of the block about to be decompressed next.
  // Returns 0 when there is no next block (i.e., when decompressor is not
  // connected, all data has been decompressed, or an input error has
  // occurred.)
  int64 GetUncompressedSizeOfNextBlock() const {
    return uncompressed_block_size_;
  }

  // Decompress a single block and write the uncompressed data to output.
  // Can only be called when HasData() is true.
  // Return true on success.
  // Return false if input is invalid or the memory budget is too low.
  //
  // In the case of too low memory budget, no input is consumed, and
  // the state of the decompressor does not change. Thus,
  // the decompression may be attempted again with a higher memory
  // budget or a different algorithm. A second failure on the same
  // block causes a crash, though. (This is to prevent an infinite
  // loop.)
  bool DecompressBlock(OutStream* output, Options options);

  int GetInputError() const { return input_error_; }

 private:
  StreamMaster<InStreamBuffer> rawin_;
  StreamDecompressor* decompressor_;
  int64 uncompressed_block_size_;
  int input_error_;
  bool low_memory_failure_;
  bool use_printable_characters_;

  void SetError(int error);
  int64 ReadUint63();
  int64 ReadUint63Printable();
};

}  // namespace dcsbwt

#endif  // DCSBWT_BWT_COMPRESS_H__
