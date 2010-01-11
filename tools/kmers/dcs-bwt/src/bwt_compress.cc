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

#include <iostream>
#include <algorithm>  // for min
#include <cstdio>

#include "bwt_compress.h"
#include "prefix_doubling.h"
//#include "inverse_bwt.h"
//#include "bbb_compress.h"
#include "stream.h"

namespace dcsbwt {

extern int verbosity;

namespace {
// Integer division that rounds up
// REQUIRE: both arguments are positive
uint64 DivideAndRoundUp(uint64 numerator, uint64 denominator) {
  assert(numerator > 0);
  assert(denominator > 0);
  return 1 + (numerator - 1) / denominator;
}
}  // namespace

//const int64 BwtCompressor::kMinMemoryBudget = (1LL << 24);
//const int64 BwtCompressor::kDefaultMemoryBudget = (1LL << 30);
//const int64 BwtCompressor::kMaxMemoryBudget = (2LL << 30);

void BwtCompressor::Compress(InStream* input, int64 input_size,
                             OutStream* output, Options options,
                             int64 max_block_size)
{
  assert(input_size >= 0);
  assert(max_block_size >= 0);
  assert(max_block_size <= options.MaxBlockSize());
  // Determine block sizes
  // There are (at most) two sizes of blocks: small and large
  int64 num_large_blocks, large_block_size;
  int64 num_small_blocks, small_block_size;
  if (0 == input_size) {
    // Special case for empty input
    num_large_blocks = num_small_blocks = 0;
    large_block_size = small_block_size = 0;
  } else if (max_block_size > 0) {
    // All blocks except possibly the last one have size max_block_size.
    large_block_size = std::min(max_block_size, input_size);
    num_large_blocks = input_size / large_block_size;
    small_block_size = input_size - num_large_blocks * large_block_size;
    num_small_blocks = (small_block_size > 0) ? 1 : 0;
  } else {
    // max_block_size == 0: use SuggestedBlockSize()
    max_block_size = options.SuggestedBlockSize();
    // Minimize the number of blocks and make block sizes as even as possible
    int64 num_blocks = DivideAndRoundUp(input_size, max_block_size);
    large_block_size = DivideAndRoundUp(input_size, num_blocks);
    small_block_size = large_block_size - 1;
    num_small_blocks = num_blocks * large_block_size - input_size;
    num_large_blocks = num_blocks - num_small_blocks;
  }

  assert(num_large_blocks >= 0);
  assert(num_small_blocks >= 0);
  assert(input_size == 0 || num_large_blocks > 0);
  assert(0 <= small_block_size);
  assert(small_block_size <= large_block_size);
  assert(large_block_size <= max_block_size);
  assert(large_block_size <= options.MaxBlockSize());
  assert(num_small_blocks * small_block_size
           + num_large_blocks * large_block_size == input_size);

  if (verbosity > 0) {
    std::clog << "Compressing " << input_size << " bytes using block size "
              << large_block_size << std::endl;
  }

  // Compress
  BwtCompressor compressor;
  compressor.Connect(output, options);
  std::vector<char> block(large_block_size, 0);
  char* block_begin = &block[0];
  while (num_large_blocks--) {
    input->Read(block_begin, large_block_size);
    // CompressBlock should not fail for this block size
    bool ok = compressor.CompressBlock(block_begin, large_block_size);
    assert(ok);
  }
  std::fill(block_begin + small_block_size, block_begin + large_block_size, 0);
  while (num_small_blocks--) {
    input->Read(block_begin, small_block_size);
    bool ok = compressor.CompressBlock(block_begin, small_block_size);
    assert(ok);
  }
  compressor.Disconnect();
}

void BwtCompressor::Connect(OutStream* output, Options options) {
  assert(compressor_ == NULL);
  rawout_.Connect(output);
  rawout_->Write("BWT", 3);
  compressor_ = StreamCompressor::Connect(output,
                                          options.GetCompressionOptions());
  assert(compressor_ != NULL);
  options_ = options;
}

OutStream* BwtCompressor::Disconnect() {
  assert(compressor_ != NULL);
  compressor_->Disconnect(); compressor_ = NULL;
  WriteEmptyBlock();
  return rawout_.Disconnect();
}

bool BwtCompressor::CompressBlock(char* block, size_t block_size) {
  assert(compressor_ != NULL);
  assert(block_size > 0);
  if (block_size > options_.MaxBlockSize()) { return false; }

  if (verbosity > 2) {
    std::clog << "Compressing block of size " << block_size << std::endl;
  }

  WriteUint63(block_size);
  compressor_->WriteBegin();
  uint64 eob_position =
      BWTransformer::Transform(block, block_size,
                               compressor_,
                               options_.GetTransformOptions(),
                               options_.GetMemoryBudget());
  compressor_->WriteEnd();
  WriteUint63(eob_position);
  return true;
}

void BwtCompressor::WriteEmptyBlock() { WriteUint63(0); }

void BwtCompressor::WriteUint63(int64 value) {
  if (options_.UsePrintableCharacters()) return WriteUint63Printable(value);
  assert(value >= 0);
  char bytes[9];
  int position = 9;
  unsigned char byte = value & 127;
  value >>= 7;
  bytes[--position] = byte;
  while (value > 0) {
    byte = (value & 127) | 128;
    value >>= 7;
    assert(position > 0);
    bytes[--position] = byte;
  }
  rawout_->Write(bytes + position, 9 - position);
}

void BwtCompressor::WriteUint63Printable(int64 value) {
  // Write value in decimal ended by newline.
  assert(value >= 0);
  char bytes[12];
  int length = sprintf(bytes, "%lld", value);
  bytes[length] = '\n';  // Replace '\0' with '\n'
  rawout_->Write(bytes, length + 1);
}


///////////////////////////// BwtDecompressor ///////////////////////////////

int64 BwtDecompressor::Decompress(InStreamBuffer* input, OutStream* output,
                                 Options options) {
  BwtDecompressor decompressor;
  bool success = decompressor.Connect(input);
  if (!success) return decompressor.GetInputError();
  while (decompressor.HasData()) {
    success = decompressor.DecompressBlock(output, options);
    if (!success) {
      if (decompressor.GetInputError()) {
        return decompressor.GetInputError();
      } else {
        int64 too_large_block_size =
            decompressor.GetUncompressedSizeOfNextBlock();
        decompressor.Disconnect();
        return too_large_block_size;
      }
    }
  }
  decompressor.Disconnect();
  return 0;
}

bool BwtDecompressor::Connect(InStreamBuffer* input) {
  assert(!low_memory_failure_);
  input_error_ = 0;  // clear previous errors
  rawin_.Connect(input);
  // Read the initial "BWT"
  char magic_number[3];
  rawin_->Read(magic_number, 3);
  if (magic_number[0] != 'B' || magic_number[1] != 'W'
      || magic_number[2] != 'T') {
    SetError(-1);
    return false;
  }
  assert(decompressor_ == NULL);
  decompressor_ = StreamDecompressor::Connect(input);
  if (!decompressor_) {
    SetError(-1);
    return false;
  }
  use_printable_characters_ = ('t' == (decompressor_->GetFormat())[0]);
  uncompressed_block_size_ = ReadUint63();
  if (uncompressed_block_size_ < 0) {
    SetError(-1);
    return false;
  }
  return true;
}

InStreamBuffer* BwtDecompressor::Disconnect() {
  assert(decompressor_ != NULL);
  decompressor_->Disconnect(); decompressor_ = NULL;
  uncompressed_block_size_ = 0;
  low_memory_failure_ = false;
  return rawin_.Disconnect();
}

bool BwtDecompressor::DecompressBlock(OutStream* output, Options options) {
  assert(uncompressed_block_size_ > 0);
  assert(decompressor_ != NULL);
  int64 max_block_size = options.MaxBlockSize(decompressor_->SizeInBytes());
  if (uncompressed_block_size_ > max_block_size) {
    assert(!low_memory_failure_);
    low_memory_failure_ = true;
    return false;
  } else {
    low_memory_failure_ = false;
  }
  std::vector<char> bwt(uncompressed_block_size_ + 1);
  decompressor_->ReadBegin();
  if (verbosity > 1) {
    std::clog << "Decompressing a block of size " <<  uncompressed_block_size_
              << std::endl;
  }
  decompressor_->Read(&bwt[0], bwt.size());
  decompressor_->ReadEnd();
  int64 eob_position = ReadUint63();
  if (eob_position < 0 || eob_position > uncompressed_block_size_) {
    SetError(-1);
    return false;
  }
  if ( ! InverseBWTransformer::Transform(&bwt[0], bwt.size(),
                                         eob_position, output,
                                         options.GetTransformOptions(),
                                         options.GetMemoryBudget())) {
    SetError(-1);
    return false;
  }

  // Get ready for the next block.
  uncompressed_block_size_ = ReadUint63();
  if (uncompressed_block_size_ < 0) {
    SetError(-1);
    return false;
  }
  return true;
}

// Read an integer value in the range [0,2^63) from the input
int64 BwtDecompressor::ReadUint63() {
  if (use_printable_characters_) return ReadUint63Printable();
  unsigned char byte;
  int64 value = 0;
  int num_bytes = 0;  // used only for checking
  do {
    assert(++num_bytes <= 9);
    byte = rawin_->ReadByte();
    value = (value << 7) + (byte & 127);
  } while (byte > 127);
  return value;
}

// Read an integer value in the range [0,2^63) from the input
int64 BwtDecompressor::ReadUint63Printable() {
  unsigned char byte;
  char bytes[32];
  int num_bytes = 0;  // used only for checking
  do {
    assert(num_bytes < 32);
    byte = rawin_->ReadByte();
    bytes[num_bytes++] = byte;
  } while ('\n' != byte);
  char* next;
  int64 value = strtoll(bytes, &next, 10);
  return value;
}

// Store the error and clean up
void BwtDecompressor::SetError(int error) {
  assert(error < 0);
  input_error_ = error;
  uncompressed_block_size_ = 0;
  low_memory_failure_ = false;
  if (decompressor_) {
    decompressor_->Disconnect(); decompressor_ = NULL;
  }
  rawin_.Disconnect();
}

}  // namespace dcsbwt
