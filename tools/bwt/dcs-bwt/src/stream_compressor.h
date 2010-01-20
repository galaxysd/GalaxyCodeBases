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

// Base classes for stream compressors and decompressors.
// A stream (de)compressor can handle an arbitrary stream
// of data in a nearly online fashion.
// The streams are based on the stream protocol described in stream.h.
//
// BASIC USAGE
//
// * Initialize compressor:              * Initialize decompressor:
// StreamCompressor* compressor          StreamDecompressor* decompressor
//  = StreamCompressor::Connect(output);  = StreamDecompressor::Connect(input);
//
// * Start compressing block:            * Start decompressing block:
// compressor->WriteBegin();             decompressor->ReadBegin();
//
// * Compress block:                     * Decompress block:
// compressor->Write(...);               decompressor->Read(...);
// compressor->Write(...);               decompressor->Read(...);
// ...                                   ...
//
// * Finish compressing block:           * Finish decompressing block
// compressor->WriteEnd();               decompressor->ReadEnd();
//
// * Compress more blocks:               * Decompress more blocks
// compressor->WriteBegin();             decompressor->ReadBegin();
// ...                                   ...
// compressor->WriteEnd();               decompressor->ReadEnd();
// ...                                   ...
//
// * End compression:                    * End decompression:
// compressor->Disconnect();             decompressor->Disconnect();
//
// SYNCHRONIZATION
//
// Synchronization between compressor and decompressor happens
// at block boundaries:
// * The sizes of individual writes and reads do not have to match.
// * The client must ensure that at the end of the block (i.e., just
//   before calling WriteEnd/ReadEnd), the compressor and
//   decompressor are synchronized on the uncompressed data, i.e.,
//   that the compressor and the decompressor have processed
//   the same amount of uncompressed data.
// * WriteEnd/ReadEnd then ensures that the compressor and decompressor
//   are synchronized on the compressed data, too. Then it is possible
//   to write some data (such as the size of the next block) to the
//   output past the compressor as long as the same amount of data
//   is read past the decompressor.
//
// OPTIONS
//
// Compression options can be provided like this:
//
// string option_string = ...;
// StreamCompressor::Options options;
// if ( ! options.Set(option_string)) { ... // error: invalid option_string }
// StreamCompressor* compressor = StreamCompressor::Connect(output, options);
// ...
//
// The compressor writes the options to the beginning of the compressed data
// and the decompressor reads them from there (both during Connect()).
//
// The first character of option_string specifies the compression algorithm;
// the rest provides options specific to the algorithm.
// See stream_compressor.cc for available algorithms and algorithm-specific
// files for the rest of the options.
//
// IMPLEMENTING (DE)COMPRESSORS
//
// New stream compression/decompression algorithms can be defined
// by deriving from StreamCompressor/StreamDecompressor.
// The simple TrivialCompressor/TrivialDecompressor
// in stream_compressor.cc is a useful example.
//
// The set of compression algorithms is hard-coded in the definition
// of the functions StreamCompressor::Options::Set() and
// StreamDecompressor::Connect() (in stream_compressor.cc).
// Those are the only places that need to be modified when adding a new
// compression algorithm. In particular, there is no need to modify
// this file stream_compressor.h.
//
// Changes to an already supported algorithm including changing
// the algorithm-specific options should need no modification of either
// stream_compressor.h or stream_compressor.cc.
//
#ifndef DCSBWT_STREAM_COMPRESSOR_H__
#define DCSBWT_STREAM_COMPRESSOR_H__

#include <string>
#include <cstddef>  // for size_t

#include "stream.h"

namespace dcsbwt {

// StreamCompressor is the base class for stream compressors.
//
// A derived class should:
// 1. define the pure virtual functions Write, WriteBegin and WriteEndPrivate
// 2. have a nested class derived from StreamCompressor::OptionsBase
// 3. add an entry to StreamCompressor::Options::Set in stream_compressor.cc
// 4. write compressed output using Emit and EmitByte (which do buffering)
//
class StreamCompressor : public OutStream {
 protected:
  class OptionsBase;

 public:
  static const int kBufferSize;

  // A class representing all options affecting compression and decompression.
  class Options {
   public:
    Options() : algorithm_specific_options_(NULL) { Set("t"); }

    // Set/Get serialization mechanism is the simplest way
    // to implement copying.
    Options(const Options& other) : algorithm_specific_options_(NULL)
    { Set(other.Get()); }
    const Options& operator=(const Options& other) {
      Set(other.Get());
      return *this;
    }
    ~Options() {
      if (algorithm_specific_options_)
        delete algorithm_specific_options_;
    }

    // Initialize options from a string.
    // Returns false in case of an invalid options_string
    bool Set(const std::string& options_string);

    // Get() returns a string that:
    // 1. as an argument to Set() sets any options object to this object's
    //    current state.
    // 2. when written to the compressed stream and read by
    //    StreamDecompressor, sets the decompressor in the the equivalent
    //    state.
    // 3. does not contain '\n' (which is used as a terminator in the stream)
    std::string Get() const {
      if (algorithm_specific_options_) {
        return compression_algorithm_ + algorithm_specific_options_->Get();
      } else {
        return std::string();
      }
    }
    int64 SizeInBytes() const {
      return kBufferSize + algorithm_specific_options_->SizeInBytes();
    }
   private:
    friend class StreamCompressor;
    StreamCompressor* GetCompressor();

    char compression_algorithm_;
    OptionsBase* algorithm_specific_options_;
  };

  StreamCompressor() : output_(kBufferSize) {}
  virtual ~StreamCompressor() {}

  // Create a compressor with the compressed output going to output.
  // The compression format specified by options is written to output.
  static StreamCompressor* Connect(OutStream* output,
                                   Options options);
  // Disconnect deletes this
  OutStream* Disconnect();

  virtual void WriteBegin() =0;

  void WriteEnd() { WriteEndPrivate(); output_.Flush(); }

  // inherits from OutStream:
  // virtual void Write() =0;

 protected:

  // A compressor class derived from StreamCompressor should
  // have a nested class derived from OptionsBase.
  class OptionsBase {
   public:
    OptionsBase() {}
    virtual ~OptionsBase() {}

    // Set/Get should behave as StreamCompressor::Options::Set/Get
    // The algorithm identifying initial character is omitted from
    // the input to Set and should be omitted from the output of Get.
    virtual bool Set(const std::string& options_string) =0;
    virtual std::string Get() const =0;
    virtual int64 SizeInBytes() const =0;

    // Returns a compressor object corresponding to the current options.
    virtual StreamCompressor* GetCompressor() =0;
  };

  // All compressed output should be written using Emit and EmitByte.
  inline void Emit(const char* bytes, size_t n) { output_.Write(bytes, n); }
  inline void EmitByte(unsigned char byte) {
    output_.WriteByte(byte);
  }

 private:
  OutStreamBuffer output_;

  virtual void WriteEndPrivate() =0;
  virtual void ConnectPrivate(OutStreamBuffer* out) {};
  virtual void DisconnectPrivate() {};

  StreamCompressor(const StreamCompressor&);
  StreamCompressor& operator=(const StreamCompressor&);
};

////////////////////////////////////////////////////////////////////
// StreamDecompressor is the base class for stream decompressors
//
// A derived class should:
// 1. define the pure virtual functions Read, ReadBegin and ReadEnd
// 2. define a static factory member function Create
// 3. add a call to Create in StreamDecompressor::Connect
// 4. read compressed data using GetCompressed and GetCompressedByte
////////////////////////////////////////////////////////////////////
class StreamDecompressor : public InStream {
 public:
  // Factory method for creating the decompressor.
  // The decompression algorithm and its options are read from input.
  // Returns NULL if (the beginning of) the input is invalid.
  static StreamDecompressor* Connect(InStreamBuffer* input);

  StreamDecompressor() { }
  virtual ~StreamDecompressor() { }

  // Disconnect deletes this.
  InStreamBuffer* Disconnect();

  // Returns the string read from input at creation.
  std::string GetFormat() const { return format_; }

  virtual void ReadBegin() =0;
  virtual void ReadEnd() =0;
  // inherits from InStream:
  // virtual void Read(char* bytes, size_t n) =0;

  virtual int64 SizeInBytes() const =0;

 protected:

  // All compressed data should be read using GetCompressed
  // and GetCompressedByte.
  inline void GetCompressed(char* bytes, size_t n) {
    input_->ReadFast(bytes, n);
  }
  inline unsigned char GetCompressedByte() {
    unsigned char byte = input_->ReadByte();
    return byte;
  }

 private:
  StreamMaster<InStreamBuffer> input_;
  std::string format_;

  virtual void ConnectPrivate(InStreamBuffer* out) {};
  virtual void DisconnectPrivate() {};

  StreamDecompressor(const StreamDecompressor&);
  StreamDecompressor& operator=(const StreamDecompressor&);
};

}  // namespace dcsbwt

#endif  // DCSBWT_STREAM_COMPRESSOR_H__
