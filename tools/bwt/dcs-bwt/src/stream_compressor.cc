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

#include "bbb_compress.h"
#include "rl_compress.h"
#include "stream_compressor.h"

#include <string>

namespace dcsbwt {

// Trivial implementations of StreamCompressor and StreamDecompressor
// They serve as an example, and can be useful as a replacement of a real
// (de)compressor in testing, debugging and development.
class TrivialCompressor : public StreamCompressor {
 public:
  class Options : public StreamCompressor::OptionsBase {
   public:
    bool Set(const std::string& options) { return 0 == options.length(); }
    std::string Get() const { return std::string(); }
    int64 SizeInBytes() const { return 1; }
    StreamCompressor* GetCompressor() { return new TrivialCompressor; }
 };

  TrivialCompressor() {}
  virtual ~TrivialCompressor() {}

  virtual void WriteBegin() {}
  virtual void WriteEndPrivate() {}
  virtual void Write(const char* bytes, size_t n) { Emit(bytes, n); }
 private:
  TrivialCompressor(const TrivialCompressor&);
  TrivialCompressor& operator=(const TrivialCompressor&);
};

class TrivialDecompressor : public StreamDecompressor {
 public:
  static StreamDecompressor* Create(const std::string&) {
    return new TrivialDecompressor;
  }
  TrivialDecompressor() {}
  virtual ~TrivialDecompressor() {}
  virtual void ReadBegin() {}
  virtual void ReadEnd() {}
  virtual void Read(char* bytes, size_t n) { GetCompressed(bytes, n); }
  virtual int64 SizeInBytes() const { return 1; }
 private:
  TrivialDecompressor(const TrivialDecompressor&);
  TrivialDecompressor& operator=(const TrivialDecompressor&);
};


/////////////////////// StreamCompressor ///////////////////////

const int StreamCompressor::kBufferSize = 1 << 16;

// Select the compression algorithm based on the first character
// of the options_string and let the algorithm specific options
// processor handle the rest of options_string.
bool StreamCompressor::Options::Set(const std::string& options_string) {
  if (options_string.length() < 1) return false;
  char compression_algorithm = options_string[0];
  OptionsBase* options;
  switch (compression_algorithm) {
    case 't':  // TrivialCompressor (see above)
      options = new TrivialCompressor::Options;
      break;
    case 'b':  // BbbCompressor (see bbb_compressor.h)
      options = new BbbCompressor::Options;
      break;
    case 'c':
      options = new ByteCompressor::Options;
      break;
    case 'r':
      options = new RunLengthCompressor::Options;
      break;
    default:
      return false;
  }
  std::string algorithm_specific_options_string(options_string, 1);
  bool success = options->Set(algorithm_specific_options_string);
  if (success) {
    compression_algorithm_ = compression_algorithm;
    if (algorithm_specific_options_) delete algorithm_specific_options_;
    algorithm_specific_options_ = options;
  } else {
    delete options;
  }
  return success;
}

StreamCompressor* StreamCompressor::Options::GetCompressor() {
  return algorithm_specific_options_->GetCompressor();
}

StreamCompressor* StreamCompressor::Connect(OutStream* output,
                                            Options options) {
  StreamCompressor* compressor = options.GetCompressor();
  if (compressor) {
    // '\n' is used as the terminator of the compression format string.
    std::string compression_format = options.Get() + '\n';
    output->Write(compression_format.data(), compression_format.length());
    compressor->output_.Connect(output);
    compressor->ConnectPrivate(&compressor->output_);
  }
  return compressor;
}

OutStream* StreamCompressor::Disconnect() {
  DisconnectPrivate();
  OutStream* oldoutput = output_.Disconnect();
  delete this;
  return oldoutput;
}

/////////////////////// StreamDecompressor ///////////////////////

StreamDecompressor* StreamDecompressor::Connect(InStreamBuffer* input) {
  // The compression format is identified by the input prefix ending
  // with the first '\n'.
  std::string format;
  char ch = input->ReadByte();
  while (ch != '\n') {
    if (format.length() > 1000) return NULL; // guard against infinite loop
    format += ch;
    ch = input->ReadByte();
  }
  if (0 == format.length()) return NULL;
  char compression_type = format[0];
  std::string options(format, 1);
  StreamDecompressor* decompressor;
  switch (compression_type) {
    case 't':
      decompressor = TrivialDecompressor::Create(options);
      break;
    case 'b':
      decompressor = BbbDecompressor::Create(options);
      break;
    case 'c':
      decompressor = ByteDecompressor::Create(options);
      break;
    case 'r':
      decompressor = RunLengthDecompressor::Create(options);
      break;
    default:
      decompressor = NULL;
  }
  if (decompressor) {
    decompressor->format_ = format;
    decompressor->input_.Connect(input);
    decompressor->ConnectPrivate(decompressor->input_.GetServant());
  }
  return decompressor;
}

InStreamBuffer* StreamDecompressor::Disconnect() {
  DisconnectPrivate();
  InStreamBuffer* oldinput = input_.Disconnect();
  delete this;
  return oldinput;
}

}  // namespace dcsbwt
