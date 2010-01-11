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

// An abstract interface for byte streams defined by two
// base classes OutStream and InStream
// useful for building byte stream processing pipelines
// with interchangeable components.
// There are also derived classes providing buffering and file access.
// TODO: Add more derived classes: string/array access etc.
//
// BASIC PROTOCOL
//
// When two parties want to pass a byte stream between them, the supplier
// obviously decides what the content of the stream is but either or
// both parties may want to have a say on other issues:
// - When to move data, how much to move at a time, and when to stop for good?
// - How to perform the actual transfer of data?
//
// In the protocol used here, one party is the master and the other party
// is the servant. The master decides on all the issues mentioned above:
// - The servant is a class derived from either OutStream or
//   InStream depending on the direction of the stream.
// - The master holds a pointer or a reference to the servant and
//   calls the virtual member function Write or Read to move bytes.
// - The data transfer is done using a memory area (a buffer) that
//   the master supplies (as an argument to the call), and the
//   servant performs the actual movement of data.
//
// The protocol provides no way for a servant to report errors such as
// an end of data, an end of capacity, or invalid stream content.
// There are also no methods for setting up or cleaning up.
// A derived class can and should offer such mechanism when necessary
// but to a master that knows only the static type OutStream or InStream,
// there is only Read or Write that never fails.
//
// EXAMPLE
//
// A simple compressor pipeline that reads from infile and
// writes to outfile might be setup up and run like this:
//
// InStreamFromFile instream(infile);
// OutStreamToFile outstream(outfile);
// Compressor compressor;
// compressor.Connect(outstream);
// char buffer[kBufferSize];
// while (...) {
//   instream.Read(buffer, kBufferSize);
//   compressor.Write(buffer, kBufferSize);
// }
// ...
// instream.Read(buffer, num_remaining_bytes);
// compressor.Write(buffer, num_remaining_bytes);
// compressor.Disconnect();
//
// Compressor might be defined like this:
//
// class Compressor : public OutStream {
//  public:
//   void Connect(OutStream* out) { output_.Connect(out); }
//   void Disconnect() { output_.Disconnect(); }
//   virtual void Write(const char* bytes, size_t n) {
//     for (;n; --n) CompressByte(*bytes++);
//   }
//  private:
//   OutStreamBuffer output_;
//   void EmitByte(unsigned char byte) { output_.WriteByte(byte); }
//   void CompressByte(unsigned char byte);
// };

#ifndef DCSBWT_STREAM_H__
#define DCSBWT_STREAM_H__

#include "inttypes.h"

#include <vector>
#include <cstdio>

namespace dcsbwt {

////////////////////////////////////////////////////////////////
// OutStream is the base class for servants on the receiving end
// of a byte stream.
//
// Typical usage:
//
// class MyOutStream : public OutStream {
//   ...
//   public: virtual void Write(const char* bytes, size_t n) { ... }
//   ...
// };
////////////////////////////////////////////////////////////////
class OutStream {
 public:
  OutStream() { }
  virtual ~OutStream() { }

  // Do something with the data in [bytes, bytes+n)
  virtual void Write(const char* bytes, size_t n) =0;

 private:
  OutStream(const OutStream&);
  OutStream& operator=(const OutStream&);
};

////////////////////////////////////////////////////////////////
// InStream is the base class for servants on the supplying end
// of a byte stream.
//
// Usage is analogous to OutStream.
////////////////////////////////////////////////////////////////
class InStream {
 public:
  InStream() { }
  virtual ~InStream() { }

  // Fill [bytes, bytes+n) with data
  virtual void Read(char* begin, size_t n) =0;

 private:
  InStream(const InStream&);
  InStream& operator=(const InStream&);
};

////////////////////////////////////////////////////////////////
// StreamMaster helps the master end of a stream connection in setting up
// and maintaining the connection. In particular, it provides some
// protection against trying to use an uninitialized servant pointer.
//
// Being a template, StreamMaster can be used with InStream and OutStream
// as well as their derivatives. All public members of the servant
// can be accessed through operator->
//
// StreamMaster is typically used as a member of a class rather than
// a base class; see the buffer classes below for examples.
//
// Typical usage:
//
// class C {
//  public:
//   void Connect(OutStream* servant) { master_.Connect(servant); }
//   void foo() { ... master_->Write(...); ... }
//   OutStream* Disconnect() { return master_.Disconnect(); }
//  private:
//   StreamMaster<OutStream> master_;
// };
////////////////////////////////////////////////////////////////
template <typename Servant>
class StreamMaster {
 public:
  typedef Servant ServantType;

  StreamMaster() : servant_(NULL) {}
  ~StreamMaster() { assert(servant_ == NULL); }

  void Connect(Servant* newservant) {
    assert(servant_ == NULL);
    assert(newservant != NULL);
    servant_ = newservant;
  }
  Servant* Disconnect() {
    assert(servant_ != NULL);
    Servant* oldservant = servant_;
    servant_ = NULL;
    return oldservant;
  }
  bool IsConnected() const { return NULL != servant_; }

  Servant* operator->() { return GetServant(); }
  operator Servant* () { return GetServant(); }
  Servant* GetServant() { assert(servant_ != NULL); return servant_; }

 private:
  Servant* servant_;
  StreamMaster(const StreamMaster&);
  StreamMaster& operator=(const StreamMaster&);
};

////////////////////////////////////////////////////////////////
// OutStreamBuffer provides buffering for an outstream master.
// The main purpose is to make small writes faster by avoiding
// a (virtual) function call for every write.
// Writing a single byte is particularly simple and fast.
//
// Typical usage:
//
// class C {
//  public:
//   void Connect(OutStream* servant) { buffer_.Connect(servant); }
//   void foo() { ... buffer_.Write(...); ... }
//   void bar() { ... buffer_.WriteByte(...); ... }
//   OutStream* Disconnect() { return buffer_.Disconnect(); }
//  private:
//   OutStreamBuffer buffer_;
// };
////////////////////////////////////////////////////////////////
class OutStreamBuffer {
 public:
  static const int kDefaultBufferSize = (1 << 14);
  explicit OutStreamBuffer(size_t buffer_size = kDefaultBufferSize)
      : buffer_(buffer_size), next_free_slot_(buffer_.begin()) {}
  ~OutStreamBuffer() { }

  void Connect(OutStream* servant) { master_.Connect(servant); }
  OutStream* Disconnect() { Flush(); return master_.Disconnect(); }
  bool IsConnected() const { return master_.IsConnected(); }

  inline void Write(const char* bytes, size_t n) {
    assert(IsConnected());
    if (n < FreeSpace()) WriteToBuffer(bytes, n);
    else FlushAndWrite(bytes, n);
    assert(FreeSpace() > 0);
  }
  inline void WriteByte(unsigned char byte) {
    assert(IsConnected());
    assert(FreeSpace() > 0);
    *next_free_slot_++ = byte;
    if (FreeSpace() == 0) Flush();
  }
  void Flush();

  // Change the size of the buffer.
  // Any data in the buffer is flushed.
  // Can be used for releasing the space taken by the buffer
  // by giving a small value as an argument.
  void Reset(size_t size = kDefaultBufferSize);

 private:
  StreamMaster<OutStream> master_;
  std::vector<char> buffer_;
  std::vector<char>::iterator next_free_slot_;

  inline size_t FreeSpace() const {
    assert(buffer_.end() - next_free_slot_ >= 0);
    return buffer_.end() - next_free_slot_;
  }
  inline void WriteToBuffer(const char* bytes, size_t n) {
    next_free_slot_ = std::copy(bytes, bytes+n, next_free_slot_);
  }
  void FlushAndWrite(const char* bytes, size_t n);

  OutStreamBuffer(const OutStreamBuffer&);
  OutStreamBuffer& operator=(const OutStreamBuffer&);
};

////////////////////////////////////////////////////////////////
// InStreamBuffer is the InStream counterpart to OutStreamBuffer (see above).
//
// Unlike OutStreamBuffer, InStreamBuffer does not support flushing.
// Flushing would push data from servant to master, which is against
// basic idea of the protocol. More appropriate would be to send
// the unused data back to where it came from, but it would be
// unreasonable to expect every InStream object to be able to
// move data backwards. For example, a decompressor cannot in general
// reverse the decompression (which is not the same as compression).
//
// Thus any unused data is kept in the buffer, and can be read even if
// the buffer is disconnected from a servant or connected to a new servant.
// Only explicit calls to Clear() or Reset() discard the data.
// The amount of unused data can be found with AvailableInBuffer().
//
// InStreamBuffer is a subclass of InStream, so that any master
// can access the data in the buffer.
// NOTE: Read(...) is the virtual function of InStream.
// ReadFast(...) is an inlined non-virtual function that does
// the same thing (but faster).
//
// If losing data in the buffer at the end is not a problem, InStreamBuffer
// can be used internally similarly to OutStreamBuffer:
//
// class C {
//  public:
//   void Connect(InStream* servant) { buffer_.Connect(servant); }
//   void foo() { ... buffer_.ReadFast(...); ... }
//   void bar() { ... buffer_.ReadByte(...); ... }
//   InStream* Disconnect() { return buffer_.Disconnect(); }
//  private:
//   InStreamBuffer buffer_;
// };
//
// If the remaining data should not be lost, one can instead do this:
//
// class C {
//  public:
//   void Connect(InStreamBuffer* buffer) { buffer_.Connect(buffer); }
//   void foo() { ... buffer_->ReadFast(...); ... }
//   void bar() { ... buffer_->ReadByte(...); ... }
//   InStreamBuffer* Disconnect() { return buffer_.Disconnect(); }
//  private:
//   StreamMaster<InStreamBuffer> buffer_;
// };
////////////////////////////////////////////////////////////////
class InStreamBuffer : public InStream {
 public:
  static const int kDefaultBufferSize = (1 << 12);
  explicit InStreamBuffer(size_t buffer_size = kDefaultBufferSize)
      : buffer_(buffer_size), next_unused_byte_(buffer_.end()) {}
  virtual ~InStreamBuffer() { }

  void Connect(InStream* servant) { master_.Connect(servant); }
  InStream* Disconnect() { return master_.Disconnect(); }
  bool IsConnected() const { return master_.IsConnected(); }

  virtual void Read(char* bytes, size_t n) { ReadFast(bytes, n); }
  inline void ReadFast(char* bytes, size_t n) {
    if (n <= AvailableInBuffer()) ReadFromBuffer(bytes, n);
    else ReadAndRefill(bytes, n);
  }
  inline unsigned char ReadByte() {
    if (AvailableInBuffer() == 0) Refill();
    return *next_unused_byte_++;
  }

  void Clear() { next_unused_byte_ = buffer_.end(); }

  // Change the size of the buffer.
  // WARNING: Any data in the buffer is lost.
  // Mainly useful for releasing the space taken by the buffer
  // by giving a small argument.
  void Reset(size_t size = kDefaultBufferSize);

  inline size_t AvailableInBuffer() const {
    assert(buffer_.end() - next_unused_byte_ >= 0);
    return buffer_.end() - next_unused_byte_;
  }

 private:
  StreamMaster<InStream> master_;
  std::vector<char> buffer_;
  std::vector<char>::iterator next_unused_byte_;

  inline void ReadFromBuffer(char* bytes, size_t n) {
    assert(n <= AvailableInBuffer());
    std::copy(next_unused_byte_, next_unused_byte_ + n, bytes);
    next_unused_byte_ += n;
  }
  void ReadAndRefill(char* bytes, size_t n);
  void Refill();

  InStreamBuffer(const InStreamBuffer&);
  InStreamBuffer& operator=(const InStreamBuffer&);
};

////////////////////////////////////////////////////////////////
// OutStreamToFile and InStreamFromFile are used for streaming
// to/from a file.
////////////////////////////////////////////////////////////////
class OutStreamToFile : public OutStream {
 public:
  explicit OutStreamToFile(FILE* file) : file_(file), no_errors_(true) {}
  virtual ~OutStreamToFile() {}

  virtual void Write(const char* bytes, size_t n) {
    if (fwrite(bytes, 1, n, file_) != n) no_errors_ = false;
    /*
      for (; n; --n) {
        unsigned char byte = *bytes++;
        std::clog << "Wrote to file: " << int(byte);
      }
    */
  }

  bool NoErrors() const { return no_errors_; }

 private:
  FILE* file_;
  bool no_errors_;

  OutStreamToFile(const OutStreamToFile&);
  OutStreamToFile& operator=(const OutStreamToFile&);
};

class InStreamFromFile : public InStream {
 public:
  explicit InStreamFromFile(FILE* file) : file_(file), bytes_read_(0) {}
  virtual ~InStreamFromFile() {}

  virtual void Read(char* bytes, size_t n) {
    int64 size = fread(bytes, 1, n, file_);
    /*
    for (; size; --size) {
      unsigned char byte = *bytes++;
      std::clog << "Read from file: " << int(byte);
    }
    */
    bytes_read_ += size;
  }

  // Over-reading due to internal buffers is acceptable behaviour.
  // Thus no error is reported even if reading failed.
  // A client may use BytesReadFromFile() instead to check that the expected
  // number of bytes was actually read from the file.
  int64 BytesReadFromFile() const { return bytes_read_; }

 private:
  FILE* file_;
  int64 bytes_read_;

  InStreamFromFile(const InStreamFromFile&);
  InStreamFromFile& operator=(const InStreamFromFile&);
};

}  // namespace dcsbwt

#endif  // DCSBWT_STREAM_H__
