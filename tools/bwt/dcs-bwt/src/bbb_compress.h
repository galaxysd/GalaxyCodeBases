// (C) 2006, Matt Mahoney, mmahoney (at) cs.fit.edu
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
// This is an algorithm for entropy compression of the Burrows-Wheeler
// transform of a text. It was originally written by Matt Mahoney as
// bbb.cpp - big block BWT compressor version 1, Aug. 31, 2006.
// http://cs.fit.edu/~mmahoney/compression/bbb.cpp
//
// Only the entropy compression of bbb.cpp is included here.
// Other parts, including the computation of the transform, is not included.
// The entropy compressor has been fitted with an interface that
// ties it with the rest of the dcs-bwt-compressor.
//
// The class Predictor has seen the following modifications:
// - static variables in member function p() has been changed into
//   data members of Predictor.
// - The APM chain has been modified a bit.
// Here is the original documentation of the entropy coding with
// the changes to bbb.cpp added in square brackets.
//
// ENTROPY CODING
//
// BWT data is best coded with an order 0 model.  The transformed text tends
// to have long runs of identical bytes (e.g. "nnbaaa").  The BWT data is
// modeled with a modified PAQ with just one context (no mixing) followed
// by a 5 stage SSE (APM) and bitwise arithmetic coding.  Modeling typically
// takes about as much time as sorting and unsorting in slow mode.
// The model uses about 5 MB memory.
// [ Now reduced to about 256KB of memory. ]
//
// The order 0 model consists of a mapping:
//
//             order 1, 2, 3 contexts ----------+
//                                              V
//  order 0 context -> bit history -> p -> APM chain -> arithmetic coder
//                  t1             sm
//
// Bits are coded one at a time.  The arithmetic coder maintains a range
// [lo, hi), initially [0, 1) and repeatedly subdivides the range in proportion
// to p(0), p(1), the next bit probabilites predicted by the model.  The final
// output is the shortest base 256 number x such that lo <= x < hi.  As the
// leading bytes of x become known, they are output.  To decompress, the model
// predictions are repeated as during compression, then the actual bit is
// determined by which half of the subrange contains x.
//
// The model inputs a bytewise order 0 context consisting of the last 0 to 7
// bits of the current byte, plus the number of bits.  There are a total of
// 255 possible bitwise contexts.  For each context, a table (t1) maintains
// an 8 bit state representing the history of 0 and 1 bits previously seen.
// This history is mapped by another table (a StateMap sm) to a probability,
// p, that the next bit will be 1. This table is adaptive: after each
// prediction, the mapping (state -> p) is adjusted to improve the last
// prediction.
//
// The output of the StateMap is passed through a series of 6 more adaptive
// tables, (Adaptive Probability Maps, or APM) each of which maps a context
// and the input probability to an output probability.  The input probability
// is interpolated between 33 bins on a nonlinear scale with smaller bins
// near 0 and 1.  After each prediction, the corresponding table entries
// on both sides of p are adjusted to improve the last prediction.
//  The APM chain is like this:
//
//      + A11 ->+            +--->---+ +--->---+
//      |       |            |       | |       |
//  p ->+       +-> A2 -> A3 +-> A4 -+-+-> A5 -+-> Encoder
//      |       |
//      + A12 ->+
//
// [ The APM chain has been modified into:
//
//                 +--->---+ +--->---+
//                 |       | |       |
//  p --> A2 -> A3 +-> A5 -+-+-> A4 -+-> Encoder
//
// ]
//
// A11 and A12 both take c0 (the preceding bits of the current byte) as
// additional context, but one is fast adapting and the other is slow
// adapting.  Their outputs are averaged.
//
// A2 is an order 1 context (previous byte and current partial byte).
// [ A2 has been modified so that it uses only two bits of information
// from the previous byte: what is the bit in the current bit position
// and whether the preceding bits are same or different from c0. ]
//
// A3 takes the previous (but not current) byte as context, plus 2 bits
// that depend on the current run length (0, 1, 2-3, or 4+), the number
// of times the last byte was repeated.
// [ A3 now only takes the two bits on run length. ]
//
// A4 takes the current byte and the low 5 bits of the second byte back.
// The output is averaged with 3/4 weight to the A3 output with 1/4 weight.
// [ A4 has been moved after A5, it takes only the current byte (not the
// 5 additional bits), and the averaging weights are 1/2 and 1/2. ]
//
// A5 takes a 14 bit hash of an order 3 context (last 3 bytes plus
// current partial byte) and is averaged with 1/2 weight to the A4 output.
// [ A5 takes now 11 bit hash of an order 4 context. ]
//
// The StateMap, state table, APM, Encoder, and associated code (Array,
// squash(), stretch()) are taken from PAQ8 with minor non-functional
// changes (e.g. removing global context).


#ifndef DCSBWT_BBB_COMPRESS_H__
#define DCSBWT_BBB_COMPRESS_H__

#include "stream_compressor.h"
#include "inttypes.h"

#include <string>
#include <cassert>

namespace dcsbwt {

namespace bbb {

//////////////////////////////////////////////////////////////////
// Array<T, ALIGN> a(n); creates n elements of T initialized to 0 bits.
// Constructors for T are not called.
// Indexing is bounds checked if assertions are on.
// a.size() returns n.
// a.resize(n) changes size to n, padding with 0 bits or truncating.
// Copy and assignment are not supported.
// Memory is aligned on a ALIGN byte boundary (power of 2), default is none.
//////////////////////////////////////////////////////////////////
template <class T, int ALIGN=0> class Array {
public:
  explicit Array(int i=0) { create(i);}
  ~Array();
  T& operator[](int i) {
    assert(i >= 0);
    assert(i < n);
    return data[i];
  }
  const T& operator[](int i) const {
    assert(i >= 0);
    assert(i < n);
    return data[i];
  }
  int size() const { return n; }
  void resize(int i);  // change size to i
private:
  int n;     // user size
  int reserved;  // actual size
  char *ptr; // allocated memory, zeroed
  T* data;   // start of n elements of aligned data

  void create(int i);  // create with size i

  Array(const Array&);  // no copy or assignment
  Array& operator=(const Array&);
};

template<class T, int ALIGN> void Array<T, ALIGN>::resize(int i) {
  if (i<=reserved) {
    n=i;
    return;
  }
  char *saveptr=ptr;
  T *savedata=data;
  int saven=n;
  create(i);
  if (savedata && saveptr) {
    memcpy(data, savedata, sizeof(T)*std::min(i, saven));
    free(saveptr);
  }
}

template<class T, int ALIGN> void Array<T, ALIGN>::create(int i) {
  n=reserved=i;
  if (i<=0) {
    data=0;
    ptr=0;
    return;
  }
  const int sz=ALIGN+n*sizeof(T);
  ptr = (char*)calloc(sz, 1);
  if (!ptr) fprintf(stderr, "Out of memory\n"), exit(1);
  data = (ALIGN ? (T*)(ptr+ALIGN-(((long)ptr)&(ALIGN-1))) : (T*)ptr);
  assert((char*)data >= ptr);
  assert((char*)data <= ptr+ALIGN);
}

template<class T, int ALIGN> Array<T, ALIGN>::~Array() {
  free(ptr);
}

//////////////////////////////////////////////////////////////////
// A StateMap maps a nonstationary counter state to a probability.
// After each mapping, the mapping is adjusted to improve future
// predictions.  Methods:
//
// sm.p(y, cx) converts state cx (0-255) to a probability (0-4095),
//   and trains by updating the previous prediction with y (0-1).
//
// Counter state -> probability * 256
//////////////////////////////////////////////////////////////////
class StateMap {
public:
  StateMap();
  int p(int y, int cx) {
    assert(cx>=0 && cx<t.size());
    t[cxt]+=((y<<16)-t[cxt]+128) >> 8;
    return t[cxt=cx] >> 4;
  }
 private:
  int cxt;  // context
  Array<uint16> t; // 256 states -> probability * 64K
};


//////////////////////////////////////////////////////////////////
// APM maps a probability and a context into a new probability
// that bit y will next be 1.  After each guess it updates
// its state to improve future guesses.  Methods:
//
// APM a(N) creates with N contexts, uses 66*N bytes memory.
// a.p(y, pr, cx, rate=8) returned adjusted probability in context cx (0 to
//   N-1).  rate determines the learning rate (smaller = faster, default 8).
//   Probabilities are scaled 12 bits (0-4095).  Update on last bit y (0-1).
//////////////////////////////////////////////////////////////////
class APM {
 public:
  explicit APM(int n);
  int p(int y, int pr=2048, int cxt=0, int rate=8);

 private:
  int index;     // last p, context
  const int N;   // number of contexts
  Array<uint16> t;  // [N][33]:  p, context -> p
};

//////////////////////////// Predictor //////////////////////////
class Predictor {
public:
  Predictor(): pr(2048), c0(1), c4(0), bpos(0), t1(256), cp(&t1[0]),
               run(0), runcxt(0),
               // a11(256), a12(256), a2(65536), a3(1024), a4(8192), a5(16384)
               a2(1024), a3(4), a5(2048), a4(256) {}

  int p() const {return pr;}
  void update(int y);

 private:
  int pr;  // next return value of p() (0-4095)
  int c0;  // bitwise context: last 0-7 bits with a leading 1 (1-255)
  uint32 c4;  // last 4 whole bytes, last is in low 8 bits
  int bpos; // number of bits in c0 (0-7)
  Array<uint8> t1; // context -> state
  StateMap sm;  // state -> pr
  uint8* cp;  // context pointer
  int run;  // count of consecutive identical bytes (0-65535)
  int runcxt;  // (0-3) if run is 0, 1, 2-3, 4+
  //APM a11, a12, a2, a3, a4, a5;
  APM a2, a3, a5, a4;
};

}  // namespace bbb

class BbbCompressor : public StreamCompressor {
 public:
  class Options : public StreamCompressor::OptionsBase {
   public:
    virtual bool Set(const std::string& options) {
      return 0 == options.length();
    }
    virtual std::string Get() const { return std::string(); }
    int64 SizeInBytes() const { return 256 << 10; }
    virtual StreamCompressor* GetCompressor() { return new BbbCompressor; }
  };

  BbbCompressor() {}
  virtual ~BbbCompressor() { }

  virtual void WriteBegin();
  virtual void WriteEndPrivate();
  virtual void Write(const char* bytes, size_t n) {
    for (; n; --n) { CompressByte(*bytes++); }
  }

 private:
  bbb::Predictor predictor_;
  uint32 x1_;
  uint32 x2_;


  void CompressByte(unsigned char byte) {
    for (int i=7; i>=0; --i) { CompressBit((byte>>i)&1); }
  }
  void CompressBit(int bit);
};


class BbbDecompressor : public StreamDecompressor {
 public:
  static StreamDecompressor* Create(const std::string&) {
    return new BbbDecompressor;
  }
  BbbDecompressor() {}
  virtual ~BbbDecompressor() {}

  virtual void ReadBegin();
  virtual void ReadEnd() {}
  virtual void Read(char* bytes, size_t n) {
    for (; n; --n) *bytes++ = GetUncompressedByte();
  }

  virtual int64 SizeInBytes() const { return 6 << 20; }

 private:
  bbb::Predictor predictor_;
  uint32 x1_;
  uint32 x2_;
  uint32 x_;

  unsigned char GetUncompressedByte() {
    int byte = GetUncompressedBit();
    for (int i = 7; i; --i) byte = (byte<<1) + GetUncompressedBit();
    return byte;
  }
  int GetUncompressedBit();
};

}  // namespace dcsbwt

#endif  // DCSBWT_BBB_COMPRESS_H__
