// Copyright 2007 Google Inc.

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#include "stream.h"

#include <algorithm>  // for copy
#include <cassert>

namespace dcsbwt {

void OutStreamBuffer::FlushAndWrite(const char* bytes, size_t n) {
  Flush();
  if (n < buffer_.size() / 2) {
    WriteToBuffer(bytes, n);
  } else {
    master_->Write(bytes, n);
  }
}

void OutStreamBuffer::Flush() {
  assert(IsConnected());
  if (buffer_.begin() < next_free_slot_) {
    master_->Write(&buffer_[0], next_free_slot_ - buffer_.begin());
    next_free_slot_ = buffer_.begin();
  }
}

void OutStreamBuffer::Reset(size_t size) {
  Flush();
  if (size <= 0) size = 1;
  std::vector<char>(size).swap(buffer_);
  next_free_slot_ = buffer_.begin();
}


void InStreamBuffer::ReadAndRefill(char* bytes, size_t n) {
  size_t available = AvailableInBuffer();
  assert(n > available);
  ReadFromBuffer(bytes, available);
  bytes += available;
  n -= available;
  if (n < buffer_.size() / 2) {
    Refill();
    ReadFromBuffer(bytes, n);
  } else {
    master_->Read(bytes, n);
  }
}

void InStreamBuffer::Refill() {
  assert(AvailableInBuffer() == 0);
  master_->Read(&buffer_[0], buffer_.size());
  next_unused_byte_ = buffer_.begin();
}

void InStreamBuffer::Reset(size_t size) {
  if (size <= 0) size = 1;
  std::vector<char>(size).swap(buffer_);
  next_unused_byte_ = buffer_.end();
}

}  // namespace dcsbwt
