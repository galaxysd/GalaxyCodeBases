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

#include "../src/stringsort-inl.h"
#include "../src/stringsort.h"
#include "../src/ternary_partition.h"  // only for random numbers
#include "../src/inttypes.h"

#include <iostream>
#include <string>
#include <functional>
#include <vector>
#include <algorithm>
#include <cassert>

namespace {

using dcsbwt::StringsortSuffixes;
using dcsbwt::ternary_partition::RandomNumberGenerator;

void PermuteRandomly(uint32* begin, uint32* end, RandomNumberGenerator* rng) {
  int size = end - begin;
  for (int i = 0; i < size; ++i) {
    std::swap(begin[i], (begin+i)[rng->Uniform(size-i)]);
  }
}

class SortedGroupVerifier {
 public:
  SortedGroupVerifier(const unsigned char* text, const unsigned char* text_end,
                      uint32* suffixes, uint32* suffixes_end,
                      int prefix_length)
      : text_(text), text_length_(text_end - text),
        prefix_length_(prefix_length),
        suffixes_(suffixes),
        suffixes_end_(suffixes_end),
        previous_begin_(suffixes), previous_end_(suffixes),
        report_count_(text_end - text + 1, 0) {}
  void Reset() {
    previous_begin_ = suffixes_;
    previous_end_ = suffixes_;
    std::fill(report_count_.begin(), report_count_.end(), 0);
  }
  void FinalCheck() {
    assert(previous_end_ - suffixes_ == suffixes_end_ - suffixes_);
    for (const uint32* it = suffixes_; it < suffixes_end_; ++it) {
      assert(report_count_[*it] == 1);
    }
  }
  void operator() (const uint32* begin, const uint32* end) {
    assert(begin - suffixes_ == previous_end_ - suffixes_);
    assert(suffixes_end_ - suffixes_ >= end - suffixes_);
    const int size = end - begin;
    assert(size > 0);
    assert(*begin <= text_length_);
    for (const uint32* it = begin; it < end; ++it) {
      ++report_count_[*it];
    }
    // Check that all group suffixes share a prefix of length prefix_length_
    if (size > 1) {
      int full_length_suffix = text_length_ - prefix_length_;
      assert(*begin <= full_length_suffix);
      for (const uint32* it = begin + 1; it < end; ++it) {
        assert(*it <= full_length_suffix);
        for (int i = 0; i < prefix_length_; ++i) {
          assert(text_[*begin + i] == text_[*it + i]);
        }
      }
    }
    // Check that this group's suffixes are larger than previous group's
    // First check whether there is a previous group.
    if (previous_begin_ < previous_end_) {
      // Compare the first elements of the groups.
      int limit = prefix_length_;
      limit = std::min<int>(limit, text_length_ - *previous_begin_);
      limit = std::min<int>(limit, text_length_ - *begin);
      int i;
      for (i = 0; i < limit; ++i)
        if (text_[*previous_begin_ + i] != text_[*begin + i]) break;
      assert(i < text_length_ - *begin);
      if (i < limit) {
        assert(text_[*previous_begin_ + i] < text_[*begin + i]);
      } else if (i < text_length_ - *previous_begin_) {
        // The two groups share a prefix of length prefix_length_.
        // Now we must compare all-against-all all the way.
        for (const uint32* it1 = previous_begin_; it1 < previous_end_; ++it1) {
          for (const uint32* it2 = begin; it2 < end; ++it2) {
            int length1 = text_length_ - *it1;
            int length2 = text_length_ - *it2;
            limit = std::min<int>(length1, length2);
            for (i = 0; i < limit; ++i)
              if (text_[*it1 + i] != text_[*it2 + i]) break;
            assert(i < length2);
            if (i < length1) {
              assert(text_[*it1 + i] < text_[*it2 + i]);
            }
          }
        }
      }
    }
    previous_begin_ = begin;
    previous_end_ = end;
  }
 private:
  const unsigned char* const text_;
  const int text_length_;
  const int prefix_length_;
  const uint32* const suffixes_;
  const uint32* const suffixes_end_;
  const uint32* previous_begin_;
  const uint32* previous_end_;
  std::vector<int> report_count_;
};

void SortAndVerify(const std::string& str,
                   uint32 target_prefix_length) {
  int text_length = str.size();
  std::string::const_iterator text = str.begin();
  std::string::const_iterator text_end = text + text_length;

  int num_suffixes = text_length + 1;
  std::vector<uint32> suffix_array(2 * num_suffixes);
  for (int i = 0; i < num_suffixes; ++i) suffix_array[i] = i;
  uint32* suffix_area = &suffix_array[0];
  uint32* suffixes_end = suffix_area + num_suffixes;
  uint32* full_area_end = suffix_area + suffix_array.size();
  uint32* small_area_end = suffixes_end + (full_area_end - suffixes_end)/4;

  const unsigned char* text1
      = reinterpret_cast<const unsigned char*>(str.data());
  const unsigned char* text_end1 = text1 + text_length;
  SortedGroupVerifier verifier(text1, text_end1, suffix_area, suffixes_end,
                               target_prefix_length);
  RandomNumberGenerator rng;
  // no workspace
  PermuteRandomly(suffix_area, suffixes_end, &rng);
  verifier.Reset();
  uint32* result =
      StringsortSuffixes(text, text_end,
                         suffix_area, suffix_area, suffixes_end,
                         0, target_prefix_length,
                         verifier);
  assert(result == suffixes_end);
  verifier.FinalCheck();
  // some workpace
  PermuteRandomly(suffix_area, suffixes_end, &rng);
  verifier.Reset();
  uint32* suffixes_in =
      std::copy_backward(suffix_area, suffixes_end, small_area_end);
  result =
      StringsortSuffixes(text, text_end,
                     suffix_area, suffixes_in, small_area_end,
                         0, target_prefix_length,
                         verifier);
  assert(result == suffixes_end);
  verifier.FinalCheck();
  // full workpace
  PermuteRandomly(suffix_area, suffixes_end, &rng);
  verifier.Reset();
  suffixes_in =
      std::copy_backward(suffix_area, suffixes_end, full_area_end);
  result =
      StringsortSuffixes(text, text_end,
                         suffix_area, suffixes_in, full_area_end,
                         0U, target_prefix_length,
                         verifier);
  assert(result == suffixes_end);
  verifier.FinalCheck();
}

void TestEmpty() {
  std::string str;
  SortAndVerify(str, 1);
  std::cout << "TestEmpty passed" << std::endl;
}

void TestLengthOne() {
  std::string str1("\0");
  SortAndVerify(str1, 1);
  std::string str2 = "\xFF";
  SortAndVerify(str2, 1);
  std::cout << "TestLengthOne passed" << std::endl;
}

void TestBanana() {
  std::string str("banana");
  SortAndVerify(str, 6);
  SortAndVerify(str, 3);
  std::cout << "TestBanana passed" << std::endl;
}

char const * const dna_512 =
"CACTGTGCCATCATCATCACCACCACTGTCATTATCACCACCACCATCATCACCAACACCACTG"
"CCATCGTCATCACCACCACTGTCATTATCACCACCACCATCACCAACATCACCACCACCATTAT"
"CACCACCATCAACACCACCACCCCCATCATCATCATCACTACTACCATCATTACCAGCACCACC"
"ACCACTATCACCACCACCACCACAATCACCATCACCACTATCATCAACATCATCACTACCACCA"
"TCACCAACACCACCATCATTATCACCACCACCACCATCACCAACATCACCACCATCATCATCAC"
"CACCATCACCAAGACCATCATCATCACCATCACCACCAACATCACCACCATCACCAACACCACC"
"ATCACCACCACCACCACCATCATCACCACCACCACCATCATCATCACCACCACCGCCATCATCA"
"TCGCCACCACCATGACCACCACCATCACAACCATCACCACCATCACAACCACCATCATCACTAT";

void TestDNA() {
  std::string str(dna_512);
  SortAndVerify(str, 512);
  SortAndVerify(str, 3);
  std::cout << "TestDNA passed" << std::endl;
}

char const * const aaa_256 =
"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

void TestAAAA() {
  std::string str(aaa_256);
  SortAndVerify(str, 256);
  SortAndVerify(str, 13);
  std::cout << "TestAAAA passed" << std::endl;
}


void TestAAAB() {
  std::string str("AAAAAAAAB");
  SortAndVerify(str, 9);
  std::cout << "TestAAAB passed" << std::endl;
}

void TestRandomLargeAlphabet() {
  std::string str(1<<14, '\0');
  RandomNumberGenerator rng;
  for (std::string::iterator i = str.begin(); i != str.end(); ++i) {
    *i = rng.Uniform(256);
  }
  SortAndVerify(str, 1<<14);
  std::cout << "TestRandomLargeAlphabet passed" << std::endl;
}

void TestRandomSmallAlphabet() {
  std::string str(1<<14, '\0');
  RandomNumberGenerator rng;
  for (std::string::iterator i = str.begin(); i != str.end(); ++i) {
    *i = 65 + rng.Uniform(2);
  }
  SortAndVerify(str, 1<<14);
  std::cout << "TestRandomSmallAlphabet passed" << std::endl;
}

void TestRandomLimited() {
  std::string str(1<<14, '\0');
  RandomNumberGenerator rng;
  for (std::string::iterator i = str.begin(); i != str.end(); ++i) {
    *i = 65 + rng.Uniform(2);
  }
  SortAndVerify(str, 2);
  std::cout << "TestRandomLimited passed" << std::endl;
}

void TestNoSuffixes() {
  std::string str(1<<14, '\0');
  RandomNumberGenerator rng;
  for (std::string::iterator i = str.begin(); i != str.end(); ++i) {
    *i = 65 + rng.Uniform(2);
  }
  SortAndVerify(str, 0);
  std::cout << "TestNoSuffixes passed" << std::endl;
}

}  // namespace

int main(int argc, char **argv) {

  TestEmpty();
  TestLengthOne();
  TestBanana();
  TestDNA();
  TestAAAA();
  TestAAAB();
  TestRandomLargeAlphabet();
  TestRandomSmallAlphabet();
  TestRandomLimited();
  TestNoSuffixes();

  return 0;
}
