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

#include <iostream>
#include <string>
#include <functional>
#include <vector>
#include <algorithm>
#include <cassert>

#include "../src/prefix_doubling-inl.h"
#include "../src/prefix_doubling.h"

namespace {

using dcsbwt::SortSuffixesByPrefixDoubling;
using dcsbwt::PrefixDoubler;

// This verification is based on Theorem 2 in:
//
//   Burkhardt, Karkkainen:
//   Fast lightweight suffix array construction and checking.
//   Proc. 14th Symposium on Combinatorial Pattern Matching,
//   LNCS 2676, Springer, 2003, pp. 55-69.
//
void Verify(std::string const& str, std::vector<int> const& rank_array,
            std::vector<int> const& suffix_array) {
  int str_length = str.length();
  assert(suffix_array[0] == str_length);
  assert(rank_array[str_length] == 0);
  for (int rank = 1; rank <= str_length; ++rank) {  // yes, I mean "<="
    int suffix = suffix_array[rank];
    assert(suffix >= 0);
    assert(suffix < str_length);
    assert(rank_array[suffix] == rank);
  }
  for (int rank = 1; rank < str_length; ++rank) {  // yes, I mean "<"
    int suffix1 = suffix_array[rank];
    int suffix2 = suffix_array[rank + 1];
    unsigned char first_char1 = static_cast<unsigned char>(str[suffix1]);
    unsigned char first_char2 = static_cast<unsigned char>(str[suffix2]);
    assert(first_char1 <= first_char2);
    if (first_char1 == first_char2) {
      assert(rank_array[suffix1 + 1] < rank_array[suffix2 + 1]);
    }
  }
}

void ConstructAndVerify(std::string const& str) {
  std::vector<int> rank_array(str.size() + 1);
  std::vector<int> suffix_array(str.size() + 1);
  SortSuffixesByPrefixDoubling(str.begin(), str.end(),
                               rank_array.begin(), rank_array.end(),
                               suffix_array.begin(), suffix_array.end());
  Verify(str, rank_array, suffix_array);
}

void TestEmpty() {
  std::string str;
  ConstructAndVerify(str);
  std::cout << "TestEmpty passed" << std::endl;
}


void TestLengthOne() {
  std::string str("\0");
  ConstructAndVerify(str);
  str[0] = '\xFF';
  ConstructAndVerify(str);
  std::cout << "TestLengthOne passed" << std::endl;
}


void TestBanana() {
  std::string str("banana");
  ConstructAndVerify(str);
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
  ConstructAndVerify(str);
  std::cout << "TestDNA passed" << std::endl;
}

char const * const aaa_256 =
"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

void TestAAAA() {
  std::string str(aaa_256);
  std::vector<int> rank_array(str.size() + 1);
  std::vector<int> suffix_array(str.size() + 1);
  SortSuffixesByPrefixDoubling(str.begin(), str.end(),
                               rank_array.begin(), rank_array.end(),
                               suffix_array.begin(), suffix_array.end());
  Verify(str, rank_array, suffix_array);
  std::vector<int>::iterator it =
      std::adjacent_find(suffix_array.begin(), suffix_array.end(),
                         std::less_equal<int>());
  assert(suffix_array.end() - it == 0);
  std::cout << "TestAAAA passed" << std::endl;
}

void TestAAAB() {
  std::string str("AAAAAAAAB");
  std::vector<int> rank_array(str.size() + 1);
  std::vector<int> suffix_array(str.size() + 1);
  SortSuffixesByPrefixDoubling(str.begin(), str.end(),
                               rank_array.begin(), rank_array.end(),
                               suffix_array.begin(), suffix_array.end());
  Verify(str, rank_array, suffix_array);
  std::vector<int>::iterator it =
      std::adjacent_find(suffix_array.begin() + 1, suffix_array.end(),
                         std::greater_equal<int>());
  assert(suffix_array.end() - it == 0);
  std::cout << "TestAAAB passed" << std::endl;
}

// Two strings with the same order of suffixes.
void TestSameOrder() {
  std::string str1("AAAAAAAAB");
  std::vector<int> rank_array(str1.size() + 1);
  std::vector<int> suffix_array(str1.size() + 1);
  SortSuffixesByPrefixDoubling(str1.begin(), str1.end(),
                               rank_array.begin(), rank_array.end(),
                               suffix_array.begin(), suffix_array.end());
  std::string str2("XXXXXXXXY");
  Verify(str2, rank_array, suffix_array);
  std::cout << "TestSameOrder passed" << std::endl;
}

// Check that '\xFF' > 'A'
void TestAscii255() {
  std::string str1("AAAAAAAAB");
  std::vector<int> rank_array(str1.size() + 1);
  std::vector<int> suffix_array(str1.size() + 1);
  SortSuffixesByPrefixDoubling(str1.begin(), str1.end(),
                               rank_array.begin(), rank_array.end(),
                               suffix_array.begin(), suffix_array.end());
  std::string str2("AAAAAAAA\xFF");
  Verify(str2, rank_array, suffix_array);
  std::cout << "TestAscii255 passed" << std::endl;
}

void TestComputeOneOrder() {
  std::string str(dna_512);
  std::vector<int> rank_array(str.size() + 1);
  std::vector<int> suffix_array(str.size() + 1);
  PrefixDoubler<int> pd(str.size());
  pd.ComputeOneOrder(str.begin(), str.end(),
                     rank_array.begin(), rank_array.end(),
                     suffix_array.begin(), suffix_array.end());
  pd.SortSuffixes(rank_array.begin(), rank_array.end(),
                  suffix_array.begin(), suffix_array.end(), 1);
  pd.ComputeFinalSuffixArray(rank_array.begin(), rank_array.end(),
                             suffix_array.begin(), suffix_array.end());
  Verify(str, rank_array, suffix_array);
  std::cout << "TestComputeOneOrder passed" << std::endl;
}

void TestDouble() {
  std::string str(dna_512);
  std::vector<int> rank_array(str.size() + 1);
  std::vector<int> suffix_array(str.size() + 1);
  PrefixDoubler<int> pd(str.size());
  pd.ComputeOneOrder(str.begin(), str.end(),
                     rank_array.begin(), rank_array.end(),
                     suffix_array.begin(), suffix_array.end());
  int order_depth = pd.Double(rank_array.begin(), rank_array.end(),
                              suffix_array.begin(), suffix_array.end(), 1);
  pd.SortSuffixes(rank_array.begin(), rank_array.end(),
                  suffix_array.begin(), suffix_array.end(), order_depth);
  pd.ComputeFinalSuffixArray(rank_array.begin(), rank_array.end(),
                             suffix_array.begin(), suffix_array.end());
  Verify(str, rank_array, suffix_array);
  std::cout << "TestDouble passed" << std::endl;
}

void TestComputeInitialSuffixArray() {
  std::string str(dna_512);
  std::vector<int> rank_array(str.size() + 1);
  std::vector<int> suffix_array(str.size() + 1);
  PrefixDoubler<int> pd(str.size());
  pd.ComputeOneOrder(str.begin(), str.end(),
                     rank_array.begin(), rank_array.end(),
                     suffix_array.begin(), suffix_array.end());
  pd.ComputeInitialSuffixArray(rank_array.begin(), rank_array.end(),
                  suffix_array.begin(), suffix_array.end(), 1);
  pd.SortSuffixes(rank_array.begin(), rank_array.end(),
                  suffix_array.begin(), suffix_array.end(), 1);
  pd.ComputeFinalSuffixArray(rank_array.begin(), rank_array.end(),
                             suffix_array.begin(), suffix_array.end());
  Verify(str, rank_array, suffix_array);
  std::cout << "TestComputeInitialSuffixArray passed" << std::endl;
}

}  // namespace

int main(int argc, char **argv) {

  TestEmpty();
  TestLengthOne();
  TestBanana();
  TestDNA();
  TestAAAA();
  TestAAAB();
  TestSameOrder();
  TestAscii255();
  TestComputeOneOrder();
  TestComputeInitialSuffixArray();

  return 0;
}
