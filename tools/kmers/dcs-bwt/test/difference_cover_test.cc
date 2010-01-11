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

#include "../src/difference_cover-inl.h"
#include "../src/difference_cover.h"
#include "../src/inttypes.h"
#include "../src/ternary_partition.h"  // for random numbers

#include <cstdio>
#include <vector>
#include <cassert>
#include <iostream>

namespace {

using dcsbwt::DifferenceCover;
using dcsbwt::DifferenceCoverSample;
using dcsbwt::ternary_partition::RandomNumberGenerator;

const int kMaxPeriod = 1024;  // includes all precomputed covers and
                              // two covers computed by Colbourn-Ling

void TestCoverConstruction() {
  for (int32 period = 1; period <= kMaxPeriod; period *= 2) {
    DifferenceCover dc(period);
    assert(dc.period() == period);
    assert(dc.size() <= period);
    assert(DifferenceCover::CoverSizeForPeriod(period) == dc.size());
    assert(DifferenceCover::SizeInBytesForPeriod(period) > period);
  }
  std::cout << "TestCoverConstruction passed" << std::endl;
}


void TestCoverDivision() {
  for (int32 period = 1; period <= kMaxPeriod; period *= 2) {
    DifferenceCover dc(period);
    for (int32 dividend = 0; dividend < period; ++dividend) {
      assert(dc.DivideByPeriod(dividend) == 0);
      assert(dc.ModuloPeriod(dividend) == dividend);
    }
    int32 dividend = 0;
    for (int32 quotient = 0; quotient < period; ++quotient) {
      dividend += period - 1;
      int32 remainder = period - quotient - 1;
      assert(dc.DivideByPeriod(dividend) == quotient);
      assert(dc.ModuloPeriod(dividend) == remainder);
    }
    // Division should work even for values larger than MAX_UINT.
    int64 large_value = 1489245697U;  // random large number
    large_value *= (large_value-4);
    assert(dc.DivideByPeriod(large_value) * static_cast<int64>(period) +
              dc.ModuloPeriod(large_value) == large_value);
    // ModuloPeriod should work even for negative values.
    int x = -1;
    assert(dc.ModuloPeriod(x) == period-1);
  }
  std::cout << "TestCoverDivision passed" << std::endl;
}


void TestCoverRankSelect() {
  for (int32 period = 1; period <= kMaxPeriod; period *= 2) {
    DifferenceCover dc(period);
    int32 i = 0;
    for (; i < period; ++i) {
      // continue up to and including the last cover element
      if (dc.Rank(i) == dc.size()) break;
      assert(dc.Rank(i) < dc.size());
      assert(dc.Select(dc.Rank(i)) >= i);
    }
    for (; i < period; ++i) {
      // after the last cover element
      assert(dc.Rank(i) == dc.size());
    }
    for (int32 j = 0; j < dc.size(); ++j) {
      assert(dc.Select(j) < period);
      assert(dc.Rank(dc.Select(j)) == j);
      assert(dc.Contains(dc.Select(j)));
    }
  }
  std::cout << "TestCoverRankSelect passed" << std::endl;
}


void TestCoverContainer() {
  for (int32 period = 1; period <= kMaxPeriod; period *= 2) {
    DifferenceCover dc(period);
    std::vector<int32> cover(dc.begin(), dc.end());
    assert(cover.size() == dc.size());
    std::vector<int32>::iterator it = cover.begin();
    for (int32 i = 0; i < period; ++i) {
      if (it != cover.end() && *it == i) {
        assert(dc.Contains(i));
        ++it;
      } else {
        assert(!dc.Contains(i));
      }
    }
  }
  std::cout << "TestCoverContainer passed" << std::endl;
}


void TestCoverCoverage() {
  for (int32 period = 1; period <= kMaxPeriod; period *= 2) {
    DifferenceCover dc(period);
    for (int32 diff = 0; diff < period; ++diff) {
      int32 k = dc.Coverer(diff);
      assert(k < period);
      assert(dc.Contains(k));
      int32 k_plus_diff = dc.ModuloPeriod(k + diff);
      assert(dc.Contains(k_plus_diff));
    }
  }
  std::cout << "TestCoverCoverage passed" << std::endl;
}


static struct PeriodRangePair{
  int32 period;
  int32 range;
} period_range_pairs[20] = {
  {    1,   1 },
  {    1,   2 },
  {    1,   7 },
  {    2,   2 },
  {    2,   3 },
  {    2,   4 },
  {    2,   7 },
  {   64,  64 },
  {   64,  65 },
  {   64,  96 },
  {   64,  127 },
  {   64,  128 },
  {   64,  129 },
  {   64,  917 },
  { 1024, 1024 },
  { 1024, 1025 },
  { 1024, 2047 },
  { 1024, 2048 },
  { 1024, 2049 },
  { 1024, 5847 }
};

void TestSampleConstruction() {
  for (int32 i = 0; i < 20; ++i) {
    int32 period = period_range_pairs[i].period;
    int32 range = period_range_pairs[i].range;
    DifferenceCoverSample<int32> dcs(period, range);
    assert(dcs.period() == period);
    assert(dcs.range() == range);
    assert(dcs.size() <= range);
    assert(dcs.period_size() <= period);
  }
  std::cout << "TestSampleConstruction passed" << std::endl;
}


void TestSampleFill() {
  for (int32 i = 0; i < 20; ++i) {
    int32 period = period_range_pairs[i].period;
    int32 range = period_range_pairs[i].range;
    DifferenceCoverSample<int32> dcs(period, range);
    std::vector<int32> sample;
    dcs.Fill(std::back_inserter(sample));
    assert(sample.size() == dcs.size());
    std::vector<int32>::iterator it = sample.begin();
    int32 sample_count = 0;
    bool death_test_done = false;     // do just one death test
    for (int32 i = 0; i < range; ++i) {
      if (it != sample.end() && *it == i) {
        assert(dcs.Contains(i));
        assert(dcs.Rank(i) == sample_count);
        if (i < range - period) {
          assert(dcs.Rank(i+period) - dcs.Rank(i) == dcs.PeriodInterval());
        }
        ++it;
        ++sample_count;
      } else {
        assert(!dcs.Contains(i));
      }
    }
  }
  std::cout << "TestSampleFill passed" << std::endl;
}


void TestSampleShift() {
  for (int32 i = 0; i < 20; ++i) {
    int32 period = period_range_pairs[i].period;
    int32 range = period_range_pairs[i].range;
    DifferenceCoverSample<int32> dcs(period, range);
    RandomNumberGenerator rng;
    for (int32 n = 20000; n; --n) {
      int32 i = rng.Uniform(range);
      int32 j = rng.Uniform(range);
      int32 shift = dcs.Shift(i,j);
      assert(shift >= 0);
      assert(shift < period);
      if (i < range - shift) {
        assert(dcs.Contains(i + shift));
      }
      if (j < range - shift) {
        assert(dcs.Contains(j + shift));
      }
    }
  }
  std::cout << "TestSampleShift passed" << std::endl;
}

}  // unnamed namespace

int main(int argc, char **argv) {

  TestCoverConstruction();
  TestCoverDivision();
  TestCoverRankSelect();
  TestCoverContainer();
  TestCoverCoverage();
  TestSampleConstruction();
  TestSampleFill();
  TestSampleShift();

  return 0;
}
