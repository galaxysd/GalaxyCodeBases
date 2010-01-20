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

// Base class for doing the Burrows-Wheeler transform.
// Algorithms the Burrows-Wheeler transform should
// have an interface derived from the class BWTransformer.
//
// New algorithms can be added by modifying the function
// BWTransformer::Options::Set in the file bw_transform.cc

#ifndef DCSBWT_BW_TRANSFORM_H__
#define DCSBWT_BW_TRANSFORM_H__

#include "stream.h"
#include "inttypes.h"

#include <string>
#include <cstddef>  // for size_t

namespace dcsbwt {

class BWTransformer {
 protected:
  class AlgorithmSpecificOptions;

 public:
  // A class representing options affecting the BW transform
  class Options {
   public:
    Options() : algorithm_specific_options_(NULL) { Set("p"); }

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

    // Get() returns a string that as an argument to Set() sets
    // any options object to this object's current state.
    std::string Get() const {
      if (algorithm_specific_options_) {
        return algorithm_id_ + algorithm_specific_options_->Get();
      } else {
        return std::string();
      }
    }

    // MaxSizeInBytes() returns the peak memory consumption
    // during the BW transform of a block of size block_size
    // under the current options when the memory budget is
    // unlimited or large enough.
    // The size includes the input block but not the output.
    int64 MaxSizeInBytes(int64 block_size) const {
      return algorithm_specific_options_->MaxSizeInBytes(block_size);
    }

    // MaxBlockSize() returns the largest block size that the algorithm
    // specified by current options can handle within memory_budget.
    int64 MaxBlockSize(int64 memory_budget) const {
      return algorithm_specific_options_->MaxBlockSize(memory_budget);
    }

    // SuggestedBlockSize() returns the "optimal" block size
    // for the given memory_budget under the current options.
    // It may be possible to use a larger block size but at
    // a (significantly) reduced speed, and the algorithm might
    // run (slightly) faster using a smaller block size.
    int64 SuggestedBlockSize(int64 memory_budget) const {
      return algorithm_specific_options_->SuggestedBlockSize(memory_budget);
    }

   private:
    friend class BWTransformer;
    BWTransformer* GetAlgorithm(int64 memory_budget) {
      return algorithm_specific_options_->GetTransformer(memory_budget);
    }

    char algorithm_id_;
    AlgorithmSpecificOptions* algorithm_specific_options_;
  };

  static int64 Transform(char* block, size_t block_size,
                         OutStream* output,
                         Options options,
                         int64 memory_budget);

  BWTransformer() {}
  virtual ~BWTransformer() {}

 protected:
  // A class derived from BWTransform should
  // have a nested class derived from AlgorithmSpecificOptions.
  class AlgorithmSpecificOptions {
   public:
    AlgorithmSpecificOptions() {}
    virtual ~AlgorithmSpecificOptions() {}

    // Set/Get should behave as BWTransform::Options::Set/Get
    // The algorithm identifying initial character is omitted from
    // the input to Set and should be omitted from the output of Get.
    virtual bool Set(const std::string& options_string) =0;
    virtual std::string Get() const =0;
    virtual int64 MaxSizeInBytes(int64 block_size) const =0;
    virtual int64 MaxBlockSize(int64 memory_budget) const =0;
    virtual int64 SuggestedBlockSize(int64 memory_budget) const =0;

    // Returns a transformer object corresponding to the current options.
    virtual BWTransformer* GetTransformer(int64 memory_budget) =0;
  };

 private:
  virtual int64 DoTransform(char* block, size_t block_size,
                            OutStream* output) =0;

  BWTransformer(const BWTransformer&);
  BWTransformer& operator=(const BWTransformer&);
};

}  // namespace dcsbwt

#endif  // DCSBWT_BW_TRANSFORM_H__
