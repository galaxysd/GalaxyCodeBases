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

// Base class for doing the inverse Burrows-Wheeler transform.
// Algorithms the inverse Burrows-Wheeler transform should
// have an interface derived from the class InverseBWTransformer.
//
// New algorithms can be added by modifying the function
// InverseBWTransformer::Options::Set in the file inverse_bwt.cc


#ifndef DCSBWT_INVERSE_BWT_H__
#define DCSBWT_INVERSE_BWT_H__

#include "stream.h"
#include "inttypes.h"

#include <string>

namespace dcsbwt {

/////////////////////////////////////////////////////////////////
// InverseBWTransformer
/////////////////////////////////////////////////////////////////
class InverseBWTransformer {
 protected:
  class AlgorithmSpecificOptions;

 public:
  class Options {
   public:
    Options() : algorithm_specific_options_(NULL) { Set("f"); }

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

    // MaxBlockSize() returns the largest block size that the algorithm
    // specified by current options can handle within memory_budget.
    int64 MaxBlockSize(int64 memory_budget) const {
      return algorithm_specific_options_->MaxBlockSize(memory_budget);
    }

   private:
    friend class InverseBWTransformer;
    InverseBWTransformer* GetAlgorithm(int64 memory_budget) {
      return algorithm_specific_options_->GetTransformer(memory_budget);
    }

    char algorithm_id_;
    AlgorithmSpecificOptions* algorithm_specific_options_;
  };

  InverseBWTransformer() {}
  virtual ~InverseBWTransformer() {}

  static bool Transform(const char* bwt_block, size_t bwt_block_size,
                        size_t eob_position, OutStream* output,
                        Options options,
                        int64 memory_budget);

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
    virtual int64 MaxBlockSize(int64 memory_budget) const =0;

    // Returns a transformer object corresponding to the current options.
    virtual InverseBWTransformer* GetTransformer(int64 memory_budget) =0;
  };

 private:
  virtual bool DoTransform(const char* bwt, size_t bwt_size,
                            size_t eob_position, OutStream* output) =0;

  InverseBWTransformer(const InverseBWTransformer&);
  InverseBWTransformer& operator=(const InverseBWTransformer&);
};

}  // namespace dcsbwt

#endif  // DCSBWT_INVERSE_BWT_H__
