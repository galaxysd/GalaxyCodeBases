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

// This is the main program for the decompressor.

#include "bwt_compress.h"
#include "stream.h"
#include "inttypes.h"

#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>

namespace dcsbwt {
int verbosity;
int statistics;
}

using dcsbwt::verbosity;

void  Fail(const char* message, const char* argument) {
  if (message)
    std::cerr << message << argument << "\n";
  std::cerr << "Usage: dcsunbwt [options] inputfile outputfile\n"
            << "  --transform=args\n"
            << "  --memory=size[k|K|m|M|g|G]\n"
            << "  --verbosity=level\n";
  std::exit(EXIT_FAILURE);
}

void  Fail(const char* message, int64 argument) {
  char buffer[15];
  snprintf(buffer, 15, "%lld", argument);
  Fail(message, buffer);
}

class CommandLineOptions {
 public:
  CommandLineOptions()
      : transform_("f"),
        memory_(900LL * (1LL<<20)),
        verbosity_(0) {}

  std::string GetTransform() const { return transform_; }
  long long int GetMemory() const { return memory_; }
  int GetVerbosity() const { return verbosity_; }
  int GetStatistics() const { return 0; }

  int ParseOptions(int argc, char** argv) {
    while (1) {
      static struct option long_options[] = {
        {"transform", 1, 0, 't'},
        {"memory", 1, 0, 'm'},
        {"verbosity", 1, 0, 'v'},
        {0, 0, 0, 0}
      };

      int option_index = 0;
      int c = getopt_long(argc, argv, "t:m:v:",
                          long_options, &option_index);

      if (-1 == c) break;

      char* next = NULL;
      long long int size = 0;
      int level;
      switch (c) {
        case 't':
          transform_ = optarg;
          break;
        case 'm':
          size = ParseSize(optarg);
          if (size < 0)
            Fail("invalid option argument: --memory=", optarg);
          memory_ = size;
          break;
        case 'v':
          level = ParseLevel(optarg);
          if (level == kInvalidLevel)
            Fail("invalid option argument: --verbosity=", optarg);
          verbosity_ = level;
          break;
        case '?':
          Fail(NULL, "");
        default:
          std::cerr << "Unexpected error: getopt returned: " << c << '\n';
          std::exit(EXIT_FAILURE);
      }
    }
    return optind;
  }

 private:
  std::string transform_;
  long long int memory_;
  int verbosity_;

  long long int ParseSize(const char* argument) {
    char* next;
    long long int size = strtoll(argument, &next, 0);
    if (*next) {
      switch (*next++) {
        case 'k':
        case 'K':
          size *= (1LL<<10);
          break;
        case 'm':
        case 'M':
          size *= (1LL<<20);
          break;
        case 'g':
        case 'G':
          size *= (1LL<<30);
          break;
        default:
          size = -1;
      }
      if (*next) size = -1;
    }
    return size;
  }

  static const int kInvalidLevel = INT_MIN;

  int ParseLevel(const char* argument) {
    char* next;
    long long int level = strtoll(argument, &next, 0);
    if (*next || level < INT_MIN || level > INT_MAX)
      level = kInvalidLevel;
    return level;
  }
};


int main(int argc, char** argv) {

  // Memory for top level code, file buffers, stream buffers, etc.
  // This is currently a rather arbitrary value.
  static const int64 kMemoryOverhead = (1LL << 20);
  static const int64 kMaxMemory = (3LL << 30);
  static const int kInputBufferSize = 1 << 16;

  // Process command line options
  CommandLineOptions command_line_options;
  int first_argument = command_line_options.ParseOptions(argc, argv);
  verbosity = command_line_options.GetVerbosity();
  dcsbwt::statistics = command_line_options.GetStatistics();

  if (argc - first_argument < 2) Fail("Too few arguments", "");
  if (argc - first_argument > 2) Fail("Too many arguments", "");

  dcsbwt::BwtDecompressor::Options decompressor_options;
  bool ok = decompressor_options.SetTransformOptions(
      command_line_options.GetTransform());
  if (!ok) Fail("Invalid option argument: --transform=",
                 command_line_options.GetTransform().c_str());

  int64 memory_budget = command_line_options.GetMemory();
  int64 available_memory = memory_budget - kMemoryOverhead;
  ok = decompressor_options.SetMemoryBudget(available_memory);
  if (!ok)
    Fail("Invalid (possibly too small) memory budget: ", memory_budget);

  // Open files and setup streams
  const char* infilename = argv[first_argument];
  if (verbosity > 0)
    std::cout << "Opening input file: " << infilename
              << std::endl;
  std::FILE* infile = std::fopen(infilename, "r");
  if (!infile) {
    std::cerr << "Cannot open input file " << infilename << '\n';
    std::exit(EXIT_FAILURE);
  }
  dcsbwt::InStreamFromFile instream(infile);
  dcsbwt::InStreamBuffer inbuffer(kInputBufferSize);
  inbuffer.Connect(&instream);

  const char* outfilename = argv[first_argument + 1];
  if (verbosity > 0) {
    std::clog << "Opening output file: " << outfilename << std::endl;
  }
  std::FILE* outfile = std::fopen(outfilename, "w");
  if (!outfile) {
    std::cerr << "Cannot open output file " << outfilename << '\n';
    std::fclose(infile);
    std::exit(EXIT_FAILURE);
  }
  dcsbwt::OutStreamToFile outstream(outfile);

  // Decompress
  int64 result = dcsbwt::BwtDecompressor::Decompress(
      &inbuffer, &outstream, decompressor_options);
  if (result < 0) {
    std::cerr << "Corrupted or non-BWT-compressed file: " << infilename
              << '\n';
    std::fclose(outfile);
    inbuffer.Disconnect();
    std::fclose(infile);
    std::exit(EXIT_FAILURE);
  } else if (result > 0) {
    std::cerr << "Cannot decompress a block of size " << result
              << " under the memory budget of "
              << memory_budget << '\n';
    std::fclose(outfile);
    inbuffer.Disconnect();
    std::fclose(infile);
    std::exit(EXIT_FAILURE);
  }

  // Close files and streams
  std::fclose(outfile);
  inbuffer.Disconnect();
  std::fclose(infile);

  return 0;
}
