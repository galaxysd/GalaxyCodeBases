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

// This is the main program for the compressor.

#include "bwt_compress.h"
#include "stream.h"
#include "inttypes.h"

#include <string>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <getopt.h>

namespace dcsbwt {
int verbosity;
int statistics;
}

using dcsbwt::verbosity;

void  Fail(const char* message, const char* argument) {
  if (message)
    std::cerr << message << argument << "\n";
  std::cerr << "Usage: dcsbwt [options] inputfile outputfile\n"
            << "  --transform=args\n"
            << "  --compression=args\n"
            << "  --memory=size[k|K|m|M|g|G]\n"
            << "  --blocksize=size[k|K|m|M|g|G]\n"
            << "  --verbosity=level\n"
            << "  --statistics=level\n";
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
      : transform_("d"),
        compression_("r"),
        memory_(900LL * (1LL<<20)),
        blocksize_(1LL<<20),
        verbosity_(0),
        statistics_(0){}

  std::string GetTransform() const { return transform_; }
  std::string GetCompression() const { return compression_; }
  long long int GetMemory() const { return memory_; }
  long long int GetBlocksize() const { return blocksize_; }
  int GetVerbosity() const { return verbosity_; }
  int GetStatistics() const { return statistics_; }

  int ParseOptions(int argc, char** argv) {
    while (1) {
      static struct option long_options[] = {
        {"transform", 1, 0, 't'},
        {"compression", 1, 0, 'c'},
        {"memory", 1, 0, 'm'},
        {"blocksize", 1, 0, 'b'},
        {"verbosity", 1, 0, 'v'},
        {"statistics", 1, 0, 's'},
        {0, 0, 0, 0}
      };

      int option_index = 0;
      int c = getopt_long(argc, argv, "t:c:m:b:v:s:",
                          long_options, &option_index);

      if (-1 == c) break;

      char* next = NULL;
      long long int size = 0;
      int level = 0;
      switch (c) {
        case 't':
          transform_ = optarg;
          break;
        case 'c':
          compression_ = optarg;
          break;
        case 'm':
          size = ParseSize(optarg);
          if (size < 0)
            Fail("invalid option argument: --memory=", optarg);
          memory_ = size;
          break;
        case 'b':
          size = ParseSize(optarg);
          if (size < 0)
            Fail("invalid option argument: --blocksize=", optarg);
          blocksize_ = size;
          break;
        case 'v':
          level = ParseLevel(optarg);
          if (level == kInvalidLevel)
            Fail("invalid option argument: --verbosity=", optarg);
          verbosity_ = level;
          break;
        case 's':
          level = ParseLevel(optarg);
          if (level == kInvalidLevel)
            Fail("invalid option argument: --statistics=", optarg);
          statistics_ = level;
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
  std::string compression_;
  long long int memory_;
  long long int blocksize_;
  int verbosity_;
  int statistics_;

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

  // Process command line options
  CommandLineOptions command_line_options;
  int first_argument = command_line_options.ParseOptions(argc, argv);
  verbosity = command_line_options.GetVerbosity();
  dcsbwt::statistics = command_line_options.GetStatistics();

  if (argc - first_argument < 2) Fail("Too few arguments", "");
  if (argc - first_argument > 2) Fail("Too many arguments", "");

  dcsbwt::BwtCompressor::Options compressor_options;
  bool ok = compressor_options.SetTransformOptions(
      command_line_options.GetTransform());
  if (!ok) Fail("Invalid option argument: --transform=",
                command_line_options.GetTransform().c_str());
  ok = compressor_options.SetCompressionOptions(
      command_line_options.GetCompression());
  if (!ok) Fail("Invalid option argument: --compression=",
                command_line_options.GetCompression().c_str());

  int64 memory_budget = command_line_options.GetMemory();
  int64 available_memory = memory_budget - kMemoryOverhead;
  ok = compressor_options.SetMemoryBudget(available_memory);
  if (!ok)
    Fail("Invalid (possibly too small) memory budget: ", memory_budget);
  if (verbosity > 0) {
    std::clog << "Compressing with memory budget of " << memory_budget
              << " bytes (" << available_memory << " after general overhead)"
              << std::endl;;
  }

  int64 blocksize = command_line_options.GetBlocksize();
  if (0 == blocksize) {
    if (verbosity > 0)
      std::clog << "No block size given: using suggested block size"
                << " (no larger than "
                << compressor_options.SuggestedBlockSize() << ")"
                << std::endl;
  } else {
    if (blocksize > compressor_options.MaxBlockSize())
      Fail("Too large blocksize", blocksize);
    if (verbosity > 0)
      std::clog << "Using block size " << blocksize
                << " (max=" << compressor_options.MaxBlockSize() << ")"
                << std::endl;
  }

  // Open files and set up streams
  const char* infilename = argv[first_argument];
  if (verbosity > 0)
    std::cout << "Opening input file: " << infilename
              << std::endl;
  std::FILE* infile = std::fopen(infilename, "r");
  if (!infile) {
    std::cerr << "Cannot open input file " << infilename << '\n';
    std::exit(EXIT_FAILURE);
  }

  std::fseek(infile, 0, SEEK_END);
  int64 inputsize = std::ftell(infile);
  std::rewind(infile);
  if (verbosity > 0) {
    std::clog << "Input: " << inputsize << " bytes from file "<< infilename
              << std::endl;
  }
  dcsbwt::InStreamFromFile instream(infile);

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

  // Compress
  dcsbwt::BwtCompressor::Compress(&instream, inputsize, &outstream,
                                  compressor_options, blocksize);
  // Close files
  std::fclose(outfile);
  std::fclose(infile);

  return 0;
}
