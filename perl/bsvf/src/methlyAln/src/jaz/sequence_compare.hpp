/***
 *  $Id$
 **
 *  File: sequence_compare.hpp
 *  Created: May 03, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  			Xiao Yang <isuyang@gmail.com>
 *  Copyright (c) 2012-2013 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *
 *  Boost Software License - Version 1.0 - August 17th, 2003
 *
 *  Permission is hereby granted, free of charge, to any person or organization
 *  obtaining a copy of the software and accompanying documentation covered by
 *  this license (the "Software") to use, reproduce, display, distribute,
 *  execute, and transmit the Software, and to prepare derivative works of the
 *  Software, and to permit third-parties to whom the Software is furnished to
 *  do so, all subject to the following:
 *
 *  The copyright notices in the Software and this entire statement, including
 *  the above license grant, this restriction and the following disclaimer,
 *  must be included in all copies of the Software, in whole or in part, and
 *  all derivative works of the Software, unless such copies or derivative
 *  works are solely in the form of machine-executable object code generated by
 *  a source language processor.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 *  SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 *  FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 *  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 *  DEALINGS IN THE SOFTWARE.
 */

#ifndef SEQUENCE_COMPARE_HPP
#define SEQUENCE_COMPARE_HPP

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <climits>
#include <map>
#include <string>
#include <vector>
#include <tuple>


namespace bio {

  namespace detail {

    // this code comes from jaz
    template <typename Iter1, typename Iter2, typename Pred>
    std::size_t intersection_size(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, Pred pred) {
	std::size_t S = 0;

	while ((first1 != last1) && (first2 != last2)) {
	    if (pred(*first1, *first2)) ++first1;
	    else if (pred(*first2, *first1)) ++first2;
	    else {
		first1++;
		first2++;
		S++;
	    }
	} // while

	return S;
    } // intersection_size


    template <typename Iter1, typename Iter2>
    int count_distance(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2) {
	int S = 0;

	while ((first1 != last1) && (first2 != last2)) {
	    if (first1->first < first2->first) {
		S += (first1->second * first1->second);
		++first1;
	    }
	    else if (first2->first < first1->first) {
		S += (first2->second * first2->second);
		++first2;
	    }
	    else {
		int d = (first1->second - first2->second);
		S += d * d;
		first1++;
		first2++;
	    }
	} // while

	return S;
    } // count_distance


    template <typename Sequence> void general_kmer_index(const std::string& s, unsigned int k, Sequence& S) {
	unsigned int l = s.size();
	unsigned int end = l - k + 1;
	S.resize(end);
	for (unsigned int i = 0; i < end; ++i) {
	    S[i] = std::string(s.begin() + i, s.begin() + i + k);
	}
	std::sort(S.begin(), S.end());
    } // general_kmer_index


    template <typename Map> void general_kmer_count(const std::string& s, unsigned int k, Map& S) {
	S.clear();
	unsigned int l = s.size();
	unsigned int end = l - k + 1;
	for (unsigned int i = 0; i < end; ++i) {
	    S[std::string(s.begin() + i, s.begin() + i + k)]++;
	}
    } // general_kmer_count


    class dna_digit {
    public:
	dna_digit() {
	    std::memset(digit_, 0, 256);
	    digit_['c'] = digit_['C'] = 1;
	    digit_['g'] = digit_['G'] = 2;
	    digit_['t'] = digit_['T'] = 3;
	    digit_['u'] = digit_['U'] = 3;
	} // dna_digit

    protected:
	char digit_[256];

    }; // dna_digit


    class dna_kmer_index : public dna_digit {
    public:
	dna_kmer_index() : dna_digit() { }

	template <typename Sequence>
	void operator()(const std::string& s, unsigned int k, Sequence& S) {
	    unsigned int l = s.size();
	    unsigned int end = l - k + 1;

	    S.resize(end);

	    // first kmer
	    unsigned long long int v = digit_[s[k - 1]];
	    for (unsigned int i = 0; i < k - 1; ++i) {
		v += digit_[s[i]] * (1ULL << ((k - i - 1) << 1));
	    }

	    S[0] = v;

	    // and then all other
	    unsigned long long int b = 1ULL << ((k - 1) << 1);

	    for (unsigned int i = 1; i < end; ++i) {
		v = (v - b * digit_[s[i - 1]]) * 4 + digit_[s[i + k - 1]];
		S[i] = v;
	    }

	    std::sort(S.begin(), S.end());
	} // operator()

    }; // class dna_kmer_index


    class dna_kmer_count : public dna_digit {
    public:
	dna_kmer_count() : dna_digit() { }

	template <typename Map>
	void operator()(const std::string& s, unsigned int k, Map& S) {
	    unsigned int l = s.size();
	    unsigned int end = l - k + 1;

	    // first kmer
	    unsigned long long int v = digit_[s[k - 1]];

	    for (unsigned int i = 0; i < k - 1; ++i) {
		v += digit_[s[i]] * (1ULL << ((k - i - 1) << 1));
	    }

	    S[v] = 1;

	    // and then all other
	    unsigned long long int b = 1ULL << ((k - 1) << 1);

	    for (unsigned int i = 1; i < end; ++i) {
		v = (v - b * digit_[s[i - 1]]) * 4 + digit_[s[i + k - 1]];
		S[v]++;
	    }
	} // operator()

    }; // class dna_kmer_count

  } // namespace detail



  template <typename Derived> struct sequence_compare {
      std::tuple<int, int, int> operator()(const std::string& s0, const std::string& s1) {
	  return static_cast<Derived*>(this)->operator()(s0, s1);
      }
  }; // struct sequence_compare


  /** Class: scoring_matrix
   *
   *  Functor encapsulating a scoring_matrix functionality.
   */
  class scoring_matrix {
  public:
      scoring_matrix() : sz_(0), matrix_(0) { std::memset(sigma_, 0, 256); }

      /** Constructor: scoring_matrix
       *
       *  Parameter:
       *  sigma -  Map of the alphabet used by the matrix.
       *  matrix - Raw-wise substitution matrix.
       */
      scoring_matrix(unsigned char sigma[256], const std::vector<char>& matrix)
	  : matrix_(matrix), sz_(std::sqrt(matrix.size())) { std::memcpy(sigma_, sigma, 256); }

      /** Function: operator()
       *
       *  Returns:
       *  Substitution score between a and b.
       */
      int operator()(char a, char b) const { return matrix_[sigma_[a] * sz_ + sigma_[b]]; }


  private:
      unsigned int sz_;
      unsigned char sigma_[256];
      std::vector<char> matrix_;

  }; // scoring_matrix


  inline scoring_matrix make_dummy_sm(int m, int s, int meth, bool strainF) {
      unsigned char sigma[256];
      for (unsigned int i = 0; i < 256; ++i) sigma[i] = i;
      std::vector<char> matrix(256 * 256, s);
      for (unsigned int i = 0; i < 256; ++i) matrix[(i << 8) + i] = m;
	if (meth != 0) {
		unsigned char methlybases[4][5] = {"CTct","Yy","GAga","Rr"};
		//unsigned char methlybases[2][7] = {"CYcy","GRgr"};
		int x;
		if (strainF) {
			x = 0;
		} else {
			x = 1;
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 2; j++) {
				if ( i != j) {
					for (int p = 0; p <= 2; p+=2) {
						matrix[(methlybases[p][i] << 8) + methlybases[p+1][j]] = meth;
						matrix[(methlybases[p+1][j] << 8) + methlybases[p][i]] = meth;
					}
				}
			}
		}
	}
      return scoring_matrix(sigma, matrix);
  } // make_dummy_sm

  inline scoring_matrix make_dna_sm(int m, int s) {
      unsigned char sigma[256];
      for (unsigned int i = 0; i < 256; ++i) sigma[i] = 4;

      sigma['a'] = sigma['A'] = 0;
      sigma['c'] = sigma['C'] = 1;
      sigma['g'] = sigma['G'] = 2;
      sigma['t'] = sigma['T'] = 3;

      std::vector<char> matrix(5 * 5, s);
      for (unsigned int i = 0; i < 5; ++i) matrix[(5 * i) + i] = m;

      return scoring_matrix(sigma, matrix);
  }; // make_dna_sm

  inline bool read_file_sm(const std::string& name, scoring_matrix& sub) {
      std::ifstream f(name.c_str());
      if (!f) return false;

      unsigned char sigma[256];
      std::string buf;

      // read comments
      while (!f.eof()) {
	  buf = "";
	  std::getline(f, buf);
	  if ((buf.empty() == true) || (buf[0] != '#')) break;
      } // while

      if (buf.empty() == true) return false;

      // parse column header
      std::string head = buf;
      head.erase(std::remove(head.begin(), head.end(), ' '), head.end());

      int len = head.size();
      if (head.at (head.length() - 1) != '*') return false;

      std::vector<char> matrix(len * len, 0);

      for (int i = 0; i < 256; ++i) sigma[i] = len - 1;
      for (int i = 0; i < len - 1; ++i) sigma[head[i]] = i;

      // read matrix
      for (unsigned int i = 0; i < len; ++i) {
	  char id = 0;
	  int val;

	  f >> id;
	  if (!f || (id != head[i])) return false;

	  for (unsigned int j = 0; j < len; ++j) {
	      f >> val;
	      if (!f) return false;
	      unsigned int pos0 = sigma[head[i]];
	      unsigned int pos1 = sigma[head[j]];
	      matrix[pos0 * len + pos1] = val;
	  } // for j
      } // for i

      f.close();

      sub = scoring_matrix(sigma, matrix);
      return true;
  } // read_scoring_matrix


  /** Class: global_alignment
   *  Functor implementing global pairwise sequence alignment
   *  with affine gap penalty. The implementation is alphabet oblivious.
   *  No thread safety guarantees.
   *
   *  set_alignment_type() function modifies the global alignment type:
   *
   *  type 0: default normal affine gap global alignment
   *  type 1: suffix-prefix or containment alignment, where end gaps
   *		are penalty free (differing from type 0: matrix initialization
   *		and traceback starts from the max value of last row/col)
   *  type 2: prefix alignment with penalty free gaps in the suffix
   *  (differing from type 0, traceback starts from the max value of the
   *  last row/col)
   *		 e.g. ACTGTTCCC
   *	    		  -CTGTT---
   *	     the first indel is penalized whereas the last 3 indels are not
   */
  class global_alignment : public sequence_compare<global_alignment> {
  public:
      /** Constructor: global_alignment
       *
       *  Parameter:
       *  m - Match score (some positive number).
       *  s - Substitution penalty (usually negative number).
       *  g - Gap opening penalty (negative number).
       *  h - Gap extension penalty (negative number).
       */
      explicit global_alignment(int m = 0, int s = 0, int g = 0, int h = 0, int meth = 0, bool strainF = true)
	  : sub_(make_dummy_sm(m, s, meth, strainF)), g_(g), h_(h) {
    	  	  alignment_type_ = 0;
    	  	  path_prefix_ = path_suffix_ = "";
      }


      /** Constructor: global_alignment
       *
       *  Parameter:
       *  sm - Substitution matrix.
       *  g -  Gap opening penalty (negative number).
       *  h -  Gap extension penalty (negative number).
       */
      global_alignment(const scoring_matrix& sm, int g, int h)
	  : sub_(sm), g_(g), h_(h) {
    	  	  alignment_type_ = 0;
	  	  path_prefix_ = path_suffix_ = "";
      }

      void set_alignment_type (int flag) {
    	  	  alignment_type_ = flag;
      }

      /** Function: operator()
       *  Compute alignment score between s0 and s1.
       *
       *  Returns:
       *  3-tuple (alignment score, alignment length without terminal gaps,
       *  			number of matches).
       */
      std::tuple<int, int, int> operator()(const std::string& s0, const std::string& s1) {
		  int n = s0.size();
		  int m = s1.size();

		  S_.clear(); I_.clear(); D_.clear(); track_.clear();
		  // S(i, j) = max{ I(i, j), D(i, j), S(i - 1, j - 1) + d(i,j) 	}
		  // D(i, j) = max{ D(i, j - 1), S(i, j - 1) + g } + h
		  // I(i, j) = max{ I(i - 1, j), S(i - 1, j) + g } + h
		  S_.resize ((n + 1) * (m + 1));
		  I_.resize ((n + 1) * (m + 1));
		  D_.resize ((n + 1) * (m + 1));
		  track_.resize((n + 1) * (m + 1));

		  // ------ matrix initialization for first row and column ------
		  S_[0] = 0;
		  if (alignment_type_ == 1) I_[0] = D_[0] = 0;
		  else I_[0] = D_[0] = g_;

		  for (int j = 1; j <= m; ++ j) { // first row
			  track_[j] = LEFT;
			  I_[j] = -100;
			  if (alignment_type_ == 1) D_[j] = S_[j] = 0;
			  else S_[j] = D_[j] = g_ + j * h_;
		  }

		  for (int i = 1; i <= n; ++ i) { // first col
			  int idx = (m + 1) * i;
			  track_[idx] = TOP;
			  D_[idx] = -100;
			  if (alignment_type_ == 1) I_[idx] = S_[idx] = 0;
			  else S_[idx] = I_[idx] = g_ + i * h_;
		  }

		  // ----------------- matrix filling --------------------
		  for (int i = 1; i <= n; ++ i) { // row
			  for (int j = 1; j <= m; ++ j) { // col
				 int idx = i * (m + 1) + j;

				 D_[idx] = std::max (D_[idx - 1], S_[idx - 1] + g_) + h_;
				 I_[idx] = std::max (I_[idx - m - 1], S_[idx - m - 1] + g_) + h_;
				 S_[idx] = S_[idx - m - 2] + sub_(s0[i - 1], s1[j - 1]);

				 track_[idx] = DIAG;
				 if (D_[idx] >= S_[idx]) {
					 S_[idx] = D_[idx];
					 track_[idx] = LEFT;
				 }
				 if (I_[idx] >= S_[idx]) {
					 S_[idx] = I_[idx];
					 track_[idx] = TOP;
				 }
			  }
		  }

		  {// debug: print out content of S_ I_ D_
			  for (int i = 0; i <= n ; ++ i) {
				  for (int j = 0; j <= m; ++ j) {
					  int idx = i * (m + 1) + j;
					  std::cout << std::setw(6) << S_[idx];
				  }
				  std::cout << "\n";
				  for (int j = 0; j <= m; ++ j) {
					  int idx = i * (m + 1) + j;
					  std::cout << std::setw(6) << I_[idx];
				  }
				  std::cout << "\n";
				  for (int j = 0; j <= m; ++ j) {
					  int idx = i * (m + 1) + j;
					  std::cout << std::setw(6) << D_[idx];
				  }
				  std::cout << "\n";
				  std::cout << "\n";
			  }
			  // trace matrix
			  for (int i = 0; i <= n; ++ i) {
				  for (int j = 0; j <= m; ++ j) {
					  int idx = i * (m + 1) + j;
					  int track = track_[idx];
					  switch (track){
					  case 0: std::cout << "Top\t";
					  break;
					  case 1: std::cout << "lef\t";
					  break;
					  case 2: std::cout << "dia\t";
					  break;
					  }
				  }
				  std::cout << "\n";
			  }
		  }

		  // ----------------- backtracking --------------------
		  int i = n;
		  int j = m;
		  int max_val = S_.back();
		  if (alignment_type_ != 0) {
			  // identify the idx of the max value of last row
			  max_val = INT_MIN;
			  int row = n, col = 1; // row col wrt to (n+1)x(m+1) matrix
			  for (; col <= m; ++ col) { // last row
				  int idx = row * (m + 1) + col;
				  if (S_[idx] > max_val) {
					  max_val = S_[idx];
					  i = row;
					  j = col;
				  }
			  }
			  // identify the idx of the max value of last col
			  col = m; row = 1;
			  for (; row <= n; ++ row) {
				  int idx = row * (m + 1) + col;
				  if (S_[idx] > max_val) {
					  max_val = S_[idx];
					  i = row;
					  j = col;
				  }
			  }

			  if (i < n) path_prefix_ = std::string (n - i, 'D');
			  else if (j < m) path_prefix_ = std::string (m - j, 'I');

		  }

		  /*{ // debug print
		  std::cout << "path_prefix = " << path_prefix_ << "\n";
		  std::cout << "i, j: " << i << ", " << j << "\n";
		  }*/

		  int match = 0;
		  int length = 0;

		  bool has_gap = false;
		  int sgap = 0;
		  int egap = 0;

		  has_path_ = false;
		  path_.clear();
		  overlap_ = 0;

		  while ((i > 0) || (j > 0)) {

			  std::cout << "i, j = " << i << ", " << j << "\n";

			  switch (track_[i * (m + 1) + j]) {
				case TOP:
					--i;
					sgap++;
					path_.push_back('D');
					break;

				case LEFT:
					--j;
					sgap++;
					path_.push_back('I');
					break;

				case DIAG:

					if (s0[i - 1] == s1[j - 1]) {
						match++;
						path_.push_back('M');
					} else if ( abs(s0[i - 1] - s1[j - 1])==32 ) {
						path_.push_back('m');
					} else path_.push_back('R');

					++overlap_;
					--i;
					--j;

					if (has_gap == false) {
						has_gap = true;
						egap = sgap;
					}

					sgap = 0;
					break;
			  } // switch

			  length ++;
		  } // while

		  {
			  std::cout << "i, j = " << i << "," << j << "\n";
		  }

		  //if (i > 0) path_suffix_ = std::string (i, 'D');
		  //else if (j > 0) path_suffix_ = std::string (j, 'I');

		  path_ = path_prefix_ + path_;
		  //path_ += path_suffix_;

		  return std::make_tuple(max_val, length - sgap - egap, match);
      } // operator()

      /** Function: path
       *  Backtrack the last computed alignment to obtain its edit path.
       *
       *  Returns:
       *  Edit path where 'I' means insert wrt s0, 'D' is deletion, 'R' is substitution,
       *  and 'M' is match.
       */
      std::string path() {
		  if (has_path_ == false) {
			  std::reverse(path_.begin(), path_.end());
			  has_path_ = true;
		  }
		  char numstr[22]; // [21] is enough to hold all numbers up to 64-bits
		  sprintf(numstr, ",%d", overlap_);
		  return path_ + numstr;
      } // path


  private:
      enum { TOP, LEFT, DIAG };

      int alignment_type_; // default 0: standard global alignment
      	  	  	  	  	  // 		 1: end gap free alignment,
        					  //            suffix-prefix/containment etc
      	  	  	  	  	  //         2: suffix gap free alignment,
      	  	  	  	  	  //            penalize gaps in the prefix
      bool has_path_;
	  int overlap_;
      std::string path_;
      std::string path_prefix_, // if alignment end within one string,
      	  	 path_suffix_;  	   // these store what needs to be attached,
      	  	  	  	  	  	   // in prefix or suffix e.g. DDDD or III
      std::vector<unsigned char> track_;

      std::vector<int> S_;
      std::vector<int> I_;
      std::vector<int> D_;
      scoring_matrix sub_;

      int g_;
      int h_;

  }; // class global_alignment


} // namespace bio

#endif // SEQUENCE_COMPARE_HPP
