## Natural sort in C/C++

<https://www.o-rho.com/naturalsort>

Default sorting is often ASCII sort, a sort based on comparing characters on their position in the ASCII table. This sort produces results awkward for humans when the strings to compare contain numbers, as in the ASCII sort of

````
file1.txt
file10.txt
file3.txt
````

How to implement a natural sort is a popular question on the web. Many solutions are presented, quite a few slow, ugly or not fully correct. To be sure, what actually defines a natural sort is a partly matter of taste (or requirements), but here we define a string to be interpreted as an alternating sequence of all-alpha and all-digit tokens (each tokens one or more chars, so numbers are limited to natural integer values) that can be compared either numerically or alphabetically to other strings by matching those tokens in-order;

````
file 1 .txt
file 10 .txt
file 3 .txt
````

As the token 'file' is the same in all strings, the second token is an all-digit token and numerically compared to the corresponding other all-digit tokens in the other strings, giving

````
file1.txt
file3.txt
file10.txt
````

whereas the algorithm produces

````
file1.txt
file.txt
filename.txt
````

as shorter all-alpha tokens < longer all-alpha tokens, and all-digit tokens < all-alpha tokens.

Quite a few algorithms presented do not produce correct results as they use string-to-number conversions to compare the numbers, causing overflow errors when the string numbers are too large to be represented with numerical datatypes. A more efficient method to compare the 'numbers' is possible by applying two basic truths about such number strings;

 * if two number strings do not have leading zeroes, the longer string is the bigger number.
 * if two number strings do not have leading zeroes and have the same length, an regular ASCII compare yields the same result as a numerical compare.
A simple state machine divides and conquers the problem.

Behavior by example:

````
$ ./natsort z a 9 0 10 99 a00a00 a000a0 a0a000 0a00 00a0 01a01 1a1a11 z1.txt z10.txt z99.txt 3.10.1 3.9.4 z3.txt "" '!' '_' a0a000000 a00a0

0
0a00
00a0
01a01
1a1a11
3.9.4
3.10.1
9
10
99
!
_
a
a0a000
a0a000000
a00a0
a00a00
a000a0
z
z1.txt
z3.txt
z10.txt
z99.txt
````

The algorithm, in STL less-than style

````C++
/**
 * STL natural less-than string compare
 * @return true when natural s1 < s2
 */
bool natstrlt( const char* s1, const char* s2 ) {
  const char* p1 = s1;
  const char* p2 = s2;
  const unsigned short st_scan = 0;
  const unsigned short st_alpha = 1;
  const unsigned short st_numeric = 2;
  unsigned short state = st_scan;
  const char* numstart1 = 0;
  const char* numstart2 = 0;
  const char* numend1 = 0;
  const char* numend2 = 0;
  unsigned long sz1 = 0;
  unsigned long sz2 = 0;
  while ( *p1 && *p2  ) {
    switch ( state ) {
      case st_scan:
        if ( !t_digit(*p1) && !t_digit(*p2) ) {
          state = st_alpha;
          if ( *p1 == *p2 ) {
            p1++;p2++;
          } else return *p1 < *p2;
        } else
        if ( t_digit(*p1) && !t_digit(*p2) ) return true;
        else if ( !t_digit(*p1) && t_digit(*p2) ) return false;
        else {
          state = st_numeric;
          if ( sz1 == 0 )
            while ( *p1 == '0' ) {p1++; sz1++;}
          else
            while ( *p1 == '0' ) p1++;
          if ( sz2 == 0 )
            while ( *p2 == '0' ) {p2++; sz2++;}
          else
            while ( *p2 == '0' ) p2++;
          if ( sz1 == sz2 ) { sz1 = 0; sz2 = 0; };
          if ( !t_digit(*p1) ) p1--;
          if ( !t_digit(*p2) ) p2--;
          numstart1 = p1;
          numstart2 = p2;
          numend1 = numstart1;
          numend2 = numstart2;
        }
        break;
      case st_alpha:
        if ( !t_digit(*p1) && !t_digit(*p2) ) {
          if ( *p1 == *p2 ) {
            p1++;p2++;
          } else return *p1 < *p2;
        } else state = st_scan;
        break;
      case st_numeric:
        while ( t_digit(*p1) ) numend1 = p1++;
        while ( t_digit(*p2) ) numend2 = p2++;
        if ( numend1-numstart1 == numend2-numstart2 &&
            !strncmp( numstart1,numstart2,numend2-numstart2+1) ) state = st_scan; else {
          if ( numend1-numstart1 != numend2-numstart2 )
            return numend1-numstart1 < numend2-numstart2;
          while ( *numstart1 && *numstart2 ) {
            if ( *numstart1 != *numstart2 ) return *numstart1 < *numstart2;
            numstart1++;
            numstart2++;
          }
        }
        break;
    }
  }
  if ( sz1 < sz2 ) return true;
  if ( sz1 > sz2 ) return false;
  if ( *p1 == 0 && *p2 != 0 ) return true;
  if ( *p1 != 0 && *p2 == 0 ) return false;  
  return false;
}
````

See source [header](https://www.o-rho.com/sites/default/files/article_files/natsort.hpp) and [implementation](https://www.o-rho.com/sites/default/files/article_files/natsort.cpp), including C int-returning variant and std::string wrapper.

Tags: Computer science, C/C++ programming

---

## Sorting for Humans : Natural Sort Order

<https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/>

It isn't called "Alphabetical sort"; it's collectively known as **natural sort**. But she's right about one thing: it's hard to find information on natural sorting, and many programmers are completely ignorant of it. None of the common computer languages (that I know of) implement anything other than ASCIIbetical sorts. There are a few places you can find natural sort algorithms, however:

 * Dave Koelle's [The Alphanum Algorithm](http://www.davekoelle.com/alphanum.html)
 * Martin Pool's [Natural Order String Comparison](http://sourcefrog.net/projects/natsort/)
 * Ian Griffiths' [Natural Sorting in C#](http://www.interact-sw.co.uk/iangblog/2007/12/13/natural-sorting)
 * Ned Batchelder's [Compact Python Human Sort](http://nedbatchelder.com/blog/200712.html#e20071211T054956), along with Jussi Salmela's [internationalized version](http://personal.inet.fi/cool/operator/Human%20Sort.py) of same.

Don't let Ned's clever Python ten-liner fool you. Implementing a natural sort is more complex than it seems, and not just for [the gnarly i18n](http://www.codinghorror.com/blog/archives/000813.html) issues I've hinted at, above. But the Python implementations are impressively succinct. One of Ned's commenters posted this version, which is even shorter:

````python
import re
def sort_nicely( l ):
""" Sort the given list in the way that humans expect.
"""
convert = lambda text: int(text) if text.isdigit() else text
alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
l.sort( key=alphanum_key )
````

I tried to come up with a clever, similarly succinct C# 3.0 natural sort implementation, but I failed. I'm not interested in a one-liner contest, necessarily, but it does seem to me that a basic natural sort shouldn't require the 40+ lines of code it takes in most languages.

As programmers, we'd do well to keep Kate's lesson in mind: **ASCIIbetical does not equal alphabetical**. ASCII sorting serves the needs of the computer and the compiler, but what about us human beings? **Perhaps a more human-friendly natural sort option should be built into mainstream programming languages**, too.

---

## The Alphanum Algorithm

<http://www.davekoelle.com/alphanum.html>

### How does it work?

The algorithm breaks strings into chunks, where a chunk contains either all alphabetic characters, or all numeric characters. These chunks are then compared against each other. If both chunks contain numbers, a numerical comparison is used. If either chunk contains characters, the ASCII comparison is used.

There is currently a glitch when it comes to periods/decimal points - specifically, periods are treated only as strings, not as decimal points. The solution to this glitch is to recognize a period surrounded by digits as a decimal point, and continue creating a numeric chunck that includes the decimal point. If a letter exists on either side of the period, or if the period is the first or last character in the string, it should be viewed as an actual period and included in an alphabetic chunk. While I have recently figured this out in theory, I have not yet implemented it into the algorithms. To be truly international, the solution shouldn't just consider periods, but should consider whatever decimal separator is used in the current language.

Currently, the algorithm isn't designed to work with negative signs or numbers expressed in scientific notation, like "5*10e-2". In this case, there are 5 chunks: 5, *, 10, e-, and 2.

The latest version of some of the code (particularly the Java version) compares numbers one at a time if those numbers are in chunks of the same size. For example, when comparing abc123 to abc184, 123 and 184 are the same size, so their values are compared digit-by-digit: 1=1, 2<8. This was done to solve the problem of numeric chunks that are too large to fit in range of values allowed by the programming language for a particular datatype: in Java, an int is limited to 2147483647. The problem with this approach is doesn't properly handle numbers that have leading zeros. For example, 0001 is seem as larger than 1 because it's the longer number. A version that does not compare leading zeros is forthcoming.

### PHP

Just use `sort(&array, SORT_STRING);`

### AlphanumComparator.java

````Java
import java.util.Comparator;

/**
 * This is an updated version with enhancements made by Daniel Migowski,
 * Andre Bogus, and David Koelle
 *
 * To convert to use Templates (Java 1.5+):
 *   - Change "implements Comparator" to "implements Comparator<String>"
 *   - Change "compare(Object o1, Object o2)" to "compare(String s1, String s2)"
 *   - Remove the type checking and casting in compare().
 *
 * To use this class:
 *   Use the static "sort" method from the java.util.Collections class:
 *   Collections.sort(your list, new AlphanumComparator());
 */
public class AlphanumComparator implements Comparator
{
    private final boolean isDigit(char ch)
    {
        return ch >= 48 && ch <= 57;
    }

    /** Length of string is passed in for improved efficiency (only need to calculate it once) **/
    private final String getChunk(String s, int slength, int marker)
    {
        StringBuilder chunk = new StringBuilder();
        char c = s.charAt(marker);
        chunk.append(c);
        marker++;
        if (isDigit(c))
        {
            while (marker < slength)
            {
                c = s.charAt(marker);
                if (!isDigit(c))
                    break;
                chunk.append(c);
                marker++;
            }
        } else
        {
            while (marker < slength)
            {
                c = s.charAt(marker);
                if (isDigit(c))
                    break;
                chunk.append(c);
                marker++;
            }
        }
        return chunk.toString();
    }

    public int compare(Object o1, Object o2)
    {
        if (!(o1 instanceof String) || !(o2 instanceof String))
        {
            return 0;
        }
        String s1 = (String)o1;
        String s2 = (String)o2;

        int thisMarker = 0;
        int thatMarker = 0;
        int s1Length = s1.length();
        int s2Length = s2.length();

        while (thisMarker < s1Length && thatMarker < s2Length)
        {
            String thisChunk = getChunk(s1, s1Length, thisMarker);
            thisMarker += thisChunk.length();

            String thatChunk = getChunk(s2, s2Length, thatMarker);
            thatMarker += thatChunk.length();

            // If both chunks contain numeric characters, sort them numerically
            int result = 0;
            if (isDigit(thisChunk.charAt(0)) && isDigit(thatChunk.charAt(0)))
            {
                // Simple chunk comparison by length.
                int thisChunkLength = thisChunk.length();
                result = thisChunkLength - thatChunk.length();
                // If equal, the first different number counts
                if (result == 0)
                {
                    for (int i = 0; i < thisChunkLength; i++)
                    {
                        result = thisChunk.charAt(i) - thatChunk.charAt(i);
                        if (result != 0)
                        {
                            return result;
                        }
                    }
                }
            } else
            {
                result = thisChunk.compareTo(thatChunk);
            }

            if (result != 0)
                return result;
        }

        return s1Length - s2Length;
    }
}
````

### alphanum.pl

````Perl
#
# TODO: Make decimal points be considered in the same class as digits
#

# usage:
#my @sorted = sort { alphanum($a,$b) } @strings;

sub alphanum {
  # split strings into chunks
  my @a = chunkify($_[0]);
  my @b = chunkify($_[1]);
  
  # while we have chunks to compare.
  while (@a && @b) {
    my $a_chunk = shift @a;
    my $b_chunk = shift @b;
    
    my $test =
        (($a_chunk =~ /\d/) && ($b_chunk =~ /\d/)) ? # if both are numeric
            $a_chunk <=> $b_chunk : # compare as numbers
            $a_chunk cmp $b_chunk ; # else compare as strings  
    
    # return comparison if not equal.
    return $test if $test != 0;
  }

  # return longer string.
  return @a <=> @b;
}

# split on numeric/non-numeric transitions
sub chunkify {
  my @chunks = split m{ # split on
    (?= # zero width
      (?<=\D)\d | # digit preceded by a non-digit OR
      (?<=\d)\D # non-digit preceded by a digit
    )
  }x, $_[0];
  return @chunks;
}
````

### alphanum.py_v2.4

````Python
test_strings = [ "1000X Radonius Maximus", "10X Radonius", "200X Radonius", "20X Radonius", "20X Radonius Prime", "30X Radonius", "40X Radonius", "Allegia 50 Clasteron", "Allegia 500 Clasteron", "Allegia 51 Clasteron", "Allegia 51B Clasteron", "Allegia 52 Clasteron", "Allegia 60 Clasteron", "Alpha 100", "Alpha 2", "Alpha 200", "Alpha 2A", "Alpha 2A-8000", "Alpha 2A-900", "Callisto Morphamax", "Callisto Morphamax 500", "Callisto Morphamax 5000", "Callisto Morphamax 600", "Callisto Morphamax 700", "Callisto Morphamax 7000", "Callisto Morphamax 7000 SE", "Callisto Morphamax 7000 SE2", "QRS-60 Intrinsia Machine", "QRS-60F Intrinsia Machine", "QRS-62 Intrinsia Machine", "QRS-62F Intrinsia Machine", "Xiph Xlater 10000", "Xiph Xlater 2000", "Xiph Xlater 300", "Xiph Xlater 40", "Xiph Xlater 5", "Xiph Xlater 50", "Xiph Xlater 500", "Xiph Xlater 5000", "Xiph Xlater 58" ]

import re

re_chunk = re.compile("([\D]+|[\d]+)")
re_letters = re.compile("\D+")
re_numbers = re.compile("\d+")

def getchunk(item):
	itemchunk = re_chunk.match(item)

	# Subtract the matched portion from the original string
	# if there was a match, otherwise set it to ""
	item = (item[itemchunk.end():] if itemchunk else "")
	# Don't return the match object, just the text
	itemchunk = (itemchunk.group() if itemchunk else "")

	return (itemchunk, item)

def alphanum(a, b):
	n = 0

	while (n == 0):
		# Get a chunk and the original string with the chunk subtracted
		(ac, a) = getchunk(a)
		(bc, b) = getchunk(b)

		# Both items contain only letters
		if (re_letters.match(ac) and re_letters.match(bc)):
			n = cmp(ac, bc)
		else:
			# Both items contain only numbers
			if (re_numbers.match(ac) and re_numbers.match(bc)):
				n = cmp(int(ac), int(bc))
			# One item has letters and one item has numbers, or one item is empty
			else:
				n = cmp(ac, bc)

				# Prevent deadlocks
				if (n == 0):
					n = 1

	return n

test_strings.sort(cmp=alphanum)

for (v) in test_strings:
	print v
````

### alphanum.js

````JavaScript
/* ********************************************************************
 * Alphanum Array prototype version
 *  - Much faster than the sort() function version
 *  - Ability to specify case sensitivity at runtime is a bonus
 *
 */
Array.prototype.alphanumSort = function(caseInsensitive) {
  for (var z = 0, t; t = this[z]; z++) {
    this[z] = new Array();
    var x = 0, y = -1, n = 0, i, j;

    while (i = (j = t.charAt(x++)).charCodeAt(0)) {
      var m = (i == 46 || (i >=48 && i <= 57));
      if (m !== n) {
        this[z][++y] = "";
        n = m;
      }
      this[z][y] += j;
    }
  }

  this.sort(function(a, b) {
    for (var x = 0, aa, bb; (aa = a[x]) && (bb = b[x]); x++) {
      if (caseInsensitive) {
        aa = aa.toLowerCase();
        bb = bb.toLowerCase();
      }
      if (aa !== bb) {
        var c = Number(aa), d = Number(bb);
        if (c == aa && d == bb) {
          return c - d;
        } else return (aa > bb) ? 1 : -1;
      }
    }
    return a.length - b.length;
  });

  for (var z = 0; z < this.length; z++)
    this[z] = this[z].join("");
}


/* ********************************************************************
 * Alphanum sort() function version - case sensitive
 *  - Slower, but easier to modify for arrays of objects which contain
 *    string properties
 *
 */
function alphanum(a, b) {
  function chunkify(t) {
    var tz = new Array();
    var x = 0, y = -1, n = 0, i, j;

    while (i = (j = t.charAt(x++)).charCodeAt(0)) {
      var m = (i == 46 || (i >=48 && i <= 57));
      if (m !== n) {
        tz[++y] = "";
        n = m;
      }
      tz[y] += j;
    }
    return tz;
  }

  var aa = chunkify(a);
  var bb = chunkify(b);

  for (x = 0; aa[x] && bb[x]; x++) {
    if (aa[x] !== bb[x]) {
      var c = Number(aa[x]), d = Number(bb[x]);
      if (c == aa[x] && d == bb[x]) {
        return c - d;
      } else return (aa[x] > bb[x]) ? 1 : -1;
    }
  }
  return aa.length - bb.length;
}


/* ********************************************************************
 * Alphanum sort() function version - case insensitive
 *  - Slower, but easier to modify for arrays of objects which contain
 *    string properties
 *
 */
function alphanumCase(a, b) {
  function chunkify(t) {
    var tz = new Array();
    var x = 0, y = -1, n = 0, i, j;

    while (i = (j = t.charAt(x++)).charCodeAt(0)) {
      var m = (i == 46 || (i >=48 && i <= 57));
      if (m !== n) {
        tz[++y] = "";
        n = m;
      }
      tz[y] += j;
    }
    return tz;
  }

  var aa = chunkify(a.toLowerCase());
  var bb = chunkify(b.toLowerCase());

  for (x = 0; aa[x] && bb[x]; x++) {
    if (aa[x] !== bb[x]) {
      var c = Number(aa[x]), d = Number(bb[x]);
      if (c == aa[x] && d == bb[x]) {
        return c - d;
      } else return (aa[x] > bb[x]) ? 1 : -1;
    }
  }
  return aa.length - bb.length;
}
````
