From <https://github.com/boris-dimitrov/dedup>.

# dedup
Help identify duplicate files that can be profitably replaced with hardlinks.

A dangerous bandaid for sorely missing Mac OS functionality.


# caution

Replacing duplicates with hardlinks is generally unsafe.

It can be particularly harmful to

  * files under version control
  
  * files that are subject to automatic updates
  
  * files whose attributes (e.g., last access time) are
just as important as their contents


# usage

Each application requires truly careful review and manual
adjustments to `dedup.py` code to minimize the possiblity of
harming important files.  Use at your own risk.
