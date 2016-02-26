#!/usr/bin/env littler

source("a1funs.r")

#a1data <- as.matrix(read.table("http://www.cs.utoronto.ca/~radford/csc120/a1data.txt"))
a1data <- as.matrix(read.table('a1data.txt'))

DEBUG <- T

set.seed(1)
find_pairings(a1data,1,TRUE)
set.seed(1)
find_pairings(a1data,10,TRUE)

set.seed(999612806)
find_pairings(a1data,1,TRUE)
set.seed(999612806)
find_pairings(a1data,10,TRUE)
