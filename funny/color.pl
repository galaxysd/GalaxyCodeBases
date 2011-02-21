#!/bin/env perl
use strict;
use warnings;

for my $bg (0..7) {
	for my $fg (0..7) {
		print "\033[",$fg+30,";",$bg+40,"m","[t$fg,${bg}t]","\033[0m";
	}
	print "\n";
}

for my $attr (0,1,5,7,4,8,9) {
	print "\nAttr:$attr\n";
	for my $bg (0..7) {
		print "before ";
		for my $fg (0..7) {
			print "\033[",$attr,";",$fg+30,";",$bg+40,"m","[t$fg,${bg}t]","\033[0m";
		}
		print " after\n";
	}
}
