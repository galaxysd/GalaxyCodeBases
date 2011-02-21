#!/bin/env perl
use strict;
use warnings;

for my $fg (0..9) {
	print "\033[",$fg,"m","[t${fg}t]","\033[0m";
}
print "\n";

for my $fg (0..9) {
	print "\033[",$fg+30,"m","[t${fg}t]","\033[0m";
}
print "\n";

for my $bg (0..9) {
	print "\033[",$bg+40,"m","[t${bg}t]","\033[0m";
}
print "\n";

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

__END__
Attr:
0	Reset
1	FG BRIGHT
5	BG BRIGHT under sh/putty
7	REVERSE FG and BG
4	UNDERLINE under gnome-terminal
8	HIDDENT (sh => fg black, gnome-terminal => fg==bg)
9	Strikethrough under gnome-terminal

http://en.wikipedia.org/wiki/ANSI_escape_code

#include <stdio.h>

#define RESET 		0
#define BRIGHT		1
#define DIM		2
#define UNDERLINE	3
#define BLINK		4
#define REVERSE		7
#define HIDDENT		8

#define BLACK		0
#define RED		1
#define GREEN		2
#define YELLOW		3
#define BLUE		4
#define MAGENTA		5	// Æ·ºì
#define CYAN		6
#define WHITE		7

void clear_screen()
{
	system("clear");
}

void color(int attr, int fg, int bg)
{
	printf("%c[%d;%d;%dm", 0x1B, attr, fg+30, bg+40);
}
