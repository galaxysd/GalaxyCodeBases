#!/usr/bin/env perl
use strict;
use warnings;

my ($i,$z)=(0,0);

while ($i < 2333) {
	++$z;
	if ( $z%2 == 0 or $z%3 == 0 ) {
		++$i;
		print "$i: $z\n";
	}
}

__END__
a(4n) mod 6 = 0 ; a(4n+1) mod 6 = 2 ; a(4n+2) mod 6 = 3 ; a(4n+3) mod 6 = 4。
2333=583*4+1，那么a(2333)=583*6+2=3500.

？。没看懂 a(4n+1) 是哪种代数表达，……
