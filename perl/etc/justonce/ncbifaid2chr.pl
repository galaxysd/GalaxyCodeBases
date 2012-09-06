#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use File::Basename;
use Galaxy::IO;

die "Usage: $0 <in> <outdir>\n" if @ARGV<2;
my $inf=shift;
my $outf=shift;

my($filename, $directories, $suffix) = fileparse($inf);
my $chr = (split /\./,$filename)[0];
print "From [$inf] to [$outf/ChrID.fa] for [$chr]\n";

my $fh = openfile($inf);
open O,'>',"$outf/$chr.fa" or die $!;
chomp(my $t=<$fh>);
die unless $t =~ s/^>//;
print O ">$chr $t\n";
while (<$fh>) {
	print O $_;
}
close $fh;
close O;

__END__
mkdir cat62
find Felis_catus-6.2/chr*.gz|while read a;do ./ncbifaid2chr.pl "$a" cat62 ;done &
