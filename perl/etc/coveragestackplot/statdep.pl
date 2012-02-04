#!/bin/env perl
use strict;
use warnings;
use Inline C;

die "Usage: $0 <filled sam> <output>\n" if @ARGV <1;
my ($in,$out)=@ARGV;
$out=$in.'o' unless $out;
warn "From [$in] to [$out]\n";

sub openfile($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.bz2$/) {
	    open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
    return $infile;
}

__END__



__C__

    /* Find percentage of vowels to letters */
    double vowel_scan(char* str) {
        int letters = 0;
        int vowels = 0;
        int i = 0;
        char c;
        char normalize = 'a' ^ 'A';
        /* normalize forces lower case in ASCII; upper in EBCDIC */
        char A = normalize | 'a';
        char E = normalize | 'e';
        char I = normalize | 'i';
        char O = normalize | 'o';
        char U = normalize | 'u';
        char Z = normalize | 'z';

        while(c = str[i++]) {
            c |= normalize;
            if (c >= A && c <= Z) {
                 letters++;
                 if (c == A || c == E || c == I || c == O || c == U)
                     vowels++;
            }
        }

        return letters ? ((double) vowels / letters) : 0.0;
    }
