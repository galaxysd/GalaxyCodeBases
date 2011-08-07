#!/bin/env perl
use strict;
use warnings;

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

my (%Count,$CountRF,$Reads);
($CountRF,$Reads,$Count{0},$Count{1},$Count{-1})=(0,0,0,0,0);
my $fh=openfile($in);
while (<$fh>) {
	next if /^@\w\w\t\w\w:/;
	chomp;
	my @read=split /\t/, $_;
	my $Alternativehits='';
	for ( @read ) {
    	if (/^XA:Z:([\w,+-;]+)$/) { #XA:Z:chrX,+1144912,100M,0;
			$Alternativehits = $1;
			last;
		}
	}
	my @Alt=split ';',$Alternativehits;
	++$Reads;
	my $flag=0;
	for (@Alt) {
		my ($chr,$pos,$CIGAR,$NM)=split /,/;
		if (substr($pos,-2) eq '01') {
			if ($pos>0) {
				$flag=1;
			} else {$flag=-1;}
		}
	}
	++$Count{$flag};
	++$CountRF if $flag;
}
close $fh;
open O,'>',"$out" or die "[x]Error opening $out: $!\n";
print O $Count{1}+$Count{-1}," ,RF:$CountRF\t+$Count{1},-$Count{-1},z$Count{0}\t$Reads\t$in\n";
close O;
print $Count{1}+$Count{-1}," ,RF:$CountRF\t+$Count{1},-$Count{-1},z$Count{0}\t$Reads\n";
__END__
find . -name '*.r'|xargs -n1 ./samrstat.pl
find . -name '*.ro'|xargs -n1 cat|sort -n
