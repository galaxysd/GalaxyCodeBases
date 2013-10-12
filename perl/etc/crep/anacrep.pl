#!/bin/env perl
use strict;
use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

die "Usage: $0 <input> <cmplist(blank to quary)>\n" if @ARGV < 1;
my ($inf,$cmplist)=@ARGV;

my $skip = 10;

sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.xz$/) {
			open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.lz4$/) {
		open( $infile,"-|","/opt/bin/lz4c -dy $filename /dev/stdout") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}

my $IN = openfile($inf);
chomp(my $tmp = <$IN>);
my @arr = split /\t/,$tmp;

my @arr0 = @arr;

print "Read List:\n",'-' x 75,"\n";
splice(@arr,0,$skip);
my (%Experiments,%Samples,%Dat);
$tmp = $skip - 1;
for (@arr) {
	++$tmp;
	/^([^(]+)\(([^)]+)\)$/ or die;
	++$Samples{$1};
	++$Experiments{$2};
	push @{$Dat{$2}}, $tmp;
	print join("|",$tmp+1,$1,$2,$arr0[$tmp]),"\n";
	die if $arr0[$tmp] ne $_;
}
my @SamplesL = sort keys %Samples;
my @ExperimentsL = sort keys %Experiments;

my (@toCmp,@CmpPlans,@CmpPairs);
unless ($cmplist) {
	print "\nPlease Select Experiments: (like: 0,3,4)\n",'-' x 75,"\n";
	for $tmp ( 0 .. $#ExperimentsL ) {
		print join("\t",$tmp,$ExperimentsL[$tmp]),"\n";
	}
	exit(2);
} else {
	@toCmp = split /,/,$cmplist;
	for (@toCmp) {
		die if $_ != int($_);
		die if $_ < 0 or $_ > $#ExperimentsL;
	}
	for my $i (0 .. $#toCmp) {
		for my $j (($i+1) .. $#toCmp) {
			push @CmpPlans,[$toCmp[$i],$toCmp[$j]];
			push @CmpPairs,[ $Dat{$ExperimentsL[$toCmp[$i]]},$Dat{$ExperimentsL[$toCmp[$j]]} ];
		}
	}
#ddx \@CmpPairs;
	print "\nCmp plans:\n",'-' x 75,"\n";
	for my $i (0 .. $#CmpPlans) {
		print $CmpPlans[$i][0],'(',join('|',@{$CmpPairs[$i][0]}),') <-> ',
			$CmpPlans[$i]->[1],'(',join('|',@{$CmpPairs[$i][1]}),"): $ExperimentsL[$CmpPlans[$i][0]] --- $ExperimentsL[$CmpPlans[$i][1]]","\n";
		for ( @{$CmpPairs[$i][0]},@{$CmpPairs[$i][1]} ) {
			print " $_: ",$arr0[$_],"\n";
		}
	}
	#print '-' x 75,"\n";
}





close $IN;

