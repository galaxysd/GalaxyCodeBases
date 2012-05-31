#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ PKU <galaxy001@gmail.com>
Version: 1.0.0 @ 20120530
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <sam1> <sam2> <out prefix>\n" if @ARGV<2;
my $in1=shift;
my $in2=shift;
my $outp=shift;

sub openfile($) {
    my ($filename)=@_;
    my ($infile,$lastline);
    if ($filename=~/.bz2$/) {
	    open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	my $line;
	while (defined($line=<$infile>)) {
		if ($line =~ /^@\w\w\t\w\w:/) {
			next;
		}
		chomp $line;
		my @items = split /\t/,$line;
#warn "[",join(" | ",@items),"]\n";
		if ($items[1] & 16) {	# reverse
			next;
		}
		return [$infile,\@items];
	}
	die "File Error.";
}

my $samin1 = openfile($in1);
my $samin2 = openfile($in2);
open O,'|-',"gzip -9c >$outp.cmpsam.gz" or die "Error opening $outp.cmpsam.gz with gzip: $!\n";
open L,'>',$outp.'.log' or die "Error opening $outp.log: $!\n";
select(L);
$|=1;
select(STDOUT);
print L "From [$in1,$in2] to [$outp.cmpsam.gz]\n";
my ($Total,$Out,$notOut)=(0,0,0);

sub getNextRead($) {
	my $FHa=$_[0];
	my $t=$FHa->[1];
	$FHa->[1] = undef;
	my $line;
	while (defined($line = readline($$FHa[0]))) {
		chomp $line;
		my @items = split /\t/,$line;
		if ($items[1] & 128) {	# second read
			next;
		}
		$FHa->[1] = \@items;
		return $t;
	}
}

my ($StatSum,$arr1,$arr2,%Stat,$XT1,$XT2)=(0);
while (1) {
	$arr1 = getNextRead($samin1) or last;
	$arr2 = getNextRead($samin2) or last;
	die "$$arr1[0] ne $$arr2[0]" if $$arr1[0] ne $$arr2[0];
	my $diff=0;
	$diff |= 1 if (($$arr1[1] & 16) != ($$arr2[1] & 16));	# strand of the query
	$diff |= 2 if ($$arr1[2] ne $$arr2[2]);	# Reference sequence NAME
	$diff |= 4 if ($$arr1[3] ne $$arr2[3]);	# POS
	$diff |= 8 if (($$arr1[4]>0) != ($$arr2[4]>0));	# MAPQ == 0 or not.
	$diff |= 16 if ($$arr1[5] ne $$arr2[5]);	# CIAGR
	$XT1 = $XT2 = '*';
	if ($$arr1[2] eq '*') {
		$diff |= 256;
	} else {
		$XT1 = (grep(/^XT:/,@$arr1))[0];
		if (defined $XT1) {
			$XT1 = (split(':',$XT1))[2];
			$diff |= 64 if $XT1 eq 'U';
		} else {
			#warn "---[",join(" | ",@$arr1),"]\n";
			$diff |= 1024;
		}
	}
	if ($$arr2[2] eq '*') {
		$diff |= 512;
	} else {
		$XT2 = (grep(/^XT:/,@$arr2))[0];
		if (defined $XT2) {
			$XT2 = (split(':',$XT2))[2];
			$diff |= 128 if $XT2 eq 'U';
		} else {
			#warn "---[",join(" | ",@$arr2),"]\n";
			$diff |= 2048;
		}
	}
	$diff |= 32 if ( $XT1 and $XT2 and ($XT1 ne $XT2) );
	print O "[$#$arr1] [$#$arr2] $$arr1[0] - $diff\n";
	++$Stat{$diff};
	++$StatSum;
}
print L '# ',join(',',reverse qw/Strand Ref Pos MapQ>0 CIAGR XT 1=U 2=U 1unmap 2unmap 1noXT 2noXT/),"\n";
close $samin1->[0];
close $samin2->[0];
close O;

print L "Stat:\n";
for (sort {$Stat{$b}<=>$Stat{$a}} keys %Stat) {
	print L "$_\t$Stat{$_}\t";
	printf L "%.6f\t%#x\t%#014b\n",$Stat{$_}/$StatSum,($_) x 2;
}
print L "\n";
close L;
