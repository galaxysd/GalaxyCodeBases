#!/bin/env perl
use strict;
use warnings;

die "usage: $0 <file prefix>\n" unless @ARGV == 1;

my (@cnt,@ref,@div);
open CNT,'<',$ARGV[0].'.sam.gc.stat' or die $!;
open REF,'<',$ARGV[0].'.ref.gc.stat' or die $!;

while(<CNT>) { last if /^#COL=higher/; }
while(<CNT>) {
	chomp;
	next if /^#/;
	#last if /^$/;
	my @t=split /\t/,$_;
	my $id=shift @t;
	$cnt[$id]=[@t];
}
close CNT;

while(<REF>) { last if /^#COL=higher/; }
while(<REF>) {
	chomp;
	next if /^#/;
	my @t=split /\t/,$_;
	my $id=shift @t;
	$ref[$id]=[@t];
}
close REF;

$cnt[$_][$_] *= 2  for (0..100);
$ref[$_][$_] *= 2  for (0..100);

sub dePepperC($$) {
	my ($h,$range)=@_;
	return 0 if $range > 2;
	my ($points,$sum)=(0,0);
	for my $i ($h-$range .. $h+$range) {
		next unless $div[$i][$i];
		next if $i==$h;
		++$points;
		$sum+=$div[$i][$i];
	}
	if ($points >= 2) {
		for my $i ($h-$range .. $h+$range) {
			print $div[$i][$i],' ';
		}
		print "\n",$sum / $points,"\n",'*' x 75,"\n";
		return $sum / $points;
	} else {
		return &dePepperC($h,++$range);
	}
}
sub dePepperA($$) {
	my ($h,$l)=@_;
	my $range=1;
	my ($points,$sum)=(0,0);
	for my $i ($h-$range .. $h+$range) {
		for my $j ($l-$range .. $l+$range) {
			next unless $div[$i][$j];
			next if $i==$h and $j==$l;
			++$points;
			$sum+=$div[$i][$j];
		}
	}
	if ($points >= 2) {
		for my $i ($h-$range .. $h+$range) {
			for my $j ($l-$range .. $l+$range) {
				print $div[$i][$j],' ';
			}
			print "\n";
		}
		print $sum / $points,"\n",'=' x 75,"\n";
		return $sum / $points;
	} else {
		return 0;
	}
}
sub dePepper($$$){
	my ($h,$l,$range)=@_;
	my ($points,$sum)=(0,0);
	for my $i ($h-$range .. $h+$range) {
		for my $j ($l-$range .. $l+$range) {
			next unless $ref[$i][$j];
			next if $i==$h and $j==$l;
			++$points;
			$sum+=$ref[$i][$j];
		}
	}
	if ($points >= 2) {
		for my $i ($h-$range .. $h+$range) {
			for my $j ($l-$range .. $l+$range) {
				print $ref[$i][$j],' ';
			}
			print "\n";
		}
		print $sum / $points,"\n";
		return $sum / $points;
	} else {
		return &dePepper($h,$l,++$range);
	}
}

my $refvalue=0;
for my $h (0..100) {
	for my $l (0..100) {
		#print "$h,$l\tCnt:$cnt[$h][$l]\tRef:$ref[$h][$l]\n";
		#next unless $cnt[$h][$l];
		$refvalue = $ref[$h][$l];
		unless ($refvalue) {
			if ($cnt[$h][$l]) {
				#$ref[$h][$l]=-1;
				$refvalue = &dePepper($h,$l,1);
				print '-' x 75,"\n";
			} else {$refvalue=1;}
		}
		$div[$h][$l] = $cnt[$h][$l] / $refvalue;
	}
}
my @divOut;
for my $h (0..100) {
	for my $l (0..100) {
		if ($h<$l) {
			$divOut[$h][$l] = 0;
		} elsif ($div[$h][$l]) {
			$divOut[$h][$l] = $div[$h][$l];
		} else {
			$divOut[$h][$l] = &dePepperA($h,$l);
			if ($h == $l and $divOut[$h][$l] == 0) {
				$divOut[$h][$l] = &dePepperC($h,1);
			}
		}
	}
}

open OUT,'>',$ARGV[0].'.div.gc.stat' or die $!;
print OUT "#COL=higher_GC%_read, ROW=lower_GC%_read\n";
print OUT join("\t",'#',0..100),"\n";
print OUT "$_\t",join("\t",@{$divOut[$_]}),"\n" for (0..100);
close OUT;
#print "$_: [",scalar @{$cnt[$_]},']',join(' ',@{$cnt[$_]}),"\n" for (0 .. $#cnt);
#print "$_: [",scalar @{$ref[$_]},']',join(' ',@{$ref[$_]}),"\n" for (0 .. $#ref);

