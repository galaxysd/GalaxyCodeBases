#!/bin/env perl
use strict;
use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

die "Usage: $0 <input> <cmplist(blank to quary)>\n" if @ARGV < 1;
my ($inf,$cmplist)=@ARGV;

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

print "Read List:\n",'-' x 75,"\n";

my (%Experiments,%Samples,%Dat);
for my $i ( 0 .. $#arr ) {
	print join("|",$i,$arr[$i]),"\t";
	$arr[$i] =~ /^([^(]+)\(([^)]+)\)$/ or (print "<=\n" and next);
	++$Samples{$1};
	++$Experiments{$2};
	push @{$Dat{$2}}, $i;
	print join("|",$1,$2),"\n";
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
			print " $_: ",$arr[$_],"\n";
		}
	}
	#print '-' x 75,"\n";
}

ddx \@CmpPlans;
ddx \@CmpPairs;

open O,'>',$inf.'.up2' or die;

while (<$IN>) {
	chomp;
	@arr = split /\t/;
	my (@cmpA,@cmpB);
	for my $p ( @CmpPairs ) {
		@cmpA = @{$p->[0]};
		@cmpB = @{$p->[1]};
		my ($i,$avgA,$avgB)=(0);
		for (@cmpA) {
			$avgA += $arr[$_];
			++$i;
		}
		$avgA /= $i;
		$i=0;
		for (@cmpB) {
			$avgB += $arr[$_];
			++$i;
		}
		$avgB /= $i;
		my $ratio = -1;
		$ratio = $avgA / $avgB if ($avgB);
		if ($ratio >= 2) {
			print O join("\t",$arr[0],$ratio,int(0.5+100*$avgA)/100,int(0.5+100*$avgB)/100,join('|',@arr[@cmpA]),join('|',@arr[@cmpB])),"\n";
			for ( @cmpA,@cmpB ) {
				print " $_: ",$arr[$_],"\t";
			}
			print "\n",'-'x5," $avgA / $avgB = $ratio\n";
		}
	}

}

close $IN;
close O;

__END__
perl anacrep.pl crep_all_tsv_new.txt 12,13 > crep_all_tsv_new.txt.out
