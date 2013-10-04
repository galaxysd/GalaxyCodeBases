#!/bin/env perl
use strict;
use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

my $GRID = 40;	# 0 .. $GRID-1, but $GRID for bigger numbers.

die "Usage: $0 <max_freq> <input1> <input2> <output>\n" if @ARGV < 4;
my ($max,$inf1,$inf2,$outf)=@ARGV;

die if $max <= 0;
my $oneGrid = $max / $GRID;	# no need to +1 as max is for main parts of data only, 不包括极端值。

sub openfile($) {
	my ($filename)=@_;
	my ($KmerSum,$infile)=(0);
	if ($filename=~/.xz$/) {
			open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.lz4$/) {
		open( $infile,"-|","/opt/bin/lz4c -dy $filename /dev/stdout") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	open T,'<',"$filename.hist" or die $!;
	while (<T>) {
		my @t = split /\s+/;
print join('|',@t),"\n" if $t[0] eq '#';
		if ( $t[1] eq 'KmerSum:' ) {
			chomp($KmerSum = $t[2]);
		}
	}
	close T;
	return [$infile,$KmerSum];
}

sub getKsize($) {
	my $fh = $_[0];
	my $tmp = <$fh>;
	unread $fh, $tmp;
	my $str = (split /\t/,$tmp)[0];
	my $len = length $str;
	#my $allcnt = 4 ** $len;
	return $len;
}

my ($IN1,$KmerSum1) = @{openfile($inf1)};
my ($IN2,$KmerSum2) = @{openfile($inf2)};
my ($klen1,$klen2) = ( getKsize($IN1),getKsize($IN2) );
if ($klen1 != $klen2) {
	die "$klen1,$klen2";
}
print "Kmer_Length: $klen1,$klen2 KmerSum: $KmerSum1,$KmerSum2\nGrid: $oneGrid x $GRID = ",$oneGrid * $GRID," of $max, ",$max/($GRID*$oneGrid),"\n";

my %TMP;

my ($flag,$lastunread,@RawArray,$la,$lb,$kmer1,$kmer2,$count1,$count2)=(3,0);
for my $i (0 .. $GRID) {
	for my $j (0 .. $GRID) {
		$RawArray[$i][$j] = 0;
	}
}

while ($flag) {
	unless ( defined($la=<$IN1>) ) {
		$flag &=2;
		$la = "@\t-1";
	}
	last unless $flag;
	unless ( defined($lb=<$IN2>) ) {
		$flag &=1;
		$lb = "@\t-1";
	}
	#chomp($la,$lb);
	($kmer1,$count1) = split /\t/,$la;
	($kmer2,$count2) = split /\t/,$lb;
	chomp($count1,$count2);
#warn "[$flag] $kmer1,$count1\t$kmer2,$count2\n";

#	$count1 = int($count1 / $oneGrid);
#	$count2 = int($count2 / $oneGrid);
	$count1 = ($count1 / $KmerSum1) / $oneGrid;
	$count2 = ($count1 / $KmerSum2) / $oneGrid;
	$count1 = $GRID if $count1 > $GRID;
	$count2 = $GRID if $count2 > $GRID;

	if ( $kmer1 lt $kmer2 ) {
		#print STDERR "$kmer1 < $kmer2 $count1,$count2\n";
		if ( $kmer1 ne '@' ) {
			++$RawArray[$count1][0];	# $kmer2 is bigger thus cannot be '@'
#++$TMP{"$count1."}{"$kmer1.\@x"};
			#print STDERR "$kmer1 < $kmer2 $count1,$count2\n";
			unread $IN2,$lb;
			$lastunread = 2;
		} else {
			if ($lastunread == 2) {
				last;
			}
		}
	} elsif ( $kmer1 gt $kmer2 ) {
		#print STDERR "$flag $kmer1 > $kmer2 $count1,$count2\n";
		if ( $kmer2 ne '@' ) {
			++$RawArray[0][$count2];
#++$TMP{".$count2"}{"x\@.$kmer2"};
			#print STDERR "$flag $kmer1 > $kmer2 $count1,$count2\n";
			unread $IN1,$la;
			$lastunread = 1;
		} else {
			if ($lastunread == 1) {
				last;
			}
		}
	} elsif ( $count1 + $count2 >= 0 ) {
		++$RawArray[$count1][$count2];
#++$TMP{"$count1.$count2"}{"$kmer1"};
		#print STDERR "$flag $kmer1 = $kmer2 $count1,$count2\n";
	}
}
#my $t=<$IN2>;
#die "[$t][$la][$lb]" if defined $t;
##############################
# The last lines will be missing !!!
##############################
print STDERR "\n";

close $IN1;
close $IN2;

#ddx \%TMP;
#ddx \@RawArray;

open OUT,'>',$outf or die;
open OUTP,'>',$outf.'.plot' or die;
print OUT "# Input: [$inf1],[$inf2]\n";
print OUT "# Kmer_Length: $klen1,$klen2\n# Grid: $oneGrid x $GRID = ",$oneGrid * $GRID," of $max, ",$max/($GRID*$oneGrid),"\n";
print OUT "# Horizontal:[$inf2], Vertical: [$inf1]\n";
print OUT join(',','.',( 0 .. $GRID )),"\n";
for my $i ( 0 .. $GRID ) {
	print OUT join(', ',$i,@{$RawArray[$i]}),"\n";
	print OUTP join(', ',$i,@{$RawArray[$i]}),"\n";
}
close OUT;
close OUTP;


__END__
find sss2/ -name *.?z* > kmerfreq.lst

perl plotkmercnt.pl 40 a.gz b.gz t >tt1 2>tt2
perl plotkmercnt.pl 40 a100k.gz b100k.gz t100k.40
