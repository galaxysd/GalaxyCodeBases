#!/bin/env perl
use strict;
use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

my $GRID = 20;
my $GRIDp1 = $GRID - 1;
my $minGRIDstep = 2;

die "Usage: $0 <max_freq> <input1> <input2> <output>\n" if @ARGV < 4;
my ($max,$inf1,$inf2,$outf)=@ARGV;

die if $max < $GRID * $minGRIDstep;
my $oneGrid = int(1 + $max/$GRID/$minGRIDstep) * $minGRIDstep;

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

sub getKsize($) {
	my $fh = $_[0];
	my $tmp = <$fh>;
	unread $fh, $tmp;
	my $str = (split /\t/,$tmp)[0];
	my $len = length $str;
	#my $allcnt = 4 ** $len;
	return $len;
}

my $IN1 = openfile($inf1);
my $IN2 = openfile($inf2);
my ($klen1,$klen2) = ( getKsize($IN1),getKsize($IN2) );
if ($klen1 != $klen2) {
	die "$klen1,$klen2";
}
print "Kmer_Length: $klen1,$klen2\nGrid: $oneGrid x $GRID = ",$oneGrid * $GRID," of $max, ",$max/($GRID*$oneGrid),"\n";

my %TMP;

my ($flag,$lastunread,@RawArray,$la,$lb,$kmer1,$kmer2,$count1,$count2)=(3,0);
for my $i (0 .. $GRIDp1) {
	for my $j (0 .. $GRIDp1) {
		$RawArray[$i][$j] = 0;
	}
}

while ($flag) {
	unless ( defined($la=<$IN1>) ) {
		$flag &=2;
		$la = "@\t0";
	}
	last unless $flag;
	unless ( defined($lb=<$IN2>) ) {
		$flag &=1;
		$lb = "@\t0";
	}
	#chomp($la,$lb);
	($kmer1,$count1) = split /\t/,$la;
	($kmer2,$count2) = split /\t/,$lb;
	chomp($count1,$count2);
	$count1 = int($count1 / $oneGrid);
	$count2 = int($count2 / $oneGrid);
	if ( $count1 > $GRIDp1 or $count2 > $GRIDp1 ) {
		warn "[!] $count1,$count2 vs $GRIDp1\n";
	}
	if ( $kmer1 lt $kmer2 ) {
		if ( $kmer1 ne '@' ) {
			++$RawArray[$count1][0];	# $kmer2 is bigger thus cannot be '@'
++$TMP{"$count1."}{"$kmer1.\@x"};
			print STDERR "$kmer1 < $kmer2 $count1,$count2\n";
			if ($lastunread == 2) {
				last;
			}
			unread $IN2,$lb;
			$lastunread = 2;
		}
	} elsif ( $kmer1 gt $kmer2 ) {
		if ( $kmer2 ne '@' ) {
			++$RawArray[0][$count2];
++$TMP{".$count2"}{"x\@.$kmer2"};
			print STDERR "$flag $kmer1 > $kmer2 $count1,$count2\n";
			if ($lastunread == 1) {
				last;
			}
			unread $IN1,$la;
			$lastunread = 1;
		}
	} elsif ( $count1 + $count2 ) {
		++$RawArray[$count1][$count2];
++$TMP{"$count1.$count2"}{"$kmer1"};
		print STDERR "$flag $kmer1 = $kmer2 $count1,$count2\n";
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

ddx \%TMP;
ddx \@RawArray;

open OUT,'>',$outf or die;
print OUT "# Input: [$inf1],[$inf2]\n";
print OUT "# Kmer_Length: $klen1,$klen2\n# Grid: $oneGrid x $GRID = ",$oneGrid * $GRID," of $max, ",$max/($GRID*$oneGrid),"\n";
print OUT "# Horizontal:[$inf2], Vertical: [$inf1]\n";
print OUT join(',','.',( 0 .. $GRIDp1 )),"\n";
for my $i ( 0 .. $GRIDp1 ) {
	print OUT join(', ',$i,@{$RawArray[$i]}),"\n";
}
close OUT;



__END__
find sss2/ -name *.?z* > kmerfreq.lst

perl plotkmercnt.pl 9376518 a.gz b.gz t >tt1 2>tt2
