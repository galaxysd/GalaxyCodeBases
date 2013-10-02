#!/bin/env perl
use strict;
use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

die "Usage: $0 <input1> <input2> <output>\n" if @ARGV < 3;
my ($inf1,$inf2,$outf)=@ARGV;

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
print "Kmer Length = $klen1,$klen2\n";
my %TMP;

my ($flag,@RawArray,$la,$lb,$kmer1,$kmer2,$count1,$count2)=(3);
while ($flag) {
	unless ( defined($la=<$IN1>) ) {
		$flag &=2;
		$la = "@\t0";
	}
CYCLEb:
	unless ( defined($lb=<$IN2>) ) {
		$flag &=1;
		$lb = "@\t0";
	}
CYCLE:
	chomp($la,$lb);
	($kmer1,$count1) = split /\t/,$la;
	($kmer2,$count2) = split /\t/,$lb;
	if ( $kmer1 lt $kmer2 ) {
		if ( $kmer1 ne '@' ) {
			++$RawArray[$count1][0];	# $kmer2 is bigger thus cannot be '@'
++$TMP{"$count1."}{"$kmer1.\@x"};
			unless ( defined($la=<$IN1>) ) {
				$flag &=2;
				$la = "@\t0";
			}
			print STDERR "$kmer1 < $kmer2 $count1,$count2\n";
			goto CYCLE;
		}
	} elsif ( $kmer1 gt $kmer2 ) {
		++$RawArray[0][$count2];
++$TMP{".$count2"}{"x\@.$kmer2"};
		print STDERR "$kmer1 > $kmer2 $count1,$count2\n";
		goto CYCLEb;
	} elsif ( $count1 + $count2 ) {
		++$RawArray[$count1][$count2];
++$TMP{"$count1.$count2"}{"$kmer1"};
		print STDERR "$kmer1 = $kmer2 $count1,$count2\n";
	}
}
die if defined <$IN2>;
print STDERR "\n";

close $IN1;
close $IN2;

ddx \%TMP;

__END__
find sss2/ -name *.?z* > kmerfreq.lst

perl plotkmercnt.pl a.gz b.gz t >tt1 2>tt2
