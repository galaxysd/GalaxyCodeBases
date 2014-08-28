#!/bin/env perl
use strict;
use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

die "Usage: $0 <input sam.gz/bam> <output>\n" if @ARGV < 2;
my ($inf,$outf)=@ARGV;

open CHR,'<','chr.lst' or die;
my %ChrGI2ID;
while (<CHR>) {
	my ($gi,$id) = split /\s+/,$_;
	$ChrGI2ID{$gi}=$id;
}
close CHR;
$ChrGI2ID{'='}='=';
ddx \%ChrGI2ID;

my $FH;
if ($inf =~ /\.sam\.gz$/i) {
	open $FH,'-|',"gzip -dc $inf" or die "Error opening $inf: $!\n";
} elsif ($inf =~ /\.bam$/i) {
	open $FH,'-|',"samtools view -h $inf" or die "Error opening $inf: $!\n";
} elsif ($inf =~ /\.sam$/i) {
	open $FH,'<',$inf or die "Error opening $inf: $!\n";
} else { die; }

my $line;
my %ChrLen;
while ($line = <$FH>) {
	unless ($line =~ /^\@/) {
		unread $FH, $line;
		last;
	}
	if ($line =~ /^\@SQ/) {
		my (undef,$gi,$len) = split /\s+/,$line;
		$gi =~ s/^SN://;
		$len =~ s/^LN://;
		$ChrLen{$gi} = $len;
	}
}
ddx \%ChrLen;

my %DatbyChrM;
#$DatbyChrM{$_}=[] for values %ChrGI2ID;
while ($line = <$FH>) {
	my ($id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/,$line;
	#print "$id, $flag, Chr$ChrGI2ID{$ref}, $pos, $mapq, $CIAGR, Chr$ChrGI2ID{$mref}, $mpos, $isize\n";
	my $commonChrID = $ChrGI2ID{$ref};
	my $posto1m = int($pos/1000000);
	++$DatbyChrM{$commonChrID}->[$posto1m];
}
#ddx \%DatbyChrM;

close $FH;

open OUT,'>',$outf or die "Error opening $outf: $!\n";
print OUT "# $inf\n";
for my $ccid (sort { "$a$b"=~/^\d+$/ ? $a<=>$b : $a cmp $b } keys %DatbyChrM) {
	print OUT "[$ccid]\n";
	my $ArrayRef = $DatbyChrM{$ccid};
	for my $i ( 0 .. $#$ArrayRef ) {
		my $t = $$ArrayRef[$i];
		print OUT (defined $t)?$t:0,"\t",1+$i,"\n";
	}
}

close OUT;

__END__
perl startplot.pl t.sam t.out
grep \\[ t.out
