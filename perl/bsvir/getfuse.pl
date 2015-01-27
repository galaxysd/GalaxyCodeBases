#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <ins_size> <hum_sorted_bam> <vir_sorted_bam> [out_prefix]\n" if @ARGV < 3;
my ($isize,$inhum,$invir,$out)=@ARGV;

unless (defined $out) {
	$out='fuse_'.$inhum;
	$out =~ s/\.[sb]am(\.gz)?//g;
}
warn "[!]IN: h=[$inhum],v=[$invir],i=[$isize] to [$out].*\n";

sub getsamChrLen($) {
	my $in = $_[0];
	my (%ChrLen,$Ref);
	open( IN,"-|","samtools view -H $in") or die "Error opening $in: $!\n";
	while (<IN>) {
		chomp;
		my ($t,$id,$ln)=split /\t/;
		if ($t eq '@SQ') {
			$id = (split /\:/,$id)[1];
			$ln = (split /\:/,$ln)[1];
			$ChrLen{$id} = $ln;
		} elsif ($t eq '@PG') {
			$Ref = $1 if /\-d\s+([.\w]+)\b/;
		}
	}
	return [$Ref,%ChrLen];
}

my ($ret,%ChrLenHum,$fRefHum,%ChrLenHBV,$fRefHBV);

$ret = getsamChrLen($inhum);
$fRefHum = shift @$ret; %ChrLenHum = @$ret;

$ret = getsamChrLen($invir);
$fRefHBV = shift @$ret; %ChrLenHBV = @$ret;

#ddx \%ChrLenHum;
ddx \%ChrLenHBV;

warn "[!]Ref: h=[$fRefHum],v=[$fRefHBV]\n";

open INHUM,"-|","samtools view $inhum" or die "Error opening $inhum: $!\n";
while (<INHUM>) {
	chomp;
}


__END__
./getfuse.pl 200 grep_s00_C.bs.virsam.sort.bam grep_s00_C.bs.sort.bam
记得把流程中的文件名改下，突出ref是哪个物种。
