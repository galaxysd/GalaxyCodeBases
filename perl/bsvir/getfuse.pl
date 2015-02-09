#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my $minLen = 5;

die "Usage: $0 <ins_size> <hum_sorted_bam> <vir_sorted_bam> [out_prefix]\n" if @ARGV < 3;
my ($isize,$inhum,$invir,$out)=@ARGV;

unless (defined $out) {
	$out='fuse_'.$inhum;
	$out =~ s/\.[sb]am(\.gz)?//g;
}
warn "[!]IN: h=[$inhum],v=[$invir],i=[$isize] to [$out].*\n";

sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.xz$/) {
		open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.bz2$/) {
		open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}

sub getsamChrLen($) {
	my $in = $_[0];
	my (%ChrLen,$Ref,$fq1,$fq2,$readlen1,$readlen2);
	open( IN,"-|","samtools view -H $in") or die "Error opening $in: $!\n";
	while (<IN>) {
		chomp;
		my ($t,$id,$ln)=split /\t/;
		if ($t eq '@SQ') {
			$id = (split /\:/,$id)[1];
			$ln = (split /\:/,$ln)[1];
			$ChrLen{$id} = $ln;
		} elsif ($t eq '@PG') {
# @PG	ID:BSMAP	VN:2.87	CL:"bsmap -u -d HomoGRCh38.fa -a F12HPCCCSZ0010_Upload/s00_C.bs_1.fq.gz -b F12HPCCCSZ0010_Upload/s00_C.bs_2.fq.gz -z 64 -p 12 -v 10 -q 2 -o s00_C.bs.bam"
			$Ref = $1 if /\-d\s+([.\w]+)\b/;
			$fq1 = $1 if /\-a\s+([^\s"]+)(\s|"?$)/;
			$fq2 = $1 if /\-b\s+([^\s"]+)(\s|"?$)/;
			my $FQ1 = openfile($fq1);
			<$FQ1>;chomp($_=<$FQ1>);
			$readlen1 = length $_;
			close $FQ1;
			my $FQ2 = openfile($fq2);
			<$FQ2>;chomp($_=<$FQ2>);
			$readlen2 = length $_;
			close $FQ2;
		}
	}
	return [[$Ref,$fq1,$fq2,$readlen1,$readlen2],%ChrLen];
}

my ($ret,%ChrLenHum,$finHum,%ChrLenHBV,$finHBV);

$ret = getsamChrLen($inhum);
$finHum = shift @$ret; %ChrLenHum = @$ret;

$ret = getsamChrLen($invir);
$finHBV = shift @$ret; %ChrLenHBV = @$ret;

#ddx \%ChrLenHum;
ddx \%ChrLenHBV;
ddx $finHum;
ddx $finHBV;

warn "[!]Ref: h=[$$finHum[0]],v=[$$finHBV[0]]\n";

open INHUM,"-|","samtools view $inhum" or die "Error opening $inhum: $!\n";
while (<INHUM>) {
	chomp;
	my ($id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/;
}


__END__
./getfuse.pl 200 grep_s00_C.bs.virsam.sort.bam grep_s00_C.bs.sort.bam
记得把流程中的文件名改下，突出ref是哪个物种。
