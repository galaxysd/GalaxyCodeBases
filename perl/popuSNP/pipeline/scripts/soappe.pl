#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <insert size> <ReadLen> <fq1> <fq2> <insert size> <ref> <output prefix>\n";
	exit;
}
# $out.soap, $out.single, $out.log
my ($insize,$readlen,$fq1,$fq2,$ref,$out)=@ARGV;
my $bin='/panfs/GAG/huxuesong/scripts/soap2.20';
my $arg0='-p 6 -t -s 40 -l 32';
open LEN,'<',$readlen or die "[x]Error opening $readlen: $!\n";
$readlen = <LEN>;
chomp $readlen;
close LEN;
my $mismatch=$readlen>70?3:1;
open SIZE,'<',$insize or die "[x]Error opening $insize: $!\n";
$insize=<SIZE>;
chomp $insize;
my ($min,$max)=split /\t/,$insize;
close SIZE;

# always put -abDo2 in the first for the poor case without break ...
my $sh="$bin -a $fq1 -b $fq2 -D $ref -o $out.soap -2 $out.single $arg0 -m $min -x $max -v $mismatch 2>$out.log";
my ($Pairs,$Paired,$Singled);

TEST:
if (-s "$out.nfo") {
	system("mv -f ${out}_soappe.sh.archive ${out}_soappe.sh.archive.old") if (-e "${out}_soappe.sh.archive");
} else {
	open OUT,'>',"${out}_soappe.sh.archive" or warn "[!]Error opening ${out}_soappe.sh.archive: $!\n";
	print OUT "#!/bin/sh\n$sh\n";
	close OUT;
	system($sh)==0 or die "[x]system [$sh] failed: $?";
	open LOG,'<',"$out.log" or die "[x]Error opening $out.log: $!\n";
	while (<LOG>) {
		$Pairs = (split)[-2] if /^Total Pairs:/;
		$Paired = (split)[1] if /^Paired:/;
		$Singled = (split)[1] if /^Singled:/;
	}
	close LOG;
	unless ($Pairs) {
		system("mv -f $out.log $out.log.0");
		goto TEST;
	}
	open NFO,'>',"$out.nfo" or die "[x]Error opening $out.nfo: $!\n";
	print NFO "Total Pairs:\t$Pairs\nPaired:\t$Paired\nSingled:\t$Singled\n";
	close NFO;
}
