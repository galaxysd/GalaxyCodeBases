#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <ReadLen> <fq> <ref> <output prefix>\n";
	exit;
}
# $out.se, $out.log
my ($readlen,$fq,$ref,$out)=@ARGV;
my $bin='/panfs/GAG/huxuesong/scripts/soap2.20';
my $arg0='-p 6 -t -s 40 -l 32';
open LEN,'<',$readlen or die "[x]Error opening $readlen: $!\n";
$readlen = <LEN>;
chomp $readlen;
close LEN;
my $mismatch=$readlen>70?3:1;

my $sh="$bin -a $fq -D $ref -o $out.se $arg0 -v $mismatch 2>$out.log";
my ($Reads,$Alignment);

TEST:
if (-s "$out.nfo") {
	system("mv -f ${out}_soapse.sh.archive ${out}_soapse.sh.archive.old") if (-e "${out}_soapse.sh.archive");
} else {
	open OUT,'>',"${out}_soapse.sh.archive" or warn "[!]Error opening ${out}_soapse.sh.archive: $!\n";
	print OUT "#!/bin/sh\n$sh\n";
	close OUT;
	system($sh)==0 or die "[x]system [$sh] failed: $?";
	open LOG,'<',"$out.log" or die "[x]Error opening $out.log: $!\n";
	while (<LOG>) {
		$Reads = (split)[-1] if /^Total Reads/;
		$Alignment = (split)[1] if /^Alignment:/;
	}
	close LOG;
	unless ($Reads) {
		system("mv -f $out.log $out.log.0");
		goto TEST;
	}
	open NFO,'>',"$out.nfo" or die "[x]Error opening $out.nfo: $!\n";
	print NFO "Total Reads:\t$Reads\nAlignment:\t$Alignment\n";
	close NFO;
}
