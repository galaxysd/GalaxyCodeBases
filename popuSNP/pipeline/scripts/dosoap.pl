#!/bin/env perl
use strict;
use warnings;

my $bin='/panfs/GAG/huxuesong/scripts/soap2.20';
my $arg0='-p 6 -t -s 40 -l 32';

unless (@ARGV){
	print "perl $0 <fq.lst> <fq.nfo> <fq.stat>\n";
	exit;
}

my ($PESE,$ins,$G,$opath,$Ref,$fqextpath,$fqname) = @ARGV;
#print '[',join('] [',@ARGV),"]\n";
#print '[',join('] [',$PESE,$ins,$G,$Ref,$fqextpath,@fqnames),"]\n";
my ($readlen,$min,$max)=split ',',$ins;
my ($ext,$path)=split ',',$fqextpath;
my @fqnames=split ',',$fqname;
my $mismatch=$readlen>70?3:1;
$arg0 .= ' -R' if $min > 1500;
my $fqcmd;
unless ($max==0) {	# PE
	$fqcmd="-a $path$fqnames[0]$ext -b $path$fqnames[-1]$ext -o $opath$fqnames[0].soap -2 $opath$fqnames[0].single";
} else {	# SE
	$fqcmd="-a $path$fqnames[0]$ext  -o $opath$fqnames[0].se";
}

# always put -abDo2 in the first for the poor case without break ...
my $sh="$bin $fqcmd -D $Ref $arg0 -m $min -x $max -v $mismatch -g $G 2>$opath$fqnames[0].log";
print "[$sh]\n";

my ($Pairs,$Paired,$Singled,$Reads,$Alignment);
TEST:
if (-s "$opath$fqnames[0].nfo") {
	system("mv -f $opath$fqnames[0]_soap.sharchive $opath$fqnames[0]_soap.sharchive.old") if (-e "$opath$fqnames[0]_soappe.sh.archive");
} else {
	open OUT,'>',"$opath$fqnames[0]_soap.sharchive" or warn "[!]Error opening $opath$fqnames[0]_soap.sharchive: $!\n";
	print OUT "#!/bin/sh\n$sh\n";
	close OUT;
	system($sh)==0 or die "[x]system [$sh] failed: $?";
	open LOG,'<',"$opath$fqnames[0].log" or die "[x]Error opening $opath$fqnames[0].log: $!\n";
	while (<LOG>) {
		$Pairs = (split)[-2] if /^Total Pairs:/;
		$Paired = (split)[1] if /^Paired:/;
		$Singled = (split)[1] if /^Singled:/;
		$Reads = (split)[-1] if /^Total Reads/;
		$Alignment = (split)[1] if /^Alignment:/;
	}
	close LOG;
	unless ($Pairs or $Reads) {
		system("mv -f $opath$fqnames[0].log $opath$fqnames[0].log.0");
		goto TEST;
	}
	open NFO,'>',"$opath$fqnames[0].nfo" or die "[x]Error opening $opath$fqnames[0].nfo: $!\n";
	if ($max==0) {
		print NFO "Total Reads:\t$Reads\nAlignment:\t$Alignment\n";
	} else {
		print NFO "Total Pairs:\t$Pairs\nPaired:\t$Paired\nSingled:\t$Singled\n";
	}
	close NFO;
}
