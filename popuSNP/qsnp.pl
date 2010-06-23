#!/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
use Galaxy::SGE::MakeJobSH 0.06;

$main::VERSION=0.0.3;
my $bin='/share/raid010/resequencing/user/huxs/release/popSNP/GLFmulti.py';

our $opts='i:o:c:f:s:bv';
our ($opt_i, $opt_c, $opt_o, $opt_f,  $opt_s, $opt_v, $opt_b);

our $desc='GLF2SNP Job file Maker';
our $help=<<EOH;
\t-i GLF path (./GLF) with SampleName_ChrID.glf, may not exists now
\t-c chrlen file (./chrlen) in format: /^ChrID\\tLen\\t?.*\$/
\t-f fragments length (1000000)
\t-s sample.list (sample.lst) in format: /^Sample\\t?.*\$/
\t-o popSNP output path (./popSNP), will mkdir if not exist
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='./GLF' if ! $opt_i;
$opt_c='./chrlen' if ! $opt_c;
$opt_f='1000000' if ! $opt_f;
$opt_f='1000000' if $opt_f<10;
$opt_s='sample.lst' if ! $opt_s;
$opt_o='./popSNP' if ! $opt_o;
$opt_i =~ s/\/$//;
$opt_o =~ s/\/$//;
print STDERR "From [$opt_i] with [$opt_c][$opt_s] frag. [$opt_f] to [$opt_o]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <STDIN>;}
--$opt_f;

my $start_time = [gettimeofday];
#BEGIN
#unless (-s $opt_c) {die "[x]chrlen file [$opt_c] is nothing !\n";}
my (%ChrLen,@Samples);
open CHRLEN,'<',$opt_c or die "[x]Error opening $opt_c: $!\n";
while (<CHRLEN>) {
	chomp;
	s/\r//g;
	my ($chrid,$len)=split /\s+/;	# I need a program to make chrlen ...
	next if ($chrid =~ /total/i);
	print "[$chrid]\t[$len]\n" if $opt_v;
	$ChrLen{$chrid}=$len;
}
close CHRLEN;
open SAMPLE,'<',$opt_s or die "[x]Error opening $opt_s: $!\n";
while (<SAMPLE>) {
	chomp;
	s/\r//g;
	my ($sample)=split /\t/;
	print "[$sample]\n" if $opt_v;
	push @Samples,$sample;
}
close SAMPLE;
system('mkdir','-p',"$opt_o/list");
system('mkdir','-p',"$opt_o/merge");
my ($pos1,$pos2,$part,$cmd,$chrlen,@parts,$partfile);
for my $chr (keys %ChrLen) {
	my $listfile="$opt_o/list/$chr.list";
	open LIST,'>',$listfile or die "[x]Error opening $listfile: $!\n";
	print LIST "$opt_i/$_/${_}_$chr.glf\n" for @Samples;
	close LIST;
	system('mkdir','-p',"$opt_o/$chr");
	$chrlen=$ChrLen{$chr};
	$part=$pos1=1;
	$pos2 = $pos1 + $opt_f;
	@parts=();
	while ($chrlen >= $pos2) {
		$partfile="$opt_o/$chr/${chr}_$part";
		push @parts,$partfile;
		$cmd="python $bin $pos1 $pos2 $listfile $partfile ERR";
		my $sh=Galaxy::SGE::MakeJobSH->new(vf=>'66M',cmd=>$cmd,name=>"${part}${chr}pSNP");
		$sh->markopt("-i $partfile");
		open SH,'>',"$opt_o/$chr/${chr}-$part.sh";
		print SH $sh->format;
		close SH;
		++$part;
		$pos1 = 1+$pos2;
		$pos2 = $pos1 + $opt_f;
	}
	if ($pos2 > $chrlen) {
		$partfile="$opt_o/$chr/${chr}_$part";
		push @parts,$partfile;
		$cmd="python $bin $pos1 $chrlen $listfile $partfile ERR";
		my $sh=Galaxy::SGE::MakeJobSH->new(vf=>'66M',cmd=>$cmd,name=>"${part}${chr}fpSNP");
		$sh->markopt("-i $partfile");
		open SH,'>',"$opt_o/$chr/${chr}-$part.sh";
		print SH $sh->format;
		close SH;
	}
	open SH,'>',"$opt_o/merge/${chr}-merge.sh";
	$cmd='cat '.join(' ',@parts)." > $opt_o/merge/${chr} LOG\nwc -l $opt_o/$chr/${chr}_* LOG\nwc -l $opt_o/merge/${chr} LOG";
	my $sh=Galaxy::SGE::MakeJobSH->new(vf=>'9M',cmd=>$cmd,name=>"${chr}mixpSNP");
	my $partsfile="$opt_o/merge/.${chr}-merge.list";
	open L,'>',$partsfile;
	print L "$_\n" for @parts;
	close L;
	$sh->waitopt("-li $partsfile");
	$sh->markopt("-i $opt_o/merge/${chr}");
	print SH $sh->format;
	close SH;
	print '.';	# dots every sample.
}
print "\nAll done !\n";



#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
print "\nPlease use the following command to batch qsub:\033[32;1m
find ${opt_o}/ -name '*.sh' |sort| while read ll; do sh \$ll; done\n\033[0;0m\n";
__END__
