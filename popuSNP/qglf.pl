#!/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
use Galaxy::SGE::MakeJobSH;

$main::VERSION=0.1.1;
my $bin = "/share/raid010/resequencing/user/zhanghao/raid11/software/SOAPsnp-v1.02/soapsnp";

our $opts='i:o:s:m:c:a:bv';
our($opt_i, $opt_s, $opt_m, $opt_c, $opt_a, $opt_o, $opt_v, $opt_b);

our $desc='Merge2GLF Job file Maker';
our $help=<<EOH;
\t-i Merged soap path (./merged) with SampleName.ChrID that can be 'find'
\t-s sample.list (sample.lst)
 in format: /^Sample\\tLib\\tRead_Length\\tMonoPloid_ChrID1,MonoPloid_ChrIDs\\t?.*\$/
\t-c fabychr path with ChrID.fa
\t-a quality calibration matrix
\t-m ChrID to use monoploid calling mode like [ChrX,ChrY,ChrM] or [FILE]
\t-o GLF output path (./GLF), will mkdir if not exist
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='./merged' if ! $opt_i;
$opt_s='sample.lst' if ! $opt_s;
$opt_o='./GLF' if ! $opt_o;
die "[x]Must specify -c 'fabychr path' !\n" unless $opt_c;
die "[x]Must specify -a 'SoapSNP.matrix' !\n" unless $opt_a;

$opt_c =~ s/\/$//;
print STDERR "From [$opt_i] of [$opt_s] with [$opt_c] to [$opt_o], ";
my ($MP_File,@MonoPloids)=(0);
if ($opt_m and $opt_m eq 'FILE') {
		$MP_File=1;
		print STDERR "MonoPloid by file:[$opt_s]\n";
	} elsif ($opt_m) {
		@MonoPloids=split /,/,$opt_m;
		$opt_m=join '][',@MonoPloids;
		print STDERR 'MonoPloid for ',1+$#MonoPloids,":[$opt_m]\n";
	} else {print STDERR "MonoPloid OFF\n"}
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
unless (-s $opt_a) {die "[x]SoapSNP matrix [$opt_a] is nothing !\n";}
my (%ReadLen,%Chrs,@chrs);
open SAMPLE,'<',$opt_s or die "[x]Error opening $opt_s: $!\n";
while (<SAMPLE>) {
	chomp;
	s/\r//g;
	my ($sample,$len,$mp)=(split /\t/)[0,2,3];
	print "[$sample]\t[$len]\n" if $opt_v;
	$ReadLen{$sample}=$len;
	if ($MP_File == 1) {
		if ($mp) {@MonoPloids=split /,/,$mp;}
		 else {@MonoPloids=();}
	}
#print "[@MonoPloids]\n";
}
close SAMPLE;
$|=1;
my @RefFiles=`find $opt_c/`;
chomp @RefFiles;
my ($flag,$cmd,$ChrID,$ChrRef,$filename);
for my $sample (keys %ReadLen) {
	@chrs=`find $opt_i/ -name $sample.*`;
	chomp @chrs;
	system('mkdir','-p',"$opt_o/$sample");
#print "$#chrs [@chrs]\n";
	for my $file (@chrs) {
		$ChrID=(split /\./,$file)[-1];
		#$ChrRef=`find $opt_c/ -name $ChrID.*`;	# better to move outside of the cycle.
		($ChrRef)=grep /$ChrID\.\w+$/,@RefFiles;
#die $ChrRef;
		next unless $ChrRef;
		#$ChrRef = (split /\n/,$ChrRef)[0];
		#chomp $ChrRef;
		$filename=(split /\//,$file)[-1];
		$flag='';	# '' or '-m ', maybe '-m -r 0.0007 '
		for (@MonoPloids) {$flag='-m ' if $file =~ /$_$/;}
		$cmd=$bin." -F 1 -i $file -d $ChrRef -o $opt_o/$sample/${sample}_$ChrID.glf -I $opt_a ".$flag.'-L '.$ReadLen{$sample}.' LOG';
#print "[$cmd]\n";
		my $sh=Galaxy::SGE::MakeJobSH->new(vf=>'270M',cmd=>$cmd,name=>"$sample${ChrID}GLF");
		open SH,'>',"$opt_o/$sample/$ChrID.sh";
		print SH $sh->format;
		close SH;
	}
	print '.';	# dots every sample.
}
print "\nAll done !\n";

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
print "\nPlease use the following command to batch qsub:\033[32;1m
find ${opt_o}/ -name '*.sh' | while read ll; do sh \$ll; done\n\033[0;0m\n";
__END__
