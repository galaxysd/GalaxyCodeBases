#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use File::Basename;
use Galaxy::ShowHelp;
use FindBin qw($RealBin);

$main::VERSION=0.0.1;
my $SCRIPTS="$RealBin/../scripts";
my $POPSPLIT=1000000;

our $opts='i:s:l:c:f:qvbd';
our($opt_i, $opt_s, $opt_l, $opt_c, $opt_f, $opt_v, $opt_b, $opt_q, $opt_d);

our $desc='1.filter fq, 2.stats';
our $help=<<EOH;
\t-i Project path
\t-s sub-group list(./subgroup.lst) in format: /^subGroupID\\tRegEx\$/
\t-l relative path of GLF.lst to project path (3GLF)
\t-c Chromosome NFO file (chr.nfo) in format: /^ChrID\\s+ChrLen\\s?.*\$/
\t-f faByChr path (./faByChr) with ChrID.fa\(s\)
\t-q run qsub automatically
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

die "[x]Must specify Project path !\n" if ! $opt_i;
$opt_s='./subgroup.lst' if ! $opt_s;
$opt_l='3GLF' if ! $opt_l;
$opt_c='chr.nfo' if ! $opt_c;
$opt_f='./faByChr' if ! $opt_f;

$opt_i=~s#/+$##;
#$opt_i .= '/';
no warnings;
$opt_v=int $opt_v;
use warnings;

$opt_f=~s#/+$##;
# `readlink -f` will be blank if target not exists.
#system('mkdir','-p',$opt_o);

my $nfoname=$opt_i.'/'.$opt_l.'/GLF.lst';
die "[x]-i $nfoname not exists !\n" unless -f $nfoname;

print STDERR "From [$opt_i]/[$opt_l]/GLF.lst with [$opt_s][$opt_f][$opt_c]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my (@ChrIDs,%ChrLen);	# Canceled, use both. ;Well, use (keys %ChrCount) instead of @ChrIDs .
open CHRLEN,'<',$opt_c or die "[x]Error opening $opt_c: $!\n";
while (<CHRLEN>) {
	chomp;
	s/\r//g;
	my ($chr,$len)=split /\s+/;
	push @ChrIDs,$chr;
	$ChrLen{$chr}=$len;
}
close CHRLEN;
if ($opt_v) {
	print '[!]ChrID(s): [',join(',',@ChrIDs),"]\n";
}

my ($GLFCount,%SampleChrGlf,%ChrCount,%SampleCount,@Samples)=(0);
open L,'<',$nfoname or die "[x]Error opening $nfoname: $!\n";
print STDERR "[!]Glf (Sample x Chr = GLF): ";
while (<L>) {
	chomp;
	my ($Sample,$Chr,$Glf)=split /\t/;
	$Sample =~ s/\s/_/g;
	$SampleChrGlf{$Sample}{$Chr}=$Glf;
	++$SampleCount{$Sample};
	++$ChrCount{$Chr};
	++$GLFCount;
}
close L;
my $SampleCount=scalar keys %SampleCount;
my $ChrCount=scalar keys %ChrCount;
die "[x]Chr. Count Error.\n" if @ChrIDs != $ChrCount;
if ($SampleCount * $ChrCount == $GLFCount) {
	print STDERR $SampleCount,' x ',$ChrCount,' = ',$GLFCount,"\n";
} else {
	die "Error !\n[x]$SampleCount x $ChrCount = ",$SampleCount * $ChrCount," <> $GLFCount\n";
}
@Samples = keys %SampleChrGlf;	# no need to sort here

my (%SampleSub,%SubSample);
open L,'<',$opt_s or die "[x]Error opening $opt_s: $!\n";
print STDERR "[!]Sub Group(s): [";
while (<L>) {
	next if /^[#;]/;
	chomp;
	my ($subID,$RegEx)=split /\t/;
	my @Smp=grep /$RegEx/,@Samples;
	++$SubSample{$subID}{$_} for @Smp;
	++$SampleSub{$_}{$subID} for @Smp;
}
close L;
print STDERR scalar keys %SubSample,"]\n";
for my $v (values %SubSample) {
	$v = [sort keys %$v];
}
print STDERR " $_ -> [",scalar @{$SubSample{$_}},'] ',join(',',@{$SubSample{$_}}),"\n" for sort keys %SubSample;

### Let's Begin ###
my $outpath=$opt_i.'/4pSNP/';
system('mkdir','-p',$outpath);
for my $Sub (keys %SubSample) {
	my $prefix=join('',$outpath,$Sub,'/');
	open GLFA,'>',$prefix.'GLF.lst' || die "$!\n";
	system('mkdir','-p',$prefix.'lst');
	system('mkdir','-p',$prefix.'sh');
	for my $Chr (@ChrIDs) {
		my $dir = $prefix.$Chr;
		my $glflist=join('',$prefix,'lst/',$Chr,'.glflst');
		open GLFL,'>',$glflist || die "$!\n";
		system('mkdir','-p',$dir);

		for my $Sample (@{$SubSample{$Sub}}) {
			print GLFA join("\t",$Sample,$Chr,$SampleChrGlf{$Sample}{$Chr}),"\n";
			print GLFL $SampleChrGlf{$Sample}{$Chr},"\n";
		}
		close GLFL;

		open LST,'>',$prefix.'lst/'.$Chr.'.psnplst' || die "$!\n";
		open CMD,'>',$dir."/$Chr.popcmd" || die "$!\n";
		my $len=$ChrLen{$Chr};
		my $lstcount=0;
		my ($i,$j,$exit)=(1,$POPSPLIT,0);
		while ($exit != -1) {
			++$exit;
			if ($j > $len) {
				$j=$len;
				$exit=-1;
			}
			print LST "$dir/${Chr}_$exit.psnp\n";
			print CMD "$i $j $glflist $dir/${Chr}_$exit.psnp >$dir/${Chr}_$exit.tag 2>$dir/${Chr}_$exit.log\n";
			++$lstcount;
			$i += $POPSPLIT;
			$j += $POPSPLIT;
		}
		close CMD;
		close LST;
		open SH,'>',$prefix.'sh/'.$Chr.'_popsnp.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"pSNP_${Chr}_${Sub}_$$\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=120M,p=1
#\$ -hold_jid glf_$Chr
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-$lstcount
SEEDFILE=${dir}/$Chr.popcmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval python $SCRIPTS/GLFmulti.py \$SEED
";
		close SH;
		open SH,'>',$prefix.'sh/'.$Chr.'_wcsnp.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"WSNP_${Chr}_${Sub}_$$\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=60M
#\$ -hold_jid pSNP_${Chr}_${Sub}_$$
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash
cat ${prefix}lst/${Chr}.psnplst|xargs wc -l > $dir.psnpwc
rm -f $dir.psnp
cat ${prefix}lst/${Chr}.psnplst|xargs cat >> $dir.psnp
wc -l $dir.psnp >> $dir.psnpwc
";	# whenever >> , remember to rm first !
		close SH;
	}
	close GLFA;
}

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
