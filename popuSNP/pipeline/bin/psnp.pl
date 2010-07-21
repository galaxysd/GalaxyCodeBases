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

our $opts='i:s:l:c:f:m:qvbd';
our($opt_i, $opt_s, $opt_l, $opt_m, $opt_c, $opt_f, $opt_v, $opt_b, $opt_q, $opt_d);

our $desc='1.filter fq, 2.stats';
our $help=<<EOH;
\t-i Project path
\t-s sub-group list(./subgroup.lst) in format: /^subGroupID\\tRegEx\$/
\t-l relative path of megred.lst to project path (2soap)
\t-m relative path of GLF.lst to project path (3GLF)
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
$opt_l='2soap' if ! $opt_l;
$opt_m='3GLF' if ! $opt_m;
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

my $spLname=$opt_i.'/'.$opt_l.'/megred.lst';
my $glfLname=$opt_i.'/'.$opt_m.'/GLF.lst';
die "[x]-i $glfLname not exists !\n" unless -f $glfLname;

print STDERR "From [$opt_i]/[$opt_l]/megred.lst,[$opt_m]/GLF.lst with [$opt_s][$opt_f][$opt_c] for 4 & 5\n";
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

my %spData;
open L,'<',$spLname or die "[x]Error opening $spLname: $!\n";
while (<L>) {
	chomp;
	my ($Sample,$Chr,$Len,$file)=split /\t/;
	$Sample =~ s/\s/_/g;
	$spData{$Sample}{$Chr}=[$Len,$file];
}
close L;

my ($GLFCount,%SampleChrGlf,%ChrCount,%SampleCount,@Samples)=(0);
open L,'<',$glfLname or die "[x]Error opening $glfLname: $!\n";
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
	++$SampleSub{$_}{$subID} for @Smp;	# was planing to check whether a Sample is used more than once, Cancled.
}
close L;
print STDERR scalar keys %SubSample,"]\n";
for my $v (values %SubSample) {
	$v = [sort keys %$v];	# This is where we sort glf lists. And %->% becomes %->@ here.
}
print STDERR " $_ -> [",scalar @{$SubSample{$_}},'] ',join(',',@{$SubSample{$_}}),"\n" for sort keys %SubSample;

### Let's Begin ###
my $outpath=$opt_i.'/4pSNP/';
#system('mkdir','-p',$outpath);
for my $Sub (keys %SubSample) {
	my $prefix=join('',$outpath,$Sub,'/');
	system('mkdir','-p',$prefix);
	open GLFA,'>',$prefix.'GLF.lst' || die "$!\n";
	open SP,'>',$prefix.'megred.lst' || die "$!\n";
	system('mkdir','-p',$prefix.'lst');
	system('mkdir','-p',$prefix.'sh');
	for my $Chr (@ChrIDs) {
		my $dir = $prefix.$Chr;
		my $glflist=join('',$prefix,'lst/',$Chr,'.glflst');
		open GLFL,'>',$glflist || die "$!\n";
		system('mkdir','-p',$dir);

		for my $Sample (@{$SubSample{$Sub}}) {
			print GLFA join("\t",$Sample,$Chr,$SampleChrGlf{$Sample}{$Chr}),"\n";
			print SP join("\t",$Sample,$Chr,@{$spData{$Sample}{$Chr}}),"\n";
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
	close SP;
	close GLFA;
}

### CN_LC_RST ###
$outpath=$opt_i.'/5FinalSNP/';
system('mkdir','-p',$outpath.'sh');
for my $Sub (keys %SubSample) {
	my $prefix=join('',$outpath,$Sub,'/1Population/');
	my $prefix2=join('',$outpath,$Sub,'/2Genotype/');
	my $prefix3=join('',$outpath,$Sub,'/3final/');
	system('mkdir','-p',$prefix);
	system('mkdir','-p',$prefix2);
	system('mkdir','-p',$prefix3);
	open CMD,'>',"${prefix}CnRstLc.cmd";
	open LST,'>',"${prefix}add_cn.lst";
	open CMDG,'>',"${prefix}GLF2GPF.cmd";
	open CMDF,'>',"${prefix}FILTER.cmd";
	open CMDH,'>',"${prefix}HW.cmd";
	open CMDA,'>',"${prefix}ADDCN.cmd";
	open CFG,'>',"${prefix}filterCRL.cfg.example";
	print CFG "# This is the CRL filtering configure for SubGroup:[$Sub]
# Please specify the following parameters for keeping data:
# depth\t[min,max]
# quality\t[min,+∞)
# rst\t[min,+∞)
# copynumber\t(-∞,max]
# REMEMBER TO ReName [filterCRL.cfg.example] to [filterCRL.cfg] BEFORE Running [${outpath}sh/step3_${Sub}_last.sh] !!!
depth\t60\t200
quality\t20
rst\t0.01
copynumber\t1.5
";
	close CFG;
	my $t=0;	# also scalar @ChrIDs
	for my $Chr (@ChrIDs) {
		print CMD "-i ${opt_i}/4pSNP/${Sub}/$Chr.psnp -r $opt_f/$Chr.fa -l $opt_c -c $Chr -m ${opt_i}/4pSNP/${Sub}/megred.lst -o ${prefix}$Chr\n";
		print LST "${prefix}${Chr}.add_cn\n";
		print CMDG "$Chr ${prefix}${Chr}.add_cn.dbsnp ${prefix}${Chr}.indsnp\n";
		print CMDF "${prefix}${Chr}.indsnp ${prefix2}${Chr}\n";
		print CMDH "-snp ${prefix}${Chr}.add_cn.filter -rm ${prefix2}${Chr}.rm -o ${prefix3}${Chr}_final.snp\n";
		print CMDA "${prefix2}${Chr}.add.filter $opt_f/$Chr.fa ${prefix3}${Chr}.add_ref\n";
		++$t;
	}
	close CMDA;
	close CMDH;
	close CMDF;
	close CMDG;
	close LST;
	close CMD;
	open SH,'>',"${outpath}sh/step1_${Sub}.sh";
	print SH "#!/bin/sh
#\$ -N \"Ps1_${Sub}_$$\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=6G
#\$ -hold_jid \"WSNP_*_${Sub}_$$\"
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-$t
SEEDFILE=${prefix}CnRstLc.cmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/copyNumLcRst.pl \$SEED
";
	close SH;
	open SH,'>',"${outpath}sh/step2_${Sub}_stat1.sh";
	print SH "#!/bin/sh
#\$ -N \"Ps2_${Sub}_stat_$$\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=100M
#\$ -hold_jid \"Ps1_${Sub}_$$\"
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash
perl $SCRIPTS/count.pl ${prefix}add_cn.lst 0addcn ${prefix}stat 2>${prefix}1stat.log
";
	close SH;

	open SH,'>',"${outpath}sh/step3_${Sub}_last.sh";
	print SH "#!/bin/sh
#\$ -N \"Ps3_${Sub}_last_$$\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=150M
#\$ -hold_jid \"Ps2_${Sub}_stat_$$\"
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-$t
SEEDFILE=${prefix}add_cn.lst
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
CMDFILE=${prefix}GLF2GPF.cmd
CMD=\$(sed -n -e \"\$SGE_TASK_ID p\" \$CMDFILE)
FILTERFILE=${prefix}FILTER.cmd
FILTER=\$(sed -n -e \"\$SGE_TASK_ID p\" \$FILTERFILE)
ADDCNFILE=${prefix}ADDCN.cmd
ADDCN=\$(sed -n -e \"\$SGE_TASK_ID p\" \$ADDCNFILE)
HWFILE=${prefix}HW.cmd
HW=\$(sed -n -e \"\$SGE_TASK_ID p\" \$HWFILE)
eval perl $SCRIPTS/filter_addcn.pl ${prefix}filterCRL.cfg \$SEED 2>${prefix}2filterCRL_\${SGE_TASK_ID}.log
eval python $SCRIPTS/GLF2GPF.py $opt_i/4pSNP/$Sub/GLF.lst \$CMD 2>${prefix}3Glf2Gpf_\${SGE_TASK_ID}.log
eval perl $SCRIPTS/gtypeHady_x2.pl ",scalar @{$SubSample{$Sub}}," \$FILTER
eval perl $SCRIPTS/H-W_filter.pl \$HW
eval perl $SCRIPTS/gtypeSnp_add_ref.pl \$ADDCN
";
	close SH;
}


#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
