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

our $opts='i:o:l:v:bqd';
our($opt_i, $opt_o, $opt_v, $opt_b, $opt_q, $opt_d, $opt_l);

our $desc='1.filter fq, 2.stats';
our $help=<<EOH;
\t-i fqs.nfo path (./1fqfilted) with [fqs.lst] and [fqs.stat]
\t-l Limit NFO (./limit.lst) in format: /^SampleID\\tFC_Lane\\tMaxBP\$/
\t-o Output path (./1fqlimited), will mkdir if not exist
\t-q run qsub automatically
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

$opt_i='./1fqfilted' if ! $opt_i;
$opt_l='./limit.lst' if ! $opt_l;
$opt_o='./1fqlimited' if ! $opt_o;

system('mkdir','-p',$opt_o);
die "[x]-i $opt_i not exists !\n" unless -d $opt_i;
die "[x]-l $opt_l error !\n" unless -s $opt_l;

print STDERR "From [$opt_i] to [$opt_o] with [$opt_l]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN


my $infile=$opt_i.'/fqs.lst';
#my ($maxRL,%SampleRL,$readlen,$min,$max)=(0);
my ($opath,@cmdlines,%SampleDAT);

my %LimitDAT;
open LST,'<',$opt_l or die "[x]Error opening $opt_l: $!\n";
while (<LST>) {
	chomp;
	my ($Sample,$FL,$BP)=split /\t/;
	$LimitDAT{$Sample}{$FL}=$BP;
}
close LST;

open NFO,'>',$opt_o.'/fqs.lst';	# Well, auto. caltulation is a future function ...
open LST,'<',$infile or die "[x]Error opening $infile: $!\n";
while (<LST>) {
	chomp;
	my ($PESE,$sample,$lib,$FL,$ins,$ext,$path,@fqs)=split /\t/;
	#$SampleDAT{$sample}{$lib}{$FL}=[$PESE,$ins,$ext,$path,@fqs];
	$opath=$opt_o.'/'.$lib.'='.$sample.'/';
	system('mkdir','-p',$opath);
	print NFO join("\t",$PESE,$sample,$lib,$FL,$ins,$ext,$opath,@fqs),"\n";
	if (defined $LimitDAT{$sample}{$FL}) {
		push @cmdlines,join(' ',$LimitDAT{$sample}{$FL},$path,$opath,$ext,@fqs);
	} else {
		system("ln -s `readlink -nf $path${_}$ext` $opath") for @fqs;
		system("ln -s `readlink -nf $path${_}.nfo` $opath") for @fqs;
		#system("find $path -type f|while read a;do b=`readlink -nf \$a`;ln -s \$b $opath;done");
	}
}
close LST;
close NFO;

if ($opt_v > 3) {
	print "[$_]\n" for @cmdlines;
	print '-' x 80,"\n";
}

open SH,'>',$opt_o.'/cmd.lst' or die "[x]Error $!\n";
print SH "$_\n" for @cmdlines;
close SH;
open SH,'>',$opt_o.'/job.sh' or die "[x]Error $!\n";
print SH "#!/bin/sh
#\$ -N \"Pfilterlim\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=276M
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-",scalar @cmdlines,"
SEEDFILE=${opt_o}/cmd.lst
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/fqlimitor.pl \$SEED
";
close SH;

open SH,'>',$opt_o.'/stat.sh' or die "[x]Error $!\n";
print SH "#!/bin/sh
#\$ -N \"Pstatfqlim\"
#\$ -hold_jid \"Pfilterlim\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=30M
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash
perl $SCRIPTS/fqsummer.pl $opt_o/fqs.lst $opt_o/fqs.nfo $opt_o/fqs.stat
";
close SH;


#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
#print "\nPlease use the following command to batch qsub:\033[32;1m
#find ${opt_o}/ -name '*.sh' | while read ll; do qsub -l vf=2G -cwd \$ll; done\n\033[0;0m\n";
__END__
