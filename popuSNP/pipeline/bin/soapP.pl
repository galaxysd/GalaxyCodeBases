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

our $opts='i:o:n:g:v:r:bqd';
our($opt_i, $opt_n, $opt_g, $opt_r, $opt_o, $opt_v, $opt_b, $opt_q, $opt_d);

our $desc='1.filter fq, 2.stats';
our $help=<<EOH;
\t-i fqs.nfo path (./1fqfilted) with fqs.nfo
\t-n fqs.nfo name (fqs.nfo)
\t-r Reference Genome for Soap2 (./Ref) with *.index.bwt
\t-o Output path (./2soap), will mkdir if not exist
\t-g the same as soap2 -g (undef=0), max is 10 for soap2.20
\t-q run qsub automatically
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

$opt_i='./1fqfilted' if ! $opt_i;
$opt_n='fqs.nfo' if ! $opt_n;
$opt_r='./Ref' if ! $opt_r;
$opt_o='./2soap' if ! $opt_o;

$opt_r=~s#/+$##;
no warnings;
$opt_v=int $opt_v;
$opt_g=int $opt_g;

my @t=`find $opt_r/ -name '*.index.???' 2>/dev/null`;
$t[0] =~ /(.+\.index)\.\w+$/;
my $ref = $1;
die "[x]*.index.bwt NOT found in [$opt_r] !\n" unless $ref;
use warnings;

$opt_i=~s#/+$##;
$opt_o=~s#/+$##;
# `readlink -f` will be blank if target not exists.
system('mkdir','-p',$opt_o);

my $nfoname=$opt_i.'/'.$opt_n;
die "[x]-i $nfoname not exists !\n" unless -f $nfoname;

print STDERR "From [$opt_i]/[$opt_n] to [$opt_o] refer to [$opt_r] as [$ref]\n";
print STDERR "Soap2 -g [$opt_g]\n" if $opt_g;
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my (%DATrbrf,%Librmx,%FQnfo,$readlen,$min,$max);
my ($opath,%cmdlines);
open LST,'<',$nfoname or die "[x]Error opening $nfoname: $!\n";
while (<LST>) {
	chomp;
	my ($PESE,$sample,$lib,$FL,$ins,$ext,$path,@fqs)=split /\t/;
	($readlen,$min,$max)=split ',',$ins;
	next if $readlen == -1;
	#$Librmx{$sample}{$lib}{$FL}=[$readlen,$min,$max];	# readlen,minIS,maxIS
	#$FQnfo{$sample}{$lib}{$FL}=[$ext,$path,@fqs];
	$opath="$opt_o/$sample/$lib/";
	system('mkdir','-p',$opath);
	push @{$cmdlines{$lib}},"$PESE $ins $opt_g '$opath' '$ref' '$ext,$path' '".join(',',@fqs)."'";
}
if ($opt_v > 3) {
	print '-' x 80,"\n";
	for my $lib (keys %cmdlines) {
		print "[$lib]\n";
		print "[$_]\n" for @{$cmdlines{$lib}};
		print "\n";
	}
	print '-' x 80,"\n";
}
$opath="$opt_o/sh/";
system('mkdir','-p',$opath);
for my $lib (keys %cmdlines) {
	open SH,'>',"$opath${lib}_soap.cmd";
	print SH join("\n",@{$cmdlines{$lib}}),"\n";
	close SH;
	open SH,'>',"$opath${lib}_soap.sh";
	print SH "#!/bin/sh
#\$ -N \"soap_$lib\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=4.1G,s_core=5
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-",scalar @{$cmdlines{$lib}},"
SEEDFILE=$opath${lib}_soap.cmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/dosoap.pl \$SEED
";
	close SH;
}

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
