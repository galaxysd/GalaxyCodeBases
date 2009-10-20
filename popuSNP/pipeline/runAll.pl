#!/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use File::Basename;
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;
my $SCRIPTS='/panfs/GAG/huxuesong/scripts';

our $opts='i:o:s:l:c:bv';
our($opt_i, $opt_s, $opt_o, $opt_v, $opt_b, $opt_c);

our $desc='1.filter fq, 2.soap, 3.rmdup';
our $help=<<EOH;
\t-i FASTAQ file path (./0rawfq) with *.adapter.list
\t-s Sample list (sample.lst) in format: /^Sample\\tLib\$/
\t-c Chromosome length list (chr.len) in format: /^ChrName\\s+ChrLen\\s?.*\$/
\t-o Project output path (.), will mkdir if not exist
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='./0rawfq' if ! $opt_i;
$opt_s='sample.lst' if ! $opt_s;
$opt_o='.' if ! $opt_o;
$opt_c='chr.len' if ! $opt_c;

$opt_i=`readlink -nf $opt_i`;
$opt_o=`readlink -nf $opt_o`;

print STDERR "From [$opt_i] to [$opt_o] refer to [$opt_s][$opt_c]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my (%SampleLib,%LibSample,%SampleMaxReadLen,%LibInsSize,%ChrLen);
open SAMPLE,'<',$opt_s or die "Error opening $opt_s: $!\n";
while (<SAMPLE>) {
	chomp;
	s/\r//g;
	my ($sample,$lib)=split /\t/;
	#print "Sample[$sample]\tLib[$lib]\n" if $opt_v;
	push @{$SampleLib{$sample}},$lib;
	#$LibSample{$lib}=$sample;
	push @{$LibSample{$lib}},$sample;
	$SampleMaxReadLen{$sample}=0;
}
close SAMPLE;
for (keys %LibSample) {
	my $t = scalar @{$LibSample{$_}};
	die "[!]Sample with Lib conflict($t) for Lib:[$_]. Check [$opt_s] !\n" if $t != 1;	# Check if a lib -> more samples
}
if ($opt_v) {
	for my $k (sort keys %SampleLib) {
		print "Sample:[$k]\tLib:[",join('],[',@{$SampleLib{$k}}),"]\n"
	}
	for my $k (sort keys %LibInsSize) {
		print "Lib:[$k]\tInsize:[",join(',',@{$LibInsSize{$k}}),"]\n";
	}
}
open CHRLEN,'<',$opt_c or die "Error opening $opt_c: $!\n";
while (<CHRLEN>) {
	chomp;
	s/\r//g;
	my ($chr,$len)=split /\s+/;
	#print "Chr:[$chr]\tLen:[$len]\n" if $opt_v;
	$ChrLen{$chr}=$len;
}
close CHRLEN;
if ($opt_v) {
	for my $k (sort keys %ChrLen) {
		print "Chr:[$k]\tLen:[$ChrLen{$k}]\n"
	}
}
system('mkdir','-p',$opt_o);
my $opath;
# 1.filter fq
$opath=$opt_o.'/1fqfilted';
system('mkdir','-p',$opath);
my @fq = `find $opt_i -name '*.fq'`;
chomp @fq;
my %fqbylib;
for my $k (keys %LibSample) {
	$fqbylib{$k}=[grep /$k/,@fq];
}
if ($opt_v) {
	for my $k (sort keys %fqbylib) {
		print "[$k]\n[",join("]\n[",@{$fqbylib{$k}}),"]\n"
	}
}
#my ($skip,$count,$copy,$lstcount,$adapter,%copy)=(0,0,0,0);
my ($sample,$cmd);
for my $k (keys %fqbylib) {
	$sample=$LibSample{$k}->[0];
	my $dir = $opath."/$sample/$k";
	system('mkdir','-p',$dir);
	open OUT,'>',$dir.'.lst' || die "$!\n";
	open SH,'>',$dir.'.sh' || die "$!\n";
	open NFO,'>',$dir.'.nfo' || die "$!\n";
	my ($skip,$count,$copy,$lstcount,$adapter,%copy)=(0,0,0,0);
	for my $fq (@{$fqbylib{$k}}) {
		++$count;
		my ($file, $path, $ext) = fileparse($fq, qr/\.fq/);
		my @adapter = `find $path -name '*.list'`;
# existance check, old version
		if (-s "${dir}/${file}.nfo" or -s "${dir}/${file}.fq.bz2") {	# We may package the data and mix *.nfo then?
			++$skip;
			system("mv -f ${dir}/${file}.sh ${dir}/${file}.oldsh") if ($#adapter<0 and -e "${dir}/${file}.sh");
			next;
		}
# existance check, old version
		if ($#adapter<0) {
			++$copy{$path};
			++$copy;
			#system('mkdir','-p',$dir);
			system('cp','-a',$fq,"${dir}/${file}.fq");
			#warn $fq,"${dir}/${file}.fq";
			open STAT,'>',"${dir}/${file}.sh" || die "$!\n";
			print STAT "#!/bin/sh
#\$ -N \"s$k\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=9M
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash
perl $SCRIPTS/fqstat.pl $fq 2>${dir}/${file}.nfo
";
			close STAT;
			next;
		}
		chomp @adapter;
		@adapter=sort @adapter;	# if _1 & _2, _1 should be first after sort.
		if ($#adapter>=1 and $file =~ /_1$/) {$adapter=$adapter[0];}
		 elsif ($#adapter>=1 and $file =~ /_2$/) {$adapter=$adapter[-1];}	# more than 2 ? last one !
		 else {$adapter=$adapter[0];}	# SE
#print "[$dir] [$file]\t[$adapter]\n";
		system('mkdir','-p',$dir);
		$cmd="$adapter $fq >$dir/${file}.fq 2>$dir/${file}.nfo\n";
		print OUT $cmd;
		++$lstcount;
	}
	close OUT;
	print SH "#!/bin/sh
#\$ -N \"f$k\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=276M
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-$lstcount
SEEDFILE=${dir}.lst
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/fqfilter.pl \$SEED
" if $lstcount > 0;
	close SH;
	if ($lstcount == 0) {
		unlink $dir.'.sh';
		unlink $dir.'.lst';
	}
	print NFO "Of $count file(s),$skip skipped, $copy with no adapter.list.\n";
	if ($copy) {
		print NFO "\nThe following path with no adapter.list\n";
		for (sort keys %copy) {
			print NFO $_,"\t",$copy{$_}," time(s)\n";
		}
	}
	close NFO;
}







#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
#print "\nPlease use the following command to batch qsub:\033[32;1m
#find ${opt_o}/ -name '*.sh' | while read ll; do qsub -l vf=2G -cwd \$ll; done\n\033[0;0m\n";
__END__
