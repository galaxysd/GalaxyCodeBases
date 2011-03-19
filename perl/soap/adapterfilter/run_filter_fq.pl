#!/usr/bin/perl -w
#use threads;
use strict;
use warnings;
#use DBI;
use File::Basename;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.2.0;
my $bin = "/share/raid010/resequencing/resequencing/human/Danmark_human/fq_filter/filter_fq.pl";
$bin='/share/raid010/resequencing/resequencing/tmp/resequencing/animal/panda20090805/filter_fq_pe.pl';

our $opts='i:o:s:d:bv';
our ($opt_i, $opt_o, $opt_v, $opt_b);

our $help=<<EOH;
\t-i raw fq dir (PE is OK) (./rawindexes)
\t-o output dir (./fqfilted)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='./rawindexes' if ! defined $opt_i;
$opt_o='./fqfilted' if ! $opt_o;
$opt_i =~ s/\/$//;	# There will be an error if $opt_i with /$ and $opt_o
$opt_o =~ s/\/$//;	# without /$. No cycle as nobody write ./out// .

print STDERR "From [$opt_i] to [$opt_o]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
my @fq = `find $opt_i -name '*.fq'`;
chomp @fq;

my ($skip,$count,$copy,$adapter,%copy)=(0,0,0);
system('mkdir','-p',"${opt_o}/sh");
system('mkdir','-p',"${opt_o}/sh_old");
system("mv -f ${opt_o}/sh/* ${opt_o}/sh_old/");
for my $fq (@fq) {
	++$count;
	my ($file, $dir, $ext) = fileparse($fq, qr/\.fq/);
	my @adapter = `find $dir -name '*.list'`;
	$dir =~ s/^$opt_i/$opt_o/;
	if (-s "${dir}/${file}.fq" or -s "${dir}/${file}.fq.bz2") {
		++$skip;
		next;
	}
	if ($#adapter<0) {
		++$copy{$dir};
		++$copy;
		system('mkdir','-p',$dir);
		system('cp','-a',$fq,"${dir}/${file}.fq");
		next;
	}
	chomp @adapter;
	@adapter=sort @adapter;	# if _1 & _2, _1 should be first after sort.
#print "$_\n" for @adapter;
#print "\n";
	if ($#adapter>=1 and $file =~ /_1$/) {$adapter=$adapter[0];}
	 elsif ($#adapter>=1 and $file =~ /_2$/) {$adapter=$adapter[-1];}	# more than 2 ? last one !
	 else {$adapter=$adapter[0];}	# SE
#print "[$dir] [$file]\t[$adapter]\n";
	system('mkdir','-p',$dir);
	open OUT, '>',"${opt_o}/sh/filter_${file}.sh" || die "$!\n";
	print OUT "#!/bin/sh\n#\$ -S /bin/bash\nperl $bin $adapter $fq > ${dir}/${file}.fq\n";
	close OUT;
}

print "Of $count file(s),$skip skipped, $copy with no adapter.list.\n";
if ($copy) {
	print "\nThe following path with no adapter.list\n";
	for (sort keys %copy) {
		print $_,"\t",$copy{$_}," time(s)\n";
	}
}

#END
my $stop_time = [gettimeofday];

$|=1;
print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";

print "Please run the command below to submit jobs:\n\033[32;1m\nfor i in `ls ${opt_o}/sh/*.sh`; do qsub -cwd -l vf=1g \$i; done\n\033[0;0m\n";
__END__
for i in `ls *.sh`; do qsub -cwd -l vf=1g $i; done
