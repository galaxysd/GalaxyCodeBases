#!/usr/bin/perl -w
#use threads;
use strict;
use warnings;
#use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;
my $bin='/share/raid010/resequencing/resequencing/tmp/resequencing/animal/panda20090805/redup/doredupmerge.pl';

our $opts='i:o:s:c:bv';
our($opt_i, $opt_s, $opt_o, $opt_v, $opt_b, $opt_c);

our $desc='SoapSort library PCR PE Duplication Remover & Merger (Qsub Edition)';
our $help=<<EOH;
\t-i SOAP result path (./soap) with ./PE and ./SE
\t-s sample.list (sample.list) in format: /^sample\\tlib\\t?.*\$/
\t-c chromosome list (chrorder) in format: /^chr_name\\n\$/
\t-o Merge output path (./merged), will mkdir if not exist
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='./soap' if ! $opt_i;
$opt_s='sample.list' if ! $opt_s;
$opt_o='./merged' if ! $opt_o;
$opt_c='chrorder' if ! $opt_c;

print STDERR "From [$opt_i] to [$opt_o] refer to [$opt_s][$opt_c]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my (%sample,$tmppath,@PEp,@PEs,@SE,@Chromosome);
opendir INDIR,$opt_i.'/PE'  or warn "[!]Error opening [${opt_i}/PE]: $!\n";
for my $file (readdir INDIR) {
	if ($file =~ /\_\d\.fq\.soap$/) {push @PEp,$file;}
	 elsif ($file =~ /\_\d\.fq\.single$/) {push @PEs,$file;}
}
closedir INDIR;
opendir INDIR,$opt_i.'/SE'  or warn "[!]Error opening [${opt_i}/SE]: $!\n";
for my $file (readdir INDIR) {
	if ($file =~ /\.fq\.soap$/) {push @SE,$file;}
}
closedir INDIR;
my @files=(@PEp,@PEs,@SE);
die "[x]No porper SOAP output files found. Quit.\n" if $#files==-1;

if ($opt_v) {
	print "[$_]\n" for @files;
}

open SAMPLE,'<',$opt_s or die "Error opening $opt_s: $!\n";
while (<SAMPLE>) {
	chomp;
	s/\r//g;
	my ($sample,$lib)=split /\t/;
	print "[$sample]\t[$lib]\n" if $opt_v;
	push @{$sample{$sample}},$lib;
}
close SAMPLE;
open CHR,'<',$opt_c or die "Error opening $opt_c: $!\n";
@Chromosome=<CHR>;
chomp @Chromosome;
close CHR;
system('mkdir','-p',$opt_o);

for my $sample (keys %sample) {
	$tmppath=$opt_o.'/'.$sample;
	mkdir $tmppath,0775;
	open OUT,'>',$tmppath.'/'.$sample.'_lst.soaplist';
	my $filecount=1;
	for my $lib (@{$sample{$sample}}) {
		for (grep /$lib/,@PEp ) {
			print OUT "PE\t${opt_i}/PE/$_\t$filecount\n";
			++$filecount;
		}
		for (grep /$lib/,@PEs ) {
			print OUT "SE\t${opt_i}/PE/$_\t$filecount\n";
			++$filecount;
		}
		for (grep /$lib/,@SE ) {
			print OUT "SE\t${opt_i}/SE/$_\t$filecount\n";
			++$filecount;
		}
	}
	close OUT;
	for (@Chromosome) {
		open OUT,'>',$tmppath.'/'.$sample.'_'.$_.'.sh';
		print OUT "#!/bin/bash\nsource /home/huxuesong/works/run.sh\n
echo running \@ \`hostname\` >'${tmppath}/${sample}_${_}.err'
date >>'${tmppath}/${sample}_${_}.err'
perl '$bin' -b -c '$_' -i '${tmppath}/${sample}_lst.soaplist' -o '${tmppath}/$sample' >'${tmppath}/${sample}_${_}.log' 2>>'${tmppath}/${sample}_${_}.err'
date >>'${tmppath}/${sample}_${_}.err'\n";
#"	# to balance the poor highlighting of gedit, >_<
		close OUT;
	}
}

# Well, we can auto qsub with caltulated vf. Adding soon, maybe.

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
print "\nPlease use the following command to batch qsub:\033[32;1m
find ${opt_o}/ -name '*.sh' | while read ll; do qsub -l vf=2G -cwd \$ll; done\n\033[0;0m\n";
__END__
GP26	PANwkgRATDXAAPEI	75
GP30	PANwkgRAXDXAAPEI	75
