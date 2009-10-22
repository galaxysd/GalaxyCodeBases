#!/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use File::Basename;
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;
my $SCRIPTS='/panfs/GAG/huxuesong/scripts';

=pod
Note:
If a PE Lib come with only one .adapter.list file, all will use that one. If more than 2, only (sort)[0,-1] is used as _1 and _2.

./0rawfq , the path to fq file cannot contain more than 1 Lib name, since Lib is searched directly on fullpath.

=cut
our $opts='i:o:s:l:c:r:bvq';
our($opt_i, $opt_s, $opt_o, $opt_v, $opt_b, $opt_c, $opt_q, $opt_r);

our $desc='1.filter fq, 2.soap, 3.rmdup';
our $help=<<EOH;
\t-i FASTAQ file path (./0rawfq) with *.adapter.list
\t-s Sample list (sample.lst) in format: /^Sample\\tLib\$/
\t-c Chromosome length list (chr.len) in format: /^ChrName\\s+ChrLen\\s?.*\$/
\t-r Reference Genome for Soap2 (./Ref) with *.index.bwt
\t-o Project output path (.), will mkdir if not exist
\t-q run qsub automatically
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='./0rawfq' if ! $opt_i;
$opt_s='sample.lst' if ! $opt_s;
$opt_o='.' if ! $opt_o;
$opt_c='chr.len' if ! $opt_c;
$opt_r='./Ref' if ! $opt_r;

$opt_i=`readlink -nf $opt_i`;
$opt_o=`readlink -nf $opt_o`;
$opt_r=`readlink -nf $opt_r`;

my @t=`find $opt_r -name '*.index.bwt'`;
$t[0] =~ /(.+\.index)\.\w+$/;
$opt_r = $1;

print STDERR "From [$opt_i] to [$opt_o] refer to [$opt_s][$opt_c]\nRef:[$opt_r]\n";
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
### 1.filter fq
$opath=$opt_o.'/1fqfilted';
system('mkdir','-p',$opath);
my @fq = `find $opt_i -name '*.fq'`;
chomp @fq;
my %fqbylib;
for my $k (keys %LibSample) {
	my @t=grep /$k/,@fq;
	$fqbylib{$k}=\@t if $#t > -1;
	#$fqbylib{$k}=[grep /$k/,@fq];
	#delete $fqbylib{$k} if $#$fqbylib{$k}==-1;
}
if ($opt_v) {
	for my $k (sort keys %fqbylib) {
		print "\n[$k]\n[",join("]\n[",@{$fqbylib{$k}}),"]\n"
	}
}
#my ($skip,$count,$copy,$lstcount,$adapter,%copy)=(0,0,0,0);
my ($sample,$cmd,@sh);
my %fqfiltedbylib;
for my $k (keys %fqbylib) {
	$sample=$LibSample{$k}->[0];
	my $dir = $opath."/$sample/$k";
	system('mkdir','-p',$dir);
	open OUT,'>',$dir.'.cmd' || die "$!\n";
	open LST,'>',$dir.'.lst' || die "$!\n";
	#open SH,'>',$dir.'_filte.sh' || die "$!\n";
	#open LOG,'>',$dir.'.log' || die "$!\n";
	my ($skip,$count,$copy,$lstcount,$adapter,%copy)=(0,0,0,0);
	for my $fq (@{$fqbylib{$k}}) {
		++$count;
		my ($file, $path, $ext) = fileparse($fq, qr/\.fq/);
		my @adapter = `find $path -name '*.list'`;
		print LST "${dir}/${file}.fq\n";
		push @{$fqfiltedbylib{$k}},"${dir}/${file}.fq";
# existance check, old version
		if (-s "${dir}/${file}.nfo" or -s "${dir}/${file}.fq.bz2") {	# We may package the data and mix *.nfo then?
			++$skip;
			system("mv -f ${dir}/${file}_filte.sh ${dir}/${file}_filte.oldsh") if ($#adapter<0 and -e "${dir}/${file}_filte.sh");
			next;
		}
# existance check, old version
		if ($#adapter<0) {
			++$copy{$path};
			++$copy;
			#system('mkdir','-p',$dir);
			system('cp','-s',$fq,"${dir}/${file}.fq");
			#warn $fq,"${dir}/${file}.fq";
			open STAT,'>',"${dir}/${file}_filte.sh" || die "$!\n";
			print STAT "#!/bin/sh
#\$ -N \"s_$k\"
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
	close LST;
	close OUT;
	open SH,'>',$dir.'_filte.sh' || die "$!\n";
	print SH "#!/bin/sh
#\$ -N \"f_$k\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=276M
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-$lstcount
SEEDFILE=${dir}.cmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/fqfilter.pl \$SEED
" if $lstcount > 0;
	close SH;
	if ($lstcount == 0) {
		unlink $dir.'_filte.sh';
		unlink $dir.'.cmd';
	}
	open LOG,'>',$dir.'.log' || die "$!\n";
	print LOG "Of $count file(s),$skip skipped, $copy with no adapter.list.\n";
	if ($copy) {
		print LOG "\nThe following path with no adapter.list\n";
		for (sort keys %copy) {
			print LOG $_,"\t",$copy{$_}," time(s)\n";
		}
	}
	close LOG;
## maxReadLen
	if (-s "${dir}.maxReadLen") {
		system("mv -f ${dir}_redlen.sh ${dir}_redlen.oldsh") if (-e "${dir}_redlen.sh");
	} else {
		open SH,'>',$dir.'_redlen.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"len_$k\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=1M
#\$ -hold_jid f_$k,s_$k
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash
grep '# MaxReadLen' $dir/*.nfo | perl -F'\\t' -lane 'END { print \$b }\$b=\$F[-1] if \$b<\$F[-1]' > $dir.maxReadLen
";
		close SH;
	}
## InsertSizing
	if (-s "${dir}.insize") {
		system("mv -f ${dir}_insize.sh ${dir}_insize.oldsh") if (-e "${dir}_insize.sh");
	} else {
		open SH,'>',$dir.'_insize.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"size_$k\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=1M
#\$ -hold_jid len_$k,f_$k,s_$k
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash

";
		close SH;
	}

}	# End for my $k (keys %fqbylib)
## Qsub
if ($opt_q) {
	@sh = `find $opath -name '*_filte.sh'`;
	chomp @sh;
	for (@sh) {
		print STDERR "[$_]\n" if $opt_v;
		system("qsub $_");
	}
	@sh = `find $opath -name '*_redlen.sh'`;
	chomp @sh;
	for (@sh) {
		print STDERR "[$_]\n" if $opt_v;
		system("qsub $_");
	}
}
## Stat



### 2.soap
# %fqfiltedbylib
# /share/raid010/resequencing/resequencing/tmp/bin/pipeline/SNPcalling/subBin/soap2.20 -p 4 -a 090811_I58_FC42C7AAAXX_L8_PANwkgRBMDXAAPEI_1.fq -b 090811_I58_FC42C7AAAXX_L8_PANwkgRBMDXAAPEI_2.fq -D Panda.merge.fa.index -o 090811_I58_FC42C7AAAXX_L8_PANwkgRBMDXAAPEI_1.fq.soap -2 090811_I58_FC42C7AAAXX_L8_PANwkgRBMDXAAPEI_1.fq.single -m 100 -x 400  -t -s 40 -l 32 -v 3




#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
#print "\nPlease use the following command to batch qsub:\033[32;1m
#find ${opt_o}/ -name '*.sh' | while read ll; do qsub -l vf=2G -cwd \$ll; done\n\033[0;0m\n";
__END__
