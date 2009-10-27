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

./0rawfq , the path to fq file cannot contain more than 1 "_{$LibName}", since it is searched directly on fullpath.

=cut
our $opts='i:o:s:l:c:r:bvq';
our($opt_i, $opt_s, $opt_o, $opt_v, $opt_b, $opt_c, $opt_q, $opt_r);

our $desc='1.filter fq, 2.soap, 3.rmdup';
our $help=<<EOH;
\t-i FASTAQ file path (./0rawfq) with *.adapter.list and *.fq
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
my @fq = `find $opt_i -name '*.fq'`;	# no need to sort
chomp @fq;
my (%fqfile2rawfp,%fq1,%fq2,%fqse,%fqpe);	# no ext.
#@fqfiles = map {[fileparse($_, qr/\.fq/)]} @fq;
for (@fq) {
	my ($file, $path, $ext) = fileparse($_, qr/\.fq/);
	next if $file =~ /IndexPooling\d+(_[12])?$/;
	$fqfile2rawfp{$file}=$path;	# $path is enough. But, who cares? Me!
	$file =~ /_([^_]+)(_[12])?$/;
	my $lib = $1;
	unless ($2) {
		push @{$fqse{$1}},[$file];	# well, same API is better ?
		next;
	}
	$fq1{$file}=$1 if $2 eq '_1';
	$fq2{$file}=$1 if $2 eq '_2';
}
for my $file (keys %fq1) {
	my $file2=$file;
	$file2=~s/_1$/_2/;
	my $lib=$fq1{$file};
	delete $fq1{$file};
	if (defined $fq2{$file2}) {
		push @{$fqpe{$lib}},[$file,$file2];
		delete $fq2{$file2};
	} else {
		warn "[!][$file.fq] is not paired with _2 !\n";
		push @{$fqse{$lib}},[$file];
	}
}
for my $file (keys %fq2) {
	my $lib=$fq2{$file};
	delete $fq2{$file};
	warn "[!][$file.fq] is not paired with _1 !\n";
	push @{$fqse{$lib}},[$file];
}
%fq1=%fq2=();	# useless now
my %fqbylib=%fqpe;
push @{$fqbylib{$_}},@{$fqse{$_}} for (keys %fqse);
if ($opt_v) {
	for my $k (sort keys %fqbylib) {
		#print "\n[$k]\n[",join("]\n[",@{$fqbylib{$k}}),"]\n"
		print "\n[$k]\n[";
		for (@{$fqbylib{$k}}) {
			print join(']-[',@$_);
		}
		print "]\n";
	}
}
#my ($skip,$count,$copy,$lstcount,$adapter,%copy)=(0,0,0,0);
sub callfqfilter($$$$$$$$$) {
	my ($k,$file,$dir,$adapter,$path,$skip,$count,$copy,$lstcount,$copy_ref)=@_;
	my $fq="${path}$file.fq";
# existance check, old version
	if (-s "${dir}/${file}.nfo" or -s "${dir}/${file}.fq.bz2") {	# We may package the data and mix *.nfo then?
		++$skip;
		system("mv -f ${dir}/${file}_filte.sh ${dir}/${file}_filte.oldsh") if ( (! $adapter) and -e "${dir}/${file}_filte.sh");
		return [$skip,$count,$copy,$lstcount,$copy_ref];
	}
# existance check, old version
	unless ($adapter) {
		++${$copy_ref}{$path};
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
		return [$skip,$count,$copy,$lstcount,$copy_ref];
	}
	chomp $adapter;
#print "[$dir] [$file]\t[$adapter]\n";
	system('mkdir','-p',$dir);
	my $cmd="$adapter $fq >$dir/${file}.fq 2>$dir/${file}.nfo\n";
	print OUT $cmd;
	++$lstcount;
	return [$skip,$count,$copy,$lstcount,$copy_ref];
}

my ($sample,$cmd,@sh);
my (%fqfile2fp);	# new path
for my $k (keys %fqbylib) {
	my $withPE=0;
	$sample=$LibSample{$k}->[0];
	my $dir = $opath."/$sample/$k";
	system('mkdir','-p',$dir);
	open OUT,'>',$dir.'.cmd' || die "$!\n";
	open LST,'>',$dir.'.lst' || die "$!\n";
	#open SH,'>',$dir.'_filte.sh' || die "$!\n";
	#open LOG,'>',$dir.'.log' || die "$!\n";
	my ($skip,$count,$copy,$lstcount,%copy)=(0,0,0,0);
	my $copy_ref=\%copy;
	for (@{$fqbylib{$k}}) {
		my ($fq1,$fq2)=@$_;
		++$count;
		my $path=$fqfile2rawfp{$fq1};
		my @adapter = `find $path -name '*.list'`;
		@adapter=sort @adapter;	# if _1 & _2, _1 should be first after sort.
		$fqfile2fp{$fq1}=$dir;
		($skip,$count,$copy,$lstcount,$copy_ref)=@{&callfqfilter($k,$fq1,$dir,$adapter[0],$path,$skip,$count,$copy,$lstcount,$copy_ref)};
		if ($fq2) {	# PE
			$withPE=1;
			print LST "PE\t${dir}/${fq1}.fq\t${dir}/${fq2}.fq\n";
			$fqfile2fp{$fq2}=$dir;
			my $path=$fqfile2rawfp{$fq2};
			my @adapter = `find $path -name '*.list'`;
			@adapter=sort @adapter;
			($skip,$count,$copy,$lstcount,$copy_ref)=@{&callfqfilter($k,$fq2,$dir,$adapter[-1],$path,$skip,$count,$copy,$lstcount,$copy_ref)};
		} else {	# SE
			print LST "SE\t${dir}/${fq1}.fq\n";
		}
	}
	%copy=%$copy_ref;
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
## ReadLen
	if (-s "${dir}.ReadLen") {
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
grep '# MaxReadLen' $dir/*.nfo | perl -F'\\t' -lane 'BEGIN {\$c=99999999} END { print \"\$b\\n\$c\" } \$b=\$F[-1] if \$b<\$F[-1];\$c=\$F[-1] if \$c>\$F[-1];' > $dir.ReadLen
";
# grep '# MaxReadLen' /panfs/GAG/huxuesong/panda/1fqfilted/GP1/PANwskRAADJBAPEI/*.nfo | perl -F'\t' -lane 'BEGIN {$c=999999} END { print "$b\n$c" } $b=$F[-1] if $b<$F[-1];$c=$F[-1] if $c>$F[-1];'
		close SH;
	}
## InsertSizing
open O,'>',"${dir}.insize";
print O "100\t400\n";
close O;
### DEBUE CODE ###
	if (-s "${dir}.insize") {
		system("mv -f ${dir}_insize.sh ${dir}_insize.oldsh") if (-e "${dir}_insize.sh");
	} else {
		open SH,'>',$dir.'_insize.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"size_$k\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=4.1G,s_core=4
#\$ -hold_jid len_$k
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash
perl $SCRIPTS/instsize.pl ${dir}.lst ${dir}.ReadLen $opt_r $dir
";
		close SH;
	}

}	# End for my $k (keys %fqbylib)
## Qsub
if ($opt_q) {
	print STDERR '-' x 75,"\n";
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
	@sh = `find $opath -name '*_insize.sh'`;
	chomp @sh;
	for (@sh) {
		print STDERR "[$_]\n" if $opt_v;
		system("qsub $_");
	}
}
### 2.soap
# %fqpe,%fqse
# /share/raid010/resequencing/resequencing/tmp/bin/pipeline/SNPcalling/subBin/soap2.20 -p 4 -a 090811_I58_FC42C7AAAXX_L8_PANwkgRBMDXAAPEI_1.fq -b 090811_I58_FC42C7AAAXX_L8_PANwkgRBMDXAAPEI_2.fq -D Panda.merge.fa.index -o 090811_I58_FC42C7AAAXX_L8_PANwkgRBMDXAAPEI_1.fq.soap -2 090811_I58_FC42C7AAAXX_L8_PANwkgRBMDXAAPEI_1.fq.single -m 100 -x 400  -t -s 40 -l 32 -v 3
my $lastopath=$opath;
$opath=$opt_o.'/2soap';
system('mkdir','-p',$opath);
for my $k (keys %fqse) {
	$sample=$LibSample{$k}->[0];
	my $dir = $opath."/$sample/$k";
	system('mkdir','-p',$dir);
	open OUT,'>',$dir.'.secmd' || die "$!\n";
	my $lstcount=0;
	for (@{$fqse{$k}}) {
		my $fq=$$_[0];
		my $path=$lastopath."/$sample/$k";
		unless (-s "${dir}${fq}.nfo" or -s "${dir}${fq}.se.bz2") {
			print OUT "${path}.ReadLen ${path}/${fq}.fq $opt_r $dir/$fq\n";
			++$lstcount;
		}
	}
	close OUT;
	if ($lstcount > 0) {
		open SH,'>',$dir.'_soapse.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"se_$k\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=4.1G,s_core=4
#\$ -hold_jid len_$k
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-$lstcount
SEEDFILE=${dir}.secmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/soapse.pl \$SEED
";
		close SH;
	} else {
		unlink $dir.'.secmd';
		unlink $dir.'_soapse.sh';
	}
}
for my $k (keys %fqpe) {
	$sample=$LibSample{$k}->[0];
	my $dir = $opath."/$sample/$k";
	system('mkdir','-p',$dir);
	open OUT,'>',$dir.'.pecmd' || die "$!\n";
	my $lstcount=0;
	for (@{$fqpe{$k}}) {
		my ($fq1,$fq2)=@$_;
		my $path=$lastopath."/$sample/$k";
		unless (-s "${dir}${fq1}.nfo" or -s "${dir}${fq1}.soap.bz2") {
			print OUT "${path}.insize ${path}.ReadLen ${path}/$fq1.fq ${path}/$fq2.fq $opt_r $dir/$fq1\n" ;
			++$lstcount;
		}
	}
	close OUT;
	if ($lstcount > 0) {
		open SH,'>',$dir.'_soappe.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"pe_$k\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=4.1G,s_core=4
#\$ -hold_jid size_$k
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-$lstcount
SEEDFILE=${dir}.pecmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/soappe.pl \$SEED
";
		close SH;
	} else {
		unlink $dir.'.pecmd';
		unlink $dir.'_soappe.sh';
	}
}
## Qsub
if ($opt_q) {
	print STDERR '-' x 75,"\n";
	@sh = `find $opath -name '*_soapse.sh'`;
	chomp @sh;
	for (@sh) {
		print STDERR "[$_]\n" if $opt_v;
		system("qsub $_");
	}
	@sh = `find $opath -name '*_soappe.sh'`;
	chomp @sh;
	for (@sh) {
		print STDERR "[$_]\n" if $opt_v;
		system("qsub $_");
	}
}


#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
#print "\nPlease use the following command to batch qsub:\033[32;1m
#find ${opt_o}/ -name '*.sh' | while read ll; do qsub -l vf=2G -cwd \$ll; done\n\033[0;0m\n";
__END__
