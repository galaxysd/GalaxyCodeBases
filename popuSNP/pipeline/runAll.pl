#!/bin/env perl
use strict;
use warnings;
use lib '/share/raid010/resequencing/soft/lib';
use Time::HiRes qw ( gettimeofday tv_interval );
use File::Basename;
use Galaxy::ShowHelp;

$main::VERSION=0.0.3;
my $SCRIPTS='/panfs/GAG/huxuesong/scripts';
my $POPSPLIT=1000000;

=pod
Note:
If a PE Lib come with only one .adapter.list file, all will use that one. If more than 2, only (sort)[0,-1] is used as _1 and _2.

./0rawfq , the path to fq file cannot contain more than 1 "_{$LibName}", since it is searched directly on fullpath.

=cut
our $opts='i:o:s:l:c:r:f:m:bvqd';
our($opt_i, $opt_s, $opt_o, $opt_m, $opt_v, $opt_b, $opt_c, $opt_q, $opt_f, $opt_r, $opt_d);

our $desc='1.filter fq, 2.soap, 3.rmdup';
our $help=<<EOH;
\t-i FASTAQ file path (./0rawfq) with *.adapter.list and *.fq
\t-s Sample list (sample.lst) in format: /^Sample\\tLib\\tReadLen\\tInsertSize\$/
\t-c Chromosome length list (chr.len) in format: /^ChrName\\s+ChrLen\\s?.*\$/
\t-m monoploid Chromosome names ('') in format: 'ChrID1,ChrID2'
\t-r Reference Genome for Soap2 (./Ref) with *.index.bwt
\t-f faByChr path (./faByChr) with ChrID.fa\(s\)
\t-o Project output path (.), will mkdir if not exist
\t-q run qsub automatically
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
\t-d Debug mode, for test only
EOH

ShowHelp();

$opt_i='./0rawfq' if ! $opt_i;
$opt_s='sample.lst' if ! $opt_s;
$opt_o='.' if ! $opt_o;
$opt_c='chr.len' if ! $opt_c;
$opt_r='./Ref' if ! $opt_r;
$opt_f='./faByChr' if ! $opt_f;

# `readlink -f` will be blank if target not exists.
system('mkdir','-p',$opt_o);
my ($lopt_i,$lopt_r,$lopt_f);

$lopt_i=`readlink -nf $opt_i`;
$opt_o=`readlink -nf $opt_o`;
$lopt_r=`readlink -nf $opt_r`;
$lopt_f=`readlink -nf $opt_f`;

warn "[x]-i $opt_i not exists !\n" unless $lopt_i;
warn "[x]-r $opt_r not exists !\n" unless $lopt_r;
warn "[x]-f $opt_f not exists !\n" unless $lopt_f;
die "\n" unless $lopt_i and $lopt_r and $lopt_f;

my @t=`find $lopt_r -name '*.index.bwt'`;
$t[0] =~ /(.+\.index)\.\w+$/;
$lopt_r = $1;
my %Monoploid;
if ($opt_m) {
	my @Monoploid=split /,/,$opt_m;
	++$Monoploid{$_} for @Monoploid;
}

print STDERR "From [$lopt_i] to [$opt_o] refer to [$opt_s][$opt_c]\nRef:[$lopt_r][$lopt_f]\n";
print STDERR "Monoploid Chr(s):[",join(',',sort keys %Monoploid),"]\n" if %Monoploid;
print STDERR "DEBUG Mode ON !\n" if $opt_d;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my (%SampleLib,%LibSample,%SampleMaxReadLen,%LibInsSize,%ChrLen,%Info);
open SAMPLE,'<',$opt_s or die "Error opening $opt_s: $!\n";
while (<SAMPLE>) {
	chomp;
	s/\r//g;
	my ($sample,$lib,$len,$min,$max)=split /\t/;
	unless (defined $max) {
		$max=int(0.5+$min*1.08);
		$min=int(0.5+$min*0.92);
	}
	#print "Sample[$sample]\tLib[$lib]\n" if $opt_v;
	push @{$SampleLib{$sample}},$lib;
	#$LibSample{$lib}=$sample;
	$Info{$lib}=[$len,$min,$max];
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
my @fq = `find $lopt_i -name '*.fq'`;	# no need to sort
chomp @fq;
my (%fqfile2rawfp,%fq1,%fq2,%fqse,%fqpe);	# no ext.
#@fqfiles = map {[fileparse($_, qr/\.fq/)]} @fq;
for (@fq) {
	my ($file, $path, $ext) = fileparse($_, qr/\.fq/);
	next if $file =~ /IndexPooling\d+(_[12])?$/i;
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
my %fqbylib;#=%fqpe; If just copy, the value of hash will be the same pointer, thus conflict.
push @{$fqbylib{$_}},@{$fqpe{$_}} for (keys %fqpe);
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
sub qsub($$) {
	my ($opath,$file)=@_;
	my @sh = `find $opath/ -name '$file'`;	# $opath may be symlink
	chomp @sh;
	for (@sh) {
		print STDERR "[$_]\n" if $opt_v;
		system("qsub $_");
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
	next unless $sample;
	my ($len,$min,$max)=@{$Info{$k}};
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
			if ($adapter[0]) {
				($skip,$count,$copy,$lstcount,$copy_ref)=@{&callfqfilter($k,$fq2,$dir,$adapter[-1],$path,$skip,$count,$copy,$lstcount,$copy_ref)};
			}	else {
				($skip,$count,$copy,$lstcount,$copy_ref)=@{&callfqfilter($k,$fq2,$dir,$adapter[0],$path,$skip,$count,$copy,$lstcount,$copy_ref)};
				# in case no list file, @adapter will be empty, '-1' will die "Modification of non-creatable array value attempted, subscript -1"
			}
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
print O "$min\t$max\n";
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
perl $SCRIPTS/instsize.pl ${dir}.lst ${dir}.ReadLen $lopt_r $dir
";
		close SH;
	}

}	# End for my $k (keys %fqbylib)
## Qsub
if ($opt_q) {
	print STDERR '-' x 75,"\n";
	&qsub($opath,'*_filte.sh');
	&qsub($opath,'*_redlen.sh');
	&qsub($opath,'*_insize.sh');
}
### 2.soap
# %fqpe,%fqse
# /share/raid010/resequencing/resequencing/tmp/bin/pipeline/SNPcalling/subBin/soap2.20 -p 4 -a 090811_I58_FC42C7AAAXX_L8_PANwkgRBMDXAAPEI_1.fq -b 090811_I58_FC42C7AAAXX_L8_PANwkgRBMDXAAPEI_2.fq -D Panda.merge.fa.index -o 090811_I58_FC42C7AAAXX_L8_PANwkgRBMDXAAPEI_1.fq.soap -2 090811_I58_FC42C7AAAXX_L8_PANwkgRBMDXAAPEI_1.fq.single -m 100 -x 400  -t -s 40 -l 32 -v 3
my $lastopath=$opath;
$opath=$opt_o.'/2soap';
system('mkdir','-p',$opath);
my %SoapCount;
system('touch',"$opath/_.soaplst");
system("find $opath/ -name '*.soaplst' | xargs rm");	# rm dies unless input
for my $k (keys %fqse) {
	$sample=$LibSample{$k}->[0];
	next unless $sample;
	my $dir = $opath."/$sample/$k";
	system('mkdir','-p',$dir);
	open LST,'>>',$dir.'.soaplst' || die "$!\n";
	open OUT,'>',$dir.'.secmd' || die "$!\n";
	my $lstcount=0;
	for (@{$fqse{$k}}) {
		my $fq=$$_[0];
		my $path=$lastopath."/$sample/$k";
		++$SoapCount{$sample};
		print LST "SE\t",$SoapCount{$sample},"\t$dir/$fq.se\n";
		unless (-s "${dir}/${fq}.nfo" or -s "${dir}/${fq}.se.bz2") {
			print OUT "${path}.ReadLen ${path}/${fq}.fq $lopt_r $dir/$fq\n";
			++$lstcount;
		}
	}
	close OUT;
	close LST;
	if ($lstcount > 0) {
		open SH,'>',$dir.'_soapse.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"se_$k\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=4.1G,s_core=5
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
	next unless $sample;
	my $dir = $opath."/$sample/$k";
	system('mkdir','-p',$dir);
	open LST,'>>',$dir.'.soaplst' || die "$!\n";
	open OUT,'>',$dir.'.pecmd' || die "$!\n";
	my $lstcount=0;
	for (@{$fqpe{$k}}) {
		my ($fq1,$fq2)=@$_;
		my $path=$lastopath."/$sample/$k";
		++$SoapCount{$sample};
		print LST "PE\t",$SoapCount{$sample},"\t$dir/$fq1.soap\n";
		++$SoapCount{$sample};
		print LST "SE\t",$SoapCount{$sample},"\t$dir/$fq1.single\n";
		unless (-s "${dir}/${fq1}.nfo" or -s "${dir}/${fq1}.soap.bz2") {
			print OUT "${path}.insize ${path}.ReadLen ${path}/$fq1.fq ${path}/$fq2.fq $lopt_r $dir/$fq1\n";
			++$lstcount;
#die $fq1 unless $fq2;
		}
	}
	close OUT;
	if ($lstcount > 0) {
		open SH,'>',$dir.'_soappe.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"pe_$k\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=4.1G,s_core=5
#\$ -hold_jid len_$k,size_$k
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
	&qsub($opath,'*_soapse.sh');
	&qsub($opath,'*_soappe.sh');
}
### 3.rmdupmerge
$lastopath=$opath;
$opath=$opt_o.'/3rmdupmerge';
system('mkdir','-p',$opath);
system('touch',"$opath/_.mglst");
system("find $opath/ -name '*.mglst' | xargs rm");
for my $k (keys %fqbylib) {
	$sample=$LibSample{$k}->[0];
	next unless $sample;
	my $dir = $opath."/$sample/$k";
	system('mkdir','-p',$dir);
	system('mkdir','-p',"$opath/$sample/lst");
	open CMD,'>',$dir.'.rdcmd' || die "$!\n";
	my $lstcount=0;
	for my $chr (keys %ChrLen) {
		open LST,'>>',$opath."/$sample/lst/${sample}_${chr}.mglst" || die "$!\n";
		print LST "$dir/$k.$chr\n" or warn "${sample}_${chr}";
		unless (-s "$dir/$k.$chr" and -s "$dir/${chr}_$k.log") {
			print CMD "-bi $lastopath/$sample/$k.soaplst -c $chr -o $dir/$k >$dir/${chr}_$k.log 2>$dir/${chr}_$k.err\n";
			++$lstcount;
		}
		#close LST;	# Let the open close it.
	}
	close CMD;
	close LST;
	if ($lstcount > 0) {
		open SH,'>',$dir.'_rmdup.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"rmdup_${sample}_$k\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=280M
#\$ -hold_jid \"pe_$k,se_$k\"
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-$lstcount
SEEDFILE=${dir}.rdcmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/rmdupbylib.pl -d \$SEED
";	# add -d for debug
		close SH;
	} else {
		unlink $dir.'.rdcmd';
		unlink $dir.'_rmdup.sh';
	}
}

for my $k (keys %SampleLib) {
	my $dir = $opath."/$k";
	my $lib=$SampleLib{$k}->[0];
	next unless defined $fqbylib{$lib};
	open CMD,'>',$dir."/$k.mgcmd" || die "$!\n";
	my $lstcount=0;
	for my $chr (keys %ChrLen) {
		unless (-s "$dir/${k}_${chr}.sp" and -s "$dir/${k}_${chr}.log") {
			print CMD "$dir/lst/${k}_${chr}.mglst $dir/${k}_${chr}\n";
			++$lstcount;
		}
	}
	close CMD;
	if ($lstcount > 0) {
		open SH,'>',$dir."/${k}_merge.sh" || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"merge_$k\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=12M
#\$ -hold_jid \"rmdup_${k}_*\"
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-$lstcount
SEEDFILE=${dir}/$k.mgcmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/merge.pl \$SEED
";
		close SH;
	} else {
		unlink $dir."/${k}.mgcmd";
		unlink $dir."/${k}_merge.sh";
	}
}
## Qsub
if ($opt_q) {
	print STDERR '-' x 75,"\n";
	&qsub($opath,'*_rmdup.sh');
	&qsub($opath,'*_merge.sh');
}
### 4.GLF
$lastopath=$opath;
$opath=$opt_o.'/4GLF';
my $dir = $opath.'/matrix';
system('mkdir','-p',$dir);
## Matrix
if (-s "$dir/all.matrix") {
	system("mv -f ${dir}/all_matrix.sh ${dir}/all_matrix.oldsh") if (-e "${dir}/all_matrix.sh");
} else {
	open SH,'>',$dir."/all_matrix.sh" || die "$!\n";
	print SH "#!/bin/sh
#\$ -N \"All_Matrix\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=266M
#\$ -hold_jid \"merge_*\"
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash
f=`find $lastopath/ -name '*.sp'|xargs ls -lH|awk '{print \$5,\$9}'|sort -nrk1|head -n1|awk '{print \$2}'`
perl $SCRIPTS/matrix.pl \$f $lopt_f $opt_o/1fqfilted $dir/all
";
	close SH;
}
## Qsub
if ($opt_q) {
	print STDERR '-' x 75,"\n";
	&qsub($dir,'all_matrix.sh');
}
## SoapSNP
for my $chr (keys %ChrLen) {
	my $dir = $opath."/$chr";
	system('mkdir','-p',$dir);
unless ($opt_d) {
	open LST,'>',$dir.'.glflst' || die "$!\n";
	open CMD,'>',$dir."/$chr.glfcmd" || die "$!\n";
	my $lstcount=0;
	for my $k (keys %SampleLib) {
		print LST "$dir/${k}_$chr.glf\n";
		unless (-s "$dir/${k}_$chr.tag" or -s "$dir/${k}_$chr.glf") {	# 'tag' not working since soapsnp return non-0 on exit
			my $Mono=' d';
			$Mono=' m' if $Monoploid{$chr};
			print CMD "$lastopath/$k/${k}_$chr.sp $dir/${k}_${chr}$Mono\n";
			++$lstcount;
		}
	}
	close CMD;
	close LST;
	if ($lstcount > 0) {
		open SH,'>',$dir.'_glf.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"GLF_$chr\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=280M,s_core=1
#\$ -hold_jid \"All_Matrix\"
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-$lstcount
SEEDFILE=${dir}/$chr.glfcmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/callglf.pl $opath/matrix/all.matrix $lopt_f $opt_o/1fqfilted \$SEED
";
		close SH;
	} else {
		unlink $dir."/$chr.glfcmd";
		unlink $dir.'_glf.sh';
	}
}	# unless ($opt_d)
}
## Qsub
if ($opt_q) {
	print STDERR '-' x 75,"\n";
	&qsub($opath,'*_glf.sh');
}

### 5.popSNP
$lastopath=$opath;
$opath=$opt_o.'/5popSNP';
$dir = $opath.'/all';
system('mkdir','-p',$dir);
for my $chr (keys %ChrLen) {
	my $dir = $opath."/$chr";
	system('mkdir','-p',$dir);
	my $lst=$lastopath."/$chr.glflst";
	my $len=$ChrLen{$chr};
	open LST,'>',$dir.'.psnplst' || die "$!\n";
	open CMD,'>',$dir."/$chr.popcmd" || die "$!\n";
	my $lstcount=0;
	my ($i,$j,$exit)=(1,$POPSPLIT,0);
	while ($exit != -1) {
		++$exit;
		if ($j > $len) {
			$j=$len;
			$exit=-1;
		}
		print LST "$dir/${chr}_$exit.psnp\n";
		unless (-s "$dir/${chr}_$exit.tag") {
			print CMD "$i $j $lst $dir/${chr}_$exit.psnp >$dir/${chr}_$exit.tag 2>$dir/${chr}_$exit.log\n";
			++$lstcount;
		}
		$i += $POPSPLIT;
		$j += $POPSPLIT;
	}
	close CMD;
	close LST;
	if ($lstcount > 0) {
		open SH,'>',$dir.'_popsnp.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"PSNP_$chr\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=70M,s_core=1
#\$ -hold_jid GLF_$chr
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-$lstcount
SEEDFILE=${dir}/$chr.popcmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval python $SCRIPTS/GLFmulti.py \$SEED
";
		close SH;
		open SH,'>',$dir.'_wcsnp.sh' || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"WSNP_$chr\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=60M
#\$ -hold_jid PSNP_$chr
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash
cat $dir.psnplst|xargs wc -l > $dir.psnpwc
rm -f $dir.psnp
cat $dir.psnplst|xargs cat >> $dir.psnp
wc -l $dir.psnp >> $dir.psnpwc
";	# whenever >> , remember to rm first !
		close SH;
	} else {
		unlink $dir."/$chr.popcmd";
		unlink $dir.'_popsnp.sh';
		unlink $dir.'_wcsnp.sh';
	}
}
## Qsub
if ($opt_q) {
	print STDERR '-' x 75,"\n";
	&qsub($opath,'*_popsnp.sh');
	&qsub($opath,'*_wcsnp.sh');
}

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
#print "\nPlease use the following command to batch qsub:\033[32;1m
#find ${opt_o}/ -name '*.sh' | while read ll; do qsub -l vf=2G -cwd \$ll; done\n\033[0;0m\n";
__END__
