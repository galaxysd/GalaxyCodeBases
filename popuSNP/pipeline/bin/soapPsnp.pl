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

our $opts='i:o:n:g:v:r:c:f:bqd';
our($opt_i, $opt_n, $opt_g, $opt_r, $opt_c, $opt_f, $opt_o, $opt_v, $opt_b, $opt_q, $opt_d);

our $desc='1.filter fq, 2.stats';
our $help=<<EOH;
\t-i fqs.nfo path (./1fqfilted) with fqs.nfo
\t-n fqs.nfo name (fqs.nfo)
\t-r Reference Genome for Soap2 (./Ref) with *.index.bwt
\t-c Chromosome NFO file (chr.nfo) in format: /^ChrID\\s+ChrLen\\s?.*\$/
\t-f faByChr path (./faByChr) with ChrID.fa\(s\)
\t-o Output path (.) for [./2soap, ./3GLF, ./4pSNP], will mkdir if not exist
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
$opt_c='chr.nfo' if ! $opt_c;
$opt_f='./faByChr' if ! $opt_f;
$opt_o='.' if ! $opt_o;

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
$opt_f=~s#/+$##;
$opt_o=~s#/+$##;
# `readlink -f` will be blank if target not exists.
system('mkdir','-p',$opt_o);

my $nfoname=$opt_i.'/'.$opt_n;
die "[x]-i $nfoname not exists !\n" unless -f $nfoname;

print STDERR "From [$opt_i]/[$opt_n] to [$opt_o] refer to [$opt_r] as [$ref], then [$opt_c][$opt_f]\n";
print STDERR "Soap2 -g [$opt_g]\n" if $opt_g;
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my (@ChrIDs,%ChrLen);
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

my ($maxRL,%DATrbrf,%SampleRL,$readlen,$min,$max)=(0);
my ($opath,%cmdlines,%Lanes);
open LST,'<',$nfoname or die "[x]Error opening $nfoname: $!\n";
while (<LST>) {
	chomp;
	my ($PESE,$sample,$lib,$FL,$ins,$ext,$path,@fqs)=split /\t/;
	($readlen,$min,$max)=split ',',$ins;
	next if $readlen == -1;
	$SampleRL{$sample} = $readlen if (! $SampleRL{$sample}) or $SampleRL{$sample} < $readlen;
	$maxRL = $readlen if$maxRL < $readlen;
	#$Librmx{$sample}{$lib}{$FL}=[$readlen,$min,$max];	# readlen,minIS,maxIS
	#$FQnfo{$sample}{$lib}{$FL}=[$ext,$path,@fqs];
	$opath="$opt_o/2soap/$sample/$lib/";
	system('mkdir','-p',$opath);
	push @{$Lanes{$sample}{$lib}},[$PESE,$opath.$fqs[0],$FL,$readlen];	# without ext
	push @{$cmdlines{$sample}},"$PESE $ins $opt_g '$opath' '$ref' '$ext,$path' '".join(',',@fqs)."'";
}
close LST;
if ($opt_v > 3) {
	print '-' x 80,"\n";
	for my $lib (keys %cmdlines) {
		print "[$lib]\n";
		print "[$_]\n" for @{$cmdlines{$lib}};
		print "\n";
	}
	print '-' x 80,"\n";
}

### 2soap ###
my (%LMSlist,%MergeOut,%LMScmdlines,%cmdlinesMerged,%mergedbychr);
$opath="$opt_o/2soap/list/";
system('mkdir','-p',$opath);
#system('mkdir','-p',"$opt_o/2soap/megred/lst/");

open SOAPL,'>',"$opt_o/2soap/soaps.lst";
open LST,'>',"$opt_o/2soap/megred.lst";
for my $sample (sort keys %Lanes) {
	my $i=1;
	#$MergedbySample{$sample}="$opt_o/2soap/megred/$sample/";
	system('mkdir','-p',"$opt_o/2soap/megred/$sample/");	# ${sample}_ChrID.sp
	system('mkdir','-p',"$opt_o/2soap/$sample/megre/");
	for my $lib (keys %{$Lanes{$sample}}) {
		my $listname="$opath${lib}.lst";
		my $mergeprefix="$opt_o/2soap/$sample/megre/$lib";	# $lib.ChrID. Planning to be $lib_ChrID.lms	lib-merged soap
		#push @{$MergeList{$sample}},$listname;
		push @{$MergeOut{$sample}},$mergeprefix;
		open L,'>',$listname;
		for (@{$Lanes{$sample}{$lib}}) {
			my ($PESE,$soapfp,$FL,$readlen)=@$_;
			if ($PESE eq 'PE') {
				print L "PE\t$i\t$soapfp.soap\n"; ++$i;
				print L "SE\t$i\t$soapfp.single\n"; ++$i;
			} else {
				print L "SE\t$i\t$soapfp.se\n"; ++$i;
			}
			print SOAPL join("\t",$PESE,$sample,$lib,$FL,$readlen,$soapfp.'.nfo'),"\n";
		}
		close L;
		push @{$LMScmdlines{$sample}},"-bi $listname -c $_ -o $mergeprefix >$mergeprefix.$_.log 2>$mergeprefix.$_.err" for @ChrIDs;
	}
	for my $chrid (@ChrIDs) {
		open L,'>',"${opath}/${sample}_$chrid.lmslst";
		print L "${_}.$chrid\n" for @{$MergeOut{$sample}};	# ${_}_$chrid.lms
		close L;
		my $spname="$opt_o/2soap/megred/$sample/${sample}_$chrid.sp";
		print LST "$sample\t$chrid\t$SampleRL{$sample}\t$spname\n";
		push @{$mergedbychr{$chrid}},$spname;	# $spname
		push @{$cmdlinesMerged{$sample}},"${opath}${sample}_$chrid.lmslst $opt_o/2soap/megred/$sample/${sample}_$chrid"
	}
}
close LST;
close SOAPL;

$opath="$opt_o/2soap/sh/";
system('mkdir','-p',$opath);
for my $sample (keys %cmdlines) {
	open SH,'>',"$opath${sample}_soap.cmd";
	print SH join("\n",@{$cmdlines{$sample}}),"\n";
	close SH;
	open SH,'>',"$opath${sample}_soap.sh";
	print SH "#!/bin/sh
#\$ -N \"sp_$sample\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=4.1G,s_core=5
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-",scalar @{$cmdlines{$sample}},"
SEEDFILE=$opath${sample}_soap.cmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/dosoap.pl \$SEED
";
	close SH;

	open SH,'>',"$opath${sample}_lms.cmd";
	print SH join("\n",@{$LMScmdlines{$sample}}),"\n";
	close SH;
	open SH,'>',"$opath${sample}_lms.sh";
	print SH "#!/bin/sh
#\$ -N \"rd_$sample\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=3.6G
#\$ -hold_jid \"sp_$sample\"
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-",scalar @{$LMScmdlines{$sample}},"
SEEDFILE=$opath${sample}_lms.cmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/rmdupbylib.pl -d \$SEED
";
	close SH;

	open SH,'>',"$opath${sample}_merge.cmd";
	print SH join("\n",@{$cmdlinesMerged{$sample}}),"\n";
	close SH;
	open SH,'>',"$opath${sample}_merge.sh";
	print SH "#!/bin/sh
#\$ -N \"mg_$sample\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=12M
#\$ -hold_jid \"rd_$sample\"
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-",scalar @{$cmdlinesMerged{$sample}},"
SEEDFILE=$opath${sample}_merge.cmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/merge.pl \$SEED
";
	close SH;

}

### 3GLF ###
$opath="$opt_o/3GLF/";
system('mkdir','-p',$opath.'sh');
system('mkdir','-p',$opath.'matrix');
## Matrix
if (-s "${opath}matrix/all.matrix") {
	system("mv -f ${opath}sh/all_matrix.sh ${opath}sh/all_matrix.oldsh") if (-e "${opath}sh/all_matrix.sh");
} else {
	open SH,'>',$opath."sh/all_matrix.sh" || die "$!\n";
	print SH "#!/bin/sh
#\$ -N \"All_Matrix\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=266M
#\$ -hold_jid \"mg_*\"
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash
#f=`find $opt_o/2soap/megred/ -name '*.sp'|xargs ls -lH|awk '{print \$5,\$9}'|sort -nrk1|head -n1|awk '{print \$2}'`
perl $SCRIPTS/getmatrix.pl $opt_o/2soap/megred.lst $opt_f ${opath}matrix/all
";
	close SH;
}
## SoapSNP
for my $chrid (@ChrIDs) {
	system('mkdir','-p',$opath.$chrid);
	open L,'>',$opath.$chrid.'.glflst';
	print L "$opath$chrid/${_}_$chrid.glf\n" for sort keys %Lanes;	# $sample
	close L;
	open L,'>',$opath.$chrid."/$chrid.cmd";
	print L "$opt_o/2soap/megred/$_/${_}_$chrid.sp $opath$chrid/${_}_$chrid d\n" for keys %Lanes;	# $sample
	close L;
	open SH,'>',$opath."sh/${chrid}_glf.sh";
	print SH "#!/bin/sh
#\$ -N \"glf_$chrid\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=280M,s_core=1
#\$ -hold_jid \"All_Matrix\"
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-",scalar keys %Lanes,"
SEEDFILE=${opath}${chrid}/$chrid.cmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval perl $SCRIPTS/callglf.pl ${opath}matrix/all.matrix $opt_f $maxRL \$SEED
";
	close SH;
}


### 4pSNP ###
$opath="$opt_o/4pSNP/";
system('mkdir','-p',$opath.'sh');
for my $chrid (@ChrIDs) {
	my $dir = $opath.$chrid;
	system('mkdir','-p',$opath.$chrid);
	my $lst=$opt_o.'/3GLF/'.$chrid.'.glflst';
	my $len=$ChrLen{$chrid};

	open LST,'>',$dir.'.psnplst' || die "$!\n";
	open CMD,'>',$dir."/$chrid.popcmd" || die "$!\n";
	my $lstcount=0;
	my ($i,$j,$exit)=(1,$POPSPLIT,0);
	while ($exit != -1) {
		++$exit;
		if ($j > $len) {
			$j=$len;
			$exit=-1;
		}
		print LST "$dir/${chrid}_$exit.psnp\n";
		print CMD "$i $j $lst $dir/${chrid}_$exit.psnp >$dir/${chrid}_$exit.tag 2>$dir/${chrid}_$exit.log\n";
		++$lstcount;
		$i += $POPSPLIT;
		$j += $POPSPLIT;
	}
	close CMD;
	close LST;

	open SH,'>',$opath.'sh/'.$chrid.'_popsnp.sh' || die "$!\n";
	print SH "#!/bin/sh
#\$ -N \"pSNP_$chrid\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=70M,s_core=1
#\$ -hold_jid glf_$chrid
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-$lstcount
SEEDFILE=${dir}/$chrid.popcmd
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval python $SCRIPTS/GLFmulti.py \$SEED
";
	close SH;
	open SH,'>',$opath.'sh/'.$chrid.'_wcsnp.sh' || die "$!\n";
	print SH "#!/bin/sh
#\$ -N \"WSNP_$chrid\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=60M
#\$ -hold_jid pSNP_$chrid
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash
cat $dir.psnplst|xargs wc -l > $dir.psnpwc
rm -f $dir.psnp
cat $dir.psnplst|xargs cat >> $dir.psnp
wc -l $dir.psnp >> $dir.psnpwc
";	# whenever >> , remember to rm first !
	close SH;
}

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
