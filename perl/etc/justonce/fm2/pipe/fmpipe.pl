#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180516
=cut
use strict;
use warnings;
use Cwd qw(abs_path cwd);
use Parse::CSV;
use Data::Dump qw(ddx);
use FindBin qw($RealBin);
if ($FindBin::VERSION < 1.51) {
	warn "[!]Your Perl is too old, thus there can only be ONE `fmpipe.pl` file in your PATH. [FindBin Version: $FindBin::VERSION < 1.51]\n\n"
}
FindBin::again();
use lib "$RealBin";
use FGI::fm2;
# ====================================
my %pPrefixs = (
	lst => '0lst',
	fq => '1fq',
	bam => '2bam',
	vcf	=> '3vcf',
	oyk => '4tsv',
	qc => '5qc',
);

my @Modes = qw(BGISEQ PROTON);
my %Mode = map { $_ => 1 } @Modes;


my ($theMode,$fninfo,$fnfam,$pchip,$pout) = @ARGV;
$pout = '.' unless $pout;
$pout =~ s/\/+$//g;
my $listFQ = "$pout/fq.lst";
my $listFamily = "$pout/family.lst";
our $rRefn = "ref/ForensicN.fa.gz";

my $Usage = "Usage: $0 <do/help/example> <info.csv> <details.csv> <chip path> [out path]\n
";

my $egInfo = <<'END_EG';
uid,sample,cell,lane,index
1,HK19M241C,V300083793,L02,385-331
2,HK19M241F,V300083793,L02,385-333
3,HK19M241M,V300083793,L02,386-335
34,NS19M342C,V300083793,L02,384-334
35,NS19M342F-N,V300083793,L02,470-335
37,NS19M342M,V300083793,L02,383-333
END_EG
my $egDetail = <<'END_EG';
sample,name,sex,idcard,mobilephone,email
HK19M241C,曹亚旭,M,110101199108120047,13712345678,CYX@qq.com
HK19M241F,钟连杰,M,232102195812116215,13812345678,zhonglj@sina.com
HK19M241M,常淑萍,F,152530196203041362,13612345678,changsp@163.com
END_EG
my $takeoutNFO = <<"END_NFO";
#### To skip `cutadapt`, run `SKIP=1 $pout/1fq.sh` manually WITHOUT qsub it.
#### If you then regret to use `cutadapt`, you must clean "$pout/$pPrefixs{fq}/", but KEEP the directory itself, before qsub it.
END_NFO
if (defined $theMode) {
	if ($theMode =~ /(\bh\b)|(help)/i) {
		warn $Usage,<<"END_MESSAGE";

### Example of "info.csv":
$egInfo
#### `UID` must be UNIQUE !

### Example of "details.csv":
$egDetail
#### Use Sample names in "details.csv".

#### All other columns will be ignored.

$takeoutNFO
#### You can run `$0 example` or `$0 eg` to overwrite "info.csv" and "details.csv" with examples above.

END_MESSAGE
		exit 0;
	} elsif ($theMode =~ /(\beg\b)|(example)/i) {
		open I,'>',"info.csv" or die $?;
		open F,'>',"details.csv" or die $?;
		print I $egInfo;
		print F $egDetail;
		close F; close I;
		exit 0;
	}
}

if (@ARGV < 2) {
	die $Usage;
}
$theMode = uc $theMode;
if (! exists $Mode{$theMode}) {
	die "[x]mode can only be:[",join(',',@Modes),"].\n";
}

my $cwd = cwd() or die $!;
die "[x]Cannot read info.csv [$fninfo].\n" unless -r -s $fninfo;
die "[x]Cannot read details.csv [$fnfam].\n" unless -r -s $fnfam;
die "[x]Cannot read BWA Reference [$RealBin/$rRefn.*].\n[!]Run `samtools faidx $RealBin/$rRefn` and `bwa index $RealBin/$rRefn` first.\n" unless -r -s "$RealBin/$rRefn.sa";
die "[x]Cannot read Chip Path [$pchip].\n" unless -r $pchip;
#die "[x]Cannot use out Path [$pout].\n" unless -r -w -x $pout;
$pchip =~ s/\/+$//g;
#$pout =~ s/\/+$//g;
$pchip = abs_path($pchip);
warn "[!]Info:[$fninfo], Fam:[$fnfam], CHIP:[$pchip] to Out:[$pout]\n[!]Current working directory:[$cwd]\n";

my $cinfo = Parse::CSV->new(file => $fninfo, names => 1);
my $cfam = Parse::CSV->new(file => $fnfam, names => 1);
system('mkdir','-p',$pout);
for (keys %pPrefixs) {
	mkdir "$pout/$pPrefixs{$_}";
}

################################
# fq.lst, q0cutadapter.sh, q1bwa.sh #
################################
my (%fqInfo,%Samples);
while ( my $value = $cinfo->fetch ) {
	next if $value->{cell} eq '';
	my $usid = join('_',$value->{sample},$value->{uid});
	$fqInfo{$usid} = [$value->{cell},$value->{lane},$value->{index}];
	push @{$Samples{$value->{sample}}},$usid;
	#die "[x]Column PESE must be either PE or SE.\n" unless $value->{PESE} =~ /(P|S)E/i;
}
die $cinfo->errstr if $cinfo->errstr;
#ddx \%fqInfo,\%Samples;
my (%SampleCnt,@sCnt);
for (keys %Samples) {
	++$SampleCnt{scalar @{$Samples{$_}}};
}
#ddx \%SampleCnt;
@sCnt = keys %SampleCnt;
die "[x]One fq, one sample, for now.\n" if @sCnt > 1;
open O,'>',$listFQ or die $?;
for (sort keys %fqInfo) {
	my @d = @{$fqInfo{$_}};
	my $fqNameP;
	if ($theMode eq 'BGISEQ') {
		$fqNameP = join('/',$pchip,$d[0],$d[1],join('_',$d[0],$d[1],$d[2]));
	} elsif ($theMode eq 'PROTON') {
		$fqNameP = join('/',$pchip,$d[0],$d[1],"basecaller_results",join('_',"IonXpress",$d[2]));
	}
	print O join("\t",$_,$fqNameP),"\n";
}
close O;
my $fSHcutadapt = "$pout/q0cutadapter.sh";
open O,'>',$fSHcutadapt or die $?;
my $FQprefix = "$pout/$pPrefixs{fq}";
print O Scutadapt($cwd,scalar(keys %fqInfo),$listFQ,$FQprefix);
close O;
chmod 0755,$fSHcutadapt;
my $fSHbwa = "$pout/q1bwa.sh";
open O,'>',$fSHbwa or die $?;
my $BAMprefix = "$pout/$pPrefixs{bam}";
my $VCFprefix = "$pout/$pPrefixs{vcf}";
my $OYKprefix = "$pout/$pPrefixs{oyk}";
print O Sbwamem($cwd,scalar(keys %fqInfo),$listFQ,$BAMprefix,$VCFprefix,$FQprefix,$OYKprefix);
close O;
################################



__END__
./fmpipe.pl BGISEQ info.csv details.csv intt outtt
./fmpipe.pl BGISEQ info.csv details.csv fq out
