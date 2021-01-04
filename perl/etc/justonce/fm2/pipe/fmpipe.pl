#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180516
=cut
use strict;
use warnings;
use Cwd qw(abs_path cwd);
use Parse::CSV;
#use Data::Dump qw(ddx);
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

my @Modes = qw(CHIP PCR);
my %Mode = map { $_ => 1 } @Modes;
my @Machines = qw(BGISEQ PROTON);
my %Machine = map { $_ => 1 } @Machines;
my @Parentages = qw(DUO TRIO);
my %Parentage = map { $_ => 1 } @Parentages;

my ($theMode,$theMachine,$theParentage,$fninfo,$fnfam,$pchip,$pout) = @ARGV;
$pout = '.' unless $pout;
$pout =~ s/\/+$//g;
my $listFQ = "$pout/fq.lst";
my $listFamily = "$pout/family.lst";
our $rRefn = "ref/ForensicN.fa.gz";

my $Usage = "Usage: $0 <mode/help/example> <BGISEQ|PROTON> <Trio|Duo> <info.csv> <details.csv> <chip path> [out path]\n
[Mode]: PCR mode require 2 repeats while Chip mode require 1.
[Machine]: PROTON for bam files from IonTorrent.
[Parentage]: A duo test involves the child and an alleged father. eg. gestational/full surrogacy。
In contrast, a trio test involves the mother, child, and the alleged father.
";

my $egInfo = <<'END_EG';
UID,Sample,Cell,Lane,Index
1,HK19M241C,V300016438,L02,520
2,HK19M241F,V300016438,L02,563
3,HK19M241M,V300016438,L02,564
34,NS19M342C,V300016438,L01,555
35,NS19M342F-N,V300016438,L01,519
37,NS19M342M,V300016438,L01,520
END_EG
my $egDetail = <<'END_EG';
Sample,Name,Sex,IDCard,MobilePhone,Email
HK19M241C,曹亚旭,M,110101199108120047,13712345678,CYX@qq.com
HK19M241F,钟连杰,M,232102195812116215,13812345678,zhonglj@sina.com
HK19M241M,常淑萍,F,152530196203041362,13612345678,changsp@163.com
END_EG
my $takeoutNFO = <<"END_NFO";
#### If the chip path is not in BGISEQ schema, you must ensure Sample names are correct, and then modify [$listFQ].
#### [BGISEQ]/nanoballs fq.gz can be SE(favored) or PE, [PROTON]/IonTorrent bams are all treated as SE.
####    If you do have PE Ion Torrent data, convert to *_[12].fq.gz pairs first.
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

if (@ARGV < 4) {
	die $Usage;
}
$theMode = uc $theMode;
if (! exists $Mode{$theMode}) {
	die "[x]mode can only be:[",join(',',@Modes),"].\n";
}
$theMachine = uc $theMachine;
if (! exists $Machine{$theMachine}) {
	die "[x]machine can only be:[",join(',',@Machines),"].\n";
}
$theParentage = uc $theParentage;
if (! exists $Parentage{$theParentage}) {
	die "[x]parentage can only be:[",join(',',@Parentages),"].\n";
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

