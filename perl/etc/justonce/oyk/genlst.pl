#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180516
=cut
use strict;
use warnings;
use Cwd qw(abs_path cwd);
use Parse::CSV;
use lib '.';
use Data::Dump qw(ddx);
use FGI::OYK;
# ====================================

#our $pI = '5';

# ====================================

my %pPrefixs = (
	lst => '0lst',
	fq => '1fq',
	bam => '2bam',
	vcf	=> '3vcf',
	oyk => '4tsv',
);

my @Modes = qw(CHIP PCR);
my %Mode = map { $_ => 1 } @Modes;
my @Machines = qw(BGISEQ PROTON);
my %Machine = map { $_ => 1 } @Machines;

my ($theMode,$theMachine,$fninfo,$fnfam,$pRef,$pchip,$pout) = @ARGV;
$pout = '.' unless $pout;
$pout =~ s/\/+$//g;
my $listFQ = "$pout/fq.lst";
my $listFamily = "$pout/family.lst";

my $Usage = "Usage: $0 <mode/help/example> <BGISEQ|PROTON> <info.csv> <fam.csv> <Ref for BWA> <chip path> [out path]\n";

my $egInfo = <<'END_EG';
UID,Sample,Cell,Lane,Index
1,HK19M241C,V300016438,L02,520
2,HK19M241F,V300016438,L02,563
3,HK19M241M,V300016438,L02,564
34,NS19M342C,V300016438,L01,555
35,NS19M342F-N,V300016438,L01,519
37,NS19M342M,V300016438,L01,520
END_EG
my $egFam = <<'END_EG';
Father,Mother,Child
HK19M241F,HK19M241M,HK19M241C
NS19M342F-N,NS19M342M,NS19M342C
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
#### For PCR mode, each sample has two lines with identical Sample name BUT different UID.

### Example of "fam.csv":
$egFam
#### Use Sample names in "fam.csv".

#### All other columns will be ignored.

$takeoutNFO
#### If your Index is not 501~599, modify `my \$pI = '5';` in this programme.
#### You can run `$0 example` or `$0 eg` to overwrite "info.csv" and "fam.csv" with examples above.

END_MESSAGE
		exit 0;
	} elsif ($theMode =~ /(\beg\b)|(example)/i) {
		open I,'>',"info.csv" or die $?;
		open F,'>',"fam.csv" or die $?;
		print I $egInfo;
		print F $egFam;
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
my $cwd = cwd() or die $!;
die "[x]Cannot read info.csv [$fninfo].\n" unless -r -s $fninfo;
die "[x]Cannot read fam.csv [$fnfam].\n" unless -r -s $fnfam;
die "[x]Cannot read BWA Reference [$pRef].\n" unless -r -s $pRef;
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
	next if $value->{Cell} eq '';
	my $usid = join('_',$value->{Sample},$value->{UID});
	$fqInfo{$usid} = [$value->{Cell},$value->{Lane},$value->{Index}];
	push @{$Samples{$value->{Sample}}},$usid;
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
die "[x]Cannot mix Chip & PCR data.\n" if @sCnt > 1;
if ($theMode eq 'CHIP') {
	die "[x]Chip mode allows only one repeat per sample.\n" if $sCnt[0] > 1;
} elsif ($theMode eq 'PCR') {
	die "[x]PCR mode allows exactly TWO repeats per sample, got $sCnt[0].\n" if $sCnt[0] != 2;
}
open O,'>',$listFQ or die $?;
for (sort keys %fqInfo) {
	my @d = @{$fqInfo{$_}};
	my $fqNameP;
	if ($theMachine eq 'BGISEQ') {
		$fqNameP = join('/',$pchip,$d[0],$d[1],join('_',$d[0],$d[1],$d[2]));
	} elsif ($theMachine eq 'PROTON') {
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
print O Sbwamem($cwd,scalar(keys %fqInfo),$listFQ,$BAMprefix,$FQprefix,$pRef);
close O;
################################
# family.lst, q2mplieup.sh #
################################
my %Families;
open O,'>',$listFamily or die $?;
while ( my $value = $cfam->fetch ) {
	next if $value->{Child} eq '';
	$Families{$value->{Child}} = [
		$value->{Mother},
		$value->{Father},
		$value->{Child},
	];
	print O join("\t",$value->{Mother},$value->{Father},$value->{Child}),"\n";
}
die $cfam->errstr if $cfam->errstr;
#ddx \%Families;
close O;
my $LSTprefix = "$pout/$pPrefixs{lst}";
my $fSHsamt = "$pout/q2mplieup.sh";
open O,'>',$fSHsamt or die $?;
my $VCFprefix = "$pout/$pPrefixs{vcf}";
my $OYKprefix = "$pout/$pPrefixs{oyk}";
print O Smpileup($cwd,scalar(keys %Families),$listFamily,$VCFprefix,$BAMprefix,$LSTprefix,$pRef,$OYKprefix,$theMode);
close O;
################################
for my $iF (keys %Families) {
	my $prefix = "$LSTprefix/p$iF";
	open P,'>',"$prefix.bams.lst" or die $?;
	open M,'>',"$prefix.M.lst" or die $?;
	open F,'>',"$prefix.F.lst" or die $?;
	open C,'>',"$prefix.C.lst" or die $?;
	print M join(' ',@{$Samples{$Families{$iF}->[0]}}),"\n";
	print F join(' ',@{$Samples{$Families{$iF}->[1]}}),"\n";
	print C join(' ',@{$Samples{$Families{$iF}->[2]}}),"\n";
	#for ( @{$Samples{$Families{$iF}->[0]}},@{$Samples{$Families{$iF}->[1]}},@{$Samples{$Families{$iF}->[2]}} ) {
	for (map {@{$Samples{$Families{$iF}->[$_]}}} (0..2) ) {
		my $nbam = "$BAMprefix/$_.bam";
		print P "$nbam\n";
	}
	close P; close M; close F; close C;
}

print "[!]Done !\nNow, check [$listFQ] and [$listFamily] first and run:\n\nqsub ",join('; qsub ',$fSHcutadapt,$fSHbwa,$fSHsamt),"\n\n";
__END__
./genlst.pl chip BGISEQ info.csv fam.csv ref/NIPPT.SNP.5538.fa.gz . ./out/

pip3 install -h cutadapt dnaio xopen
