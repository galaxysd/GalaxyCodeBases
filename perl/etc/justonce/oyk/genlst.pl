#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180516
=cut
use strict;
use warnings;
use Parse::CSV;
use lib '.';
use Data::Dump qw(ddx);
#use FGI::GetCPI;
# ====================================

my $pI = '5';

# ====================================

my %pPrefixs = (
	lst => '0lst',
	fq => '1fq',
	bam => '2bam',
	vcf	=> '3vcf',
	oyki => '4tsv',
);

my @Modes = qw(CHIP PCR);
my %Mode = map { $_ => 1 } @Modes;

my ($theMode,$fninfo,$fnfam,$pchip,$pout) = @ARGV;
$pout = '.' unless $pout;
my $listFQ = "$pout/fq.lst";

my $Usage = "Usage: $0 <mode/help/example> <info.csv> <fam.csv> <chip path> [out path]\n";

my $egInfo = <<'END_EG';
UID,Sample,Cell,Lane,Index
1,HK19M241C,V300016438,L02,20
2,HK19M241F,V300016438,L02,63
3,HK19M241M,V300016438,L02,64
34,NS19M342C,V300016438,L01,55
35,NS19M342F-N,V300016438,L01,19
37,NS19M342M,V300016438,L01,20
END_EG
my $egFam = <<'END_EG';
Father,Mother,Child
HK19M241F,HK19M241M,HK19M241C
NS19M342F-N,NS19M342M,NS19M342C
END_EG
my $takeoutNFO = <<"END_NFO";
#### If the chip path is not in BGISEQ schema, you must ensure Sample names are correct, and then modify [$listFQ].
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
die "[x]Cannot read info.csv [$fninfo].\n" unless -r -s $fninfo;
die "[x]Cannot read fam.csv [$fnfam].\n" unless -r -s $fnfam;
die "[x]Cannot read Chip Path [$pchip].\n" unless -r $pchip;
#die "[x]Cannot use out Path [$pout].\n" unless -r -w -x $pout;
$pchip =~ s/\/+$//g;
warn "[1]Info:[$fninfo], Fam:[$fnfam], CHIP:[$pchip] to Out:[$pout]\n";

my $cinfo = Parse::CSV->new(file => $fninfo, names => 1);
my $cfam = Parse::CSV->new(file => $fnfam, names => 1);
system('mkdir','-p',$pout);
for (keys %pPrefixs) {
	mkdir "$pout/$pPrefixs{$_}";
}

open O,'>',$listFQ or die $?;
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
for (sort keys %fqInfo) {
	my @d = @{$fqInfo{$_}};
	print O join("\t",$_,join('/',$pchip,$d[0],$d[1],join('_',$d[0],$d[1],$pI.$d[2]))),"\n";
}
close O;

my %Families;
while ( my $value = $cfam->fetch ) {
	next if $value->{Child} eq '';
	$Families{$value->{Child}} = [
		$Samples{$value->{Mother}},
		$Samples{$value->{Father}},
		$Samples{$value->{Child}}
	];
}
die $cfam->errstr if $cfam->errstr;
ddx \%Families;
for (keys %Families) {
	my $prefix = "$pout/$pPrefixs{lst}/p$_";
	open P,'>',"$prefix.bams.lst" or die $?;
	open M,'>',"$prefix.M.lst" or die $?;
	open F,'>',"$prefix.F.lst" or die $?;
	open C,'>',"$prefix.C.lst" or die $?;
	;
	close P; close M; close F; close C;
}

__END__
./genlst.pl chip info.csv fam.csv ./chip ./out/

