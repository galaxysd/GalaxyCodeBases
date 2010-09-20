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

our $opts='i:o:m:z:c:v:g:q:e:f:n:bd';
our($opt_i, $opt_o, $opt_m, $opt_z, $opt_c, $opt_v, $opt_g, $opt_q, $opt_e, $opt_f, $opt_n, $opt_b, $opt_d);

#our $desc='';
our $help=<<EOH;
\t-i RIL pSNP list (ril.lst) in format: /^ChrID\\tpath toadd_ref\$/
\t-m Parent A list (a.lst) in format: /^ChrID\\tpath to CNS\$/
\t-z Parent B list (b.lst)
\t-c ChrID to parse
\t-o Raw Genotype Output prefix (ril).ChrID.rgt
\t-g Dump Ref. Genotype to (ref).ChrID.gt
\t-q MinQual (15)
\t-n MaxCopyNum (1)
\t-e MinDepth (2)
\t-f MaxDepth (100)
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

$opt_i='ril.lst' if ! $opt_i;
$opt_m='a.lst' if ! $opt_m;
$opt_z='b.lst' if ! $opt_z;
$opt_g='ref' if ! $opt_g;
$opt_o='ril' if ! $opt_o;
$opt_q=15 if ! $opt_q;
$opt_e=2 if ! $opt_e;
$opt_f=100 if ! $opt_f;
$opt_n=1 if ! $opt_n;
die "[x]Must specify -c ChrID !\n" unless defined $opt_c;

no warnings;
$opt_v=int $opt_v;
$opt_q=int $opt_q;
$opt_e=int $opt_e;
$opt_f=int $opt_f;
use warnings;
die "[x]-i $opt_i not exists !\n" unless -f $opt_i;
die "[x]-m $opt_m not exists !\n" unless -f $opt_m;
die "[x]-z $opt_z not exists !\n" unless -f $opt_z;
my ($fileM,$fileZ,$fileRIL);

print STDERR "From [$opt_i][$opt_c] with [$opt_m][$opt_z],\n     [minQ:$opt_q,Depth:$opt_e,$opt_f,maxCopyNum:$opt_n]\n To  [$opt_g].$opt_c.gt, [$opt_o].$opt_c.rgt\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
# A1 T2 C4 G8
our %bIUB = ( A => 1,
		C => 4,
		G => 8,
		T => 2,
		M => 5,#[qw(A C)],
		R => 9,#[qw(A G)],
		W => 3,#[qw(A T)],
		S => 12,#[qw(C G)],
		Y => 6,#[qw(C T)],
		K => 10,#[qw(G T)],
		V => 13,#[qw(A C G)],
		H => 7,#[qw(A C T)],
		D => 11,#[qw(A G T)],
		B => 14,#[qw(C G T)],
		X => 15,#[qw(G A T C)],
		N => 15,#[qw(G A T C)]
		);
our %REV_IUB = (1	=> 'A',
		2	=> 'T',
		4	=> 'C',
		8 	=> 'G',
		5	=> 'M,AC',
		9	=> 'R,AG',
		3	=> 'W,AT',
		12	=> 'S,CG',
		6	=> 'Y,CT',
		10	=> 'K,GT',
		13	=> 'V,ACG',
		7	=> 'H,ACT',
		11	=> 'D,AGT',
		14	=> 'B,CGT',
		15	=> 'N',
		);
=pod
#http://doc.bioperl.org/bioperl-live/Bio/Tools/IUPAC.html#BEGIN1
our %IUB = ( A => [qw(A)],
	     C => [qw(C)],
	     G => [qw(G)],
	     T => [qw(T)],
	     U => [qw(U)],
	     M => [qw(A C)],
	     R => [qw(A G)],
	     W => [qw(A T)],
	     S => [qw(C G)],
	     Y => [qw(C T)],
	     K => [qw(G T)],
	     V => [qw(A C G)],
	     H => [qw(A C T)],
	     D => [qw(A G T)],
	     B => [qw(C G T)],
	     X => [qw(G A T C)],
	     N => [qw(G A T C)]
	     );
my %x = ( A => 1,
	     C => 4,
	     G => 8,
	     T => 2,);
for my $k (sort keys %IUB) {
	next if $k eq 'U';
	my $r=0;
	$r += $x{$_} for (@{$IUB{$k}});
	print join("\t",$k,$r,$bIUB{$k}),"\n";
	warn "!" if $r != $bIUB{$k};
}
=cut

=pod
SoapSNP格式：
参见：http://soap.genomics.org.cn/soapsnp.html#output2
chromosome_10	1226	A	M	36	C	28	4	5	A	33	3	6	11	0.0285714	1.54545	0	166
1)	Chromosome name
2)	Position of locus
3)	Nucleotide at corresponding locus of reference sequence
4)	Genotype of sequencing sample
5)	Quality value
6)	nucleotide with the highest probability(first nucleotide)
7)	Quality value of the nucleotide with the highest probability
8)	Number of supported reads that can only be aligned to this locus
9)	Number of all supported reads that can be aligned to this locus
10)	Nucleotide with higher probability
11)	Quality value of nucleotide with higher probability
12)	Number of supported reads that can only be aligned to this locus
13)	Number of all supported reads that can be aligned to this locus
14)	Total number of reads that can be aligned to this locus
15)	RST p-value
16)	Estimated copy number for this locus
17)	Presence of this locus in the dbSNP database. 1 refers to presence and 0 refers to inexistence
18)	The distance between this locus and another closest SNP
=cut
open M,'<',$opt_m  or die "[x]Error opening $opt_m: $!\n";
while (<M>) {
	chomp;
	my ($ChrID,$file)=split /\t/;
	next if $ChrID ne $opt_c;
	$fileM=$file and last;
}
close M;
open RIL,'<',$opt_i  or die "[x]Error opening $opt_i: $!\n";
while (<RIL>) {
	chomp;
	my ($ChrID,$file)=split /\t/;
	next if $ChrID ne $opt_c;
	$fileRIL=$fileZ=$file and last;
}
close RIL;
open M,'<',$opt_z  or die "[x]Error opening $opt_z: $!\n";
while (<M>) {
	chomp;
	my ($ChrID,$file)=split /\t/;
	next if $ChrID ne $opt_c;
	$fileZ=$file and last;
}
close M;

warn "[!]Using [$fileRIL]<-[$fileM][$fileZ].\n";

my ($ChrLen,$InZonea,$InZoneb,%SNPRefContent,%HeteroHomoContent,%DatS,$t,$NumType,$NumTypez,@OutTmp)=(0,0,0);
%SNPRefContent=%HeteroHomoContent=(1=>0, 2=>0, 3=>0);
open G,'>',$opt_g.".$opt_c.gt" or die "[x]Error opening $opt_g.$opt_c.gt: $!\n";
print G "# A1 T2 C4 G8
# minQ:$opt_q, Depth:$opt_e,$opt_f, maxCopyNum:$opt_n
# ChrID=$opt_c
# A=$fileM
# B=$fileZ
",join("\t",'Pos','A','B','a','b'),"\n";
open M,'<',$fileM  or die "[x]Error opening $fileM: $!\n";
open Z,'<',$fileZ  or die "[x]Error opening $fileZ: $!\n";
while (<M>) {
	#chomp;
	my ($ChrID,$Pos,$Ref,$Type,$Q,$Dep1,$Dep2,$DepA,$CN)=(split /\t/)[0,1,2,3,4,8,12,13,15];
	$t=<Z>;
	my ($ChrIDz,$Posz,$Refz,$Typez,$Qz,$Dep1z,$Dep2z,$DepAz,$CNz)=(split /\t/,$t)[0,1,2,3,4,8,12,13,15];
	++$ChrLen;
	print STDERR "=>\b" unless $ChrLen % 1_000_000;
	next if ($DepA < $opt_e or $DepA > $opt_f) or ($DepAz < $opt_e or $DepAz > $opt_f);
	die "[x]CNS Files NOT match !\n" if $ChrID ne $ChrIDz or $Pos ne $Posz or $Ref ne $Refz;
	my ($flag,$het)=(3,0);
	if ( $Q < $opt_q or $CN > $opt_n or ($Type !~ /[ATCG]/i and ($het |= 1) and ($Dep1<$opt_e or $Dep2<$opt_e)) ) {
		$Type = $Ref;
		$flag &= 2;
	} else { ++$InZonea; }
	if ( $Qz < $opt_q or $CNz > $opt_n or ($Typez !~ /[ATCG]/i and ($het |= 2) and ($Dep1z<$opt_e or $Dep2z<$opt_e)) ) {
		$Typez = $Refz;
		$flag &= 1;
	} else { ++$InZoneb; }
	next if $Type eq $Typez;
	$NumType=$bIUB{$Type}; $NumTypez=$bIUB{$Typez};
	next if $NumType & $NumTypez;	# keep only A ^ B = NULL
	$DatS{$Pos}=[$NumType,$NumTypez,15^($NumType | $NumTypez)];
	print G join("\t",$Pos,$Type,$Typez,$NumType,$NumTypez),"\n";
	++$SNPRefContent{$flag};
	++$HeteroHomoContent{$het};
}
@OutTmp=();
$t=0;
for (sort keys %SNPRefContent) {
	push @OutTmp,"$_:$SNPRefContent{$_}";
	$t += $SNPRefContent{$_};
}
print G "# ChrLen:$ChrLen\n# InZone: $InZonea, $InZoneb\n# SNP-Ref: ",join(', ',@OutTmp,"All:$t"),' -> ~',int(.5+$ChrLen/$t),' bp';
print STDERR "|\n[!] ChrLen:$ChrLen\n[!] InZone: $InZonea, $InZoneb\n[!] SNP-Ref: ",join(', ',@OutTmp,"All:$t"),' -> ~',int(.5+$ChrLen/$t),' bp';
@OutTmp=();
my $t2=0;
for (sort keys %HeteroHomoContent) {
	push @OutTmp,"$_:$HeteroHomoContent{$_}";
	$t2 += $HeteroHomoContent{$_} if $_;
}
print G "\n# Hetero-Homo: ",join(', ',@OutTmp,"All_but_0:$t2"),' -> ',100*$t2/$t," %\n";
print STDERR "\n[!] Hetero-Homo: ",join(', ',@OutTmp,"All_but_0:$t2"),' -> ',100*$t2/$t," %\n";
close Z;
close M;
close G;
$t2=join('','-' x 24,'+','-' x 24,'|','-' x 24,'+---',"\n");
warn "[!]Parents' SNP done !\n",$t2;

sub GetGT($$) {
	my ($Pos,$binBase)=@_;
	my $Result=0;
	#return 0 unless exists $DatS{$Pos};
	my ($a,$b,$c)=@{$DatS{$Pos}} or return 0;
	$Result |=1 if $binBase & $a;	# A
	$Result |=2 if $binBase & $b;	# B
	$Result |=4 if $binBase & $c;	# Extra GenoType
	return $Result;
}	# %DatS

my ($C_PSNP,$C_iSNP,$C_Pos,%GT)=(0,0,0);
open IN,'<',$fileRIL or die "[x]Error opening $fileRIL: $!\n";
chomp($t=<IN>);
$t=(split /\t/,$t)[3];
my @Samples=split / /,$t;
warn '[!]Sample Order: ',(scalar @Samples),':[',join('] [',@Samples),']',"\n";
my ($chr,$pos,$ref,$tail);
while (<IN>) {
	chomp;
	($chr,$pos,$ref,$tail)=split /\t/;
	next if $chr ne $opt_c;
	++$C_PSNP;
	next unless $DatS{$pos};
	++$C_Pos;
	$t=0;
	my @indSNP=split / /,$tail;	# /[ACGTRYMKSWHBVDNX-]/
#die join "\t",$chr,$pos,$ref,$tail,scalar @indSNP if @indSNP < 135;
	for my $s (@Samples) {
		unless ($tail=shift @indSNP) {
			warn "[!].add_ref Error at ($chr,$pos: $ref) !\n" unless $t;
			$t=1;
			next;
		}
		if ($tail ne '-') {
			$GT{$pos}{$s}=&GetGT($pos,$bIUB{$tail});	# if $DatBoth{$pos};
			++$C_iSNP;
		} else {
			next;
		}
	}
}
close IN;
warn "[!]GenoTyping done as ${C_PSNP}->$C_Pos,$C_iSNP,",$C_iSNP/$C_Pos,"\n";

my $Deleted=0;
unless ($opt_d) {
	for my $Pos (keys %GT) {
		my $flag=0;
		for my $s (keys %{$GT{$Pos}}) {
			if ($GT{$Pos}{$s} & 4) {
				delete $GT{$Pos}{$s};
				++$Deleted;
			}
		}
	}
	warn "[!]GenoType filtered out $Deleted\n";
}

my ($PosC,$TypeC,%C_GT)=(0,0);
open O,'>',$opt_o.".$opt_c.rgt" or die "[x]Error opening $opt_o: $!\n";
my @PosList=sort {$a<=>$b} keys %GT;
$PosC = scalar @PosList;

my %FormatGT=(1=>0, 2=>1, 3=>0.5, 0=>'NA');	# 0 for M(1), 1 for Z(2), 0.5 for mixture(3)
$fileM=`readlink -nf $fileM`;
$fileZ=`readlink -nf $fileZ`;
$fileRIL=`readlink -nf $fileRIL`;
print O "#A=0, $fileM\n#B=1, $fileZ\n#RIL $fileRIL\n",join("\t",@Samples),"\n";
for my $pos (@PosList) {
	print O $pos;
	for my $s (@Samples) {
		if (exists $GT{$pos}{$s}) {
			$t=$FormatGT{$GT{$pos}{$s}};
			print O "\t",$t;
			++$TypeC;
			++$C_GT{$t};
		} else {
			print O "\t",'NA';
			++$C_GT{'NA'};
		}
	}
	print O "\n";
}
=pod
print O 'PosList',"\t",join(',',@PosList),"\n";
for my $s (@Samples) {
	print O $s;
	for my $pos (@PosList) {
		if (exists $GT{$pos}{$s}) {
			print O "\t",$GT{$pos}{$s};
			++$TypeC;
			++$C_GT{$GT{$pos}{$s}};
		} else {
			print O "\t",'-';
			++$C_GT{'-'};
		}
		#print O "\t",exists($GT{$pos}{$s})?$GT{$pos}{$s}:'-';	# exists is faster than defined ?
	}
	print O "\n";
}
=cut
warn "[!]GenoType written $PosC,$TypeC,",$TypeC/$PosC,"\n\n[!]GenoType Counts:\n";
for (sort keys %C_GT) {
	warn "\t$_: $C_GT{$_}\n";
	print O "# $_: $C_GT{$_}\n";
}
print O "#Filtered out: $Deleted\n#GenoType written $PosC,$TypeC -> ",$TypeC/$PosC,"\n";
close O;

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
cat ./9311/chrorder |xargs -n1 ./genotyping.pl -bz lpa64cns.lst -m l9311cns.lst -o ./20100916/ril -g ./20100916/ref -c
