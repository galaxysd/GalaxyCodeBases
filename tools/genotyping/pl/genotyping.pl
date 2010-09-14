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

our $opts='i:o:m:z:c:v:g:bd';
our($opt_i, $opt_o, $opt_m, $opt_z, $opt_c, $opt_v, $opt_g, $opt_b, $opt_d);

#our $desc='';
our $help=<<EOH;
\t-i RIL pSNP list (ril.lst) in format: /^ChrID\\tpath toadd_ref\$/
\t-m Parent M list (m.lst) in format: /^ChrID\\tpath to SNP\$/
\t-z Parent Z list (undef), undef for using the Reference sequence
\t-c ChrID to parse
\t-o Raw Genotype Output prefix (ril).ChrID.rgt
\t-g Dump Ref. Genotype to (undef)
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

$opt_i='ril.lst' if ! $opt_i;
$opt_m='m.lst' if ! $opt_m;
$opt_o='ril' if ! $opt_o;
die "[x]Must specify -c ChrID !\n" unless defined $opt_c;

no warnings;
$opt_v=int $opt_v;
use warnings;
die "[x]-i $opt_i not exists !\n" unless -f $opt_i;
die "[x]-m $opt_m not exists !\n" unless -f $opt_m;
my ($fileM,$fileZ,$fileRIL);
$fileZ=$opt_z?$opt_z:'_Ref_';

print STDERR "From [$opt_i][$opt_c] with [$opt_m],[$fileZ] to [$opt_o].$opt_c.rgt\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Ref GenoType Dump to [$opt_g]\n" if $opt_g;
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
chromosome_10	1226	A	M	36	C	28	4	5	A	33	3	6	11	0.0285714	1.54545	0	166
1)	Chromosome name
2)	Position of locus
3)	Nucleotide at corresponding locus of reference sequence
4)	Genotype of sequencing sample
5)	Quality value
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
if (defined $opt_z) {
	open M,'<',$opt_z  or die "[x]Error opening $opt_z: $!\n";
	while (<M>) {
		chomp;
		my ($ChrID,$file)=split /\t/;
		next if $ChrID ne $opt_c;
		$fileZ=$file and last;
	}
	close M;
}

warn "[!]Using [$fileRIL]->[$fileM][$fileZ].\n";

my (%DatM,%DatZ,%DatRef,%DatBoth);
open M,'<',$fileM  or die "[x]Error opening $fileM: $!\n";
while (<M>) {
	#chomp;
	my ($ChrID,$Pos,$Ref,$Type,$Q)=split /\t/;
	$DatM{$Pos}=$DatBoth{$Pos}=$bIUB{$Type};
	$DatRef{$Pos}=$bIUB{$Ref};
}
close M;
if (defined $opt_z) {
	open Z,'<',$fileZ  or die "[x]Error opening $fileZ: $!\n";
	while (<Z>) {
		#chomp;
		my ($ChrID,$Pos,$Ref,$Type,$Q)=split /\t/;
		$DatZ{$Pos}=$bIUB{$Type};
		#$DatBoth{$Pos} += 0;
		$DatBoth{$Pos} |= $bIUB{$Type};
		$DatRef{$Pos}=$bIUB{$Ref};
	}
	close Z;
} else { %DatZ = %DatRef; }
=pod
# & AND, | OR, ^ XOR
# $M = $DatM{$Pos} ^ ($DatM{$Pos} & $DatRef{$Pos});	# also $DatM{$Pos} & ( ~$DatRef{$Pos} )
my (%DatC,$Common,$M,$Z);
for my $Pos (keys %DatRef) {
	if (exists $DatM{$Pos}) {
		if (exists $DatZ{$Pos}) {
			# Both exists.
			$Common = $DatM{$Pos} & $DatZ{$Pos};
			$M = $DatM{$Pos} ^ $Common;
			$Z = $DatZ{$Pos} ^ $Common;
		} else {
			# Only $DatM{$Pos}
			$M = $DatM{$Pos} ^ ($DatM{$Pos} & $DatRef{$Pos});	# also $DatM{$Pos} & ( ~$DatRef{$Pos} )
			$Z = $DatRef{$Pos};
		}
	} else {
		# Only $DatZ{$Pos}, since (keys %DatRef) = (keys %DaM) U (keys %DatZ)
		$Z = $DatZ{$Pos} ^ ($DatZ{$Pos} & $DatRef{$Pos});
		$M = $DatRef{$Pos};
	}
	$DatGT{$Pos} = [$M,$Z];
}
(%DatM,%DatZ,%DatRef)=();
=cut
%DatRef=();
my $C_SameGT=0;
for my $Pos (keys %DatBoth) {
	if ($DatM{$Pos} && $DatZ{$Pos} && $DatM{$Pos} == $DatZ{$Pos}) {
		delete $DatM{$Pos};
		delete $DatZ{$Pos};
		delete $DatBoth{$Pos};
		++$C_SameGT;
	} else { $DatBoth{$Pos} = ~$DatBoth{$Pos}; }
}
#$DatBoth{$_} = ~$DatBoth{$_} for keys %DatBoth;
warn "[!]pSNP loaded. [$C_SameGT] positions are identical.\n";

if ($opt_g) {
	open D,'>',$opt_g or die "[x]Error opening $opt_g: $!\n";
	print D join("\t",'Pos','Ref','M','Z','Both'),"\n";
	for my $Pos (sort {$a<=>$b} keys %DatRef) {
		print D join("\t",$Pos,$REV_IUB{$DatRef{$Pos}},$DatM{$Pos}?($REV_IUB{$DatM{$Pos}}):'-',$DatZ{$Pos}?($REV_IUB{$DatZ{$Pos}}):'-',$REV_IUB{~$DatBoth{$Pos}}),"\n";
		#print D join("\t",$Pos,$DatRef{$Pos},$REV_IUB{$DatM{$Pos}},$DatZ{$Pos},$DatBoth{$Pos}),"\n";
	}
	close D;
}

sub GetGT($$) {
	my ($Pos,$binBase)=@_;
	my $Result=0;
	$Result |=1 if $DatM{$Pos} && ($binBase & $DatM{$Pos});	# M
	$Result |=2 if $DatZ{$Pos} && ($binBase & $DatZ{$Pos});	# Z
	$Result |=4 if $binBase & $DatBoth{$Pos};	# Extra GenoType
	return $Result;
}

my ($C_PSNP,$C_iSNP,$C_Pos,%GT)=(0,0,0);
open IN,'<',$fileRIL or die "[x]Error opening $fileRIL: $!\n";
my $line;
chomp($line=<IN>);
$line=(split /\t/,$line)[3];
my @Samples=split / /,$line;
warn '[!]Sample Order: ',(scalar @Samples),':[',join('] [',@Samples),']',"\n";
my ($chr,$pos,$ref,$tail,$t);
while (<IN>) {
	chomp;
	($chr,$pos,$ref,$tail)=split /\t/;
	next if $chr ne $opt_c;
	++$C_PSNP;
	next unless $DatBoth{$pos};
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

unless ($opt_d) {
	my $Deleted=0;
	for my $Pos (keys %GT) {
		my $flag=0;
		for my $s (keys %{$GT{$Pos}}) {
			delete $GT{$Pos}{$s} if $GT{$Pos}{$s} & 4;
			++$Deleted;
		}
	}
	warn "[!]GenoType filtered out $Deleted\n";
}

my ($PosC,$TypeC,%C_GT)=(0,0);
open O,'>',$opt_o.".$opt_c.rgt" or die "[x]Error opening $opt_o: $!\n";
my @PosList=sort {$a<=>$b} keys %GT;
$PosC = scalar @PosList;

my %FormatGT=(1=>0, 2=>1, 3=>0.5);	# 0 for M(1), 1 for Z(2), 0.5 for mixture(3)
$fileM=`readlink -nf $fileM`;
$fileZ=`readlink -nf $fileZ`;
$fileRIL=`readlink -nf $fileRIL`;
print O "#M=0, $fileM\n#Z=1, $fileZ\n#Chimeric_SNP=$C_SameGT\n#RIL $fileRIL\n",join("\t",@Samples),"\n";
for my $pos (@PosList) {
	print O $pos;
	for my $s (@Samples) {
		if (exists $GT{$pos}{$s}) {
			$line=$FormatGT{$GT{$pos}{$s}};
			print O "\t",$line;
			++$TypeC;
			++$C_GT{$line};
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
close O;

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
cat ./9311/chrorder |xargs -n1 ./genotyping.pl -bz z.lst -c
