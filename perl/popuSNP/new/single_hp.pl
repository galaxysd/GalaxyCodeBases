#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.3;

our $opts='i:c:w:l:o:vb';
our ($opt_i, $opt_o, $opt_c, $opt_v, $opt_b, $opt_l, $opt_w);

our $help=<<EOH;
\t-i Population SNP list (./psnp.lst) for *.individual.finalSNPs
\t-c Chr Info file (./chr.nfo) in format [ChrID\\tLen]
\t-w window size (40000) bp
\t-l step length (20000) bp
\t-o Output Prefix (./Hp_result)
\t   Output Files are {./Hp_result}.dat \& {./Hp_result}.stat
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='./Hp_result' if ! defined $opt_o;
$opt_i='./psnp.lst' if ! $opt_i;
$opt_c='./chr.nfo' if ! $opt_c;
$opt_w=40000 if ! $opt_w;
$opt_l=20000 if ! $opt_l;

print STDERR "From [$opt_i] to [$opt_o].{dat,stat}, with [$opt_c][$opt_l][$opt_w]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

#http://doc.bioperl.org/bioperl-live/Bio/Tools/IUPAC.html#BEGIN1
our %IUB = ( A => [qw(A A)],
	     C => [qw(C C)],
	     G => [qw(G G)],
	     T => [qw(T T)],
	     U => [qw(U U)],
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

my $start_time = [gettimeofday];

system('mkdir','-p',$opt_o);
system('rmdir',$opt_o);	# if $opt_o =~ m#/[\s.]*[^/\s.]+[^/]*$#;

my (%ChrLen,@ChrID);
print STDERR "[!]Load ChrNFO:\t";
open C,'<',$opt_c or die "[x]Error opening $opt_c: $!\n";
while (<C>) {
	chomp;
	my ($chr,$len)=split /\t/;
	$ChrLen{$chr}=$len;
	push @ChrID,$chr;
	print STDERR "$chr,$len\t"
}
warn "\n";

my $win=$opt_w-1;

my %SNP;
print STDERR "[!]Parsing PSNP: .";
open P,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";
while (my $file=<P>) {
	chomp $file;
	open SNP,'<',$file or (warn "\n[!]Error opening $file: $!\n" and next);
	print STDERR ".\b";
	my ($chr,$pos,$ref,$tail,$i,@SNPCounts);
	while (<SNP>) {
		($chr,$pos,$ref,$tail)=split /\t/;
		my %SNPcount;
		for (split / /,$tail) {
			next unless /[ACGTRYMKSWHBVDNX]/;
			++$SNPcount{$_} for @{$IUB{$_}};
		}
		@SNPCounts=sort {$a <=> $b} values %SNPcount;
		next if $#SNPCounts < 1;	# at least 2-1=1
		#next if $SNPCounts[0] == $SNPCounts[-1];	# min != max (needed ?)
		$SNP{$chr}{$pos}=[$SNPCounts[0],$SNPCounts[-1]];	# [min,max]
	}
	print STDERR "-";
	close SNP;
}
close P;
warn ".\n";

my $file=$opt_o.'.dat.tmp';
open O,'>',$file or die "[x]Error opening $file: $!\n";

my ($POSa,$N,%Hpc,%Hpa,$Hpr,$sX,$sXX)=(1,0);
print STDERR "[!]Caltulating Hp ";
for my $chr (@ChrID) {
	print STDERR ".\b";
	my ($Cmax,$st,$ed,$pos)=($ChrLen{$chr},1,0);
	my ($Smin,$Smax,$Hp,$sum,$snpcount);
	while ($st <= $Cmax) {
		$ed = $st + $win;	# [$st,$ed]
		$ed = $Cmax if $ed > $Cmax;
		($Smin,$Smax,$Hp,$sum,$snpcount)=(0,0,0,0,0);
		for $pos ($st .. $ed) {
			next unless exists $SNP{$chr}{$pos};
			$Smin += $SNP{$chr}{$pos}->[0];
			$Smax += $SNP{$chr}{$pos}->[1];
			++$snpcount;
		}
		$sum=$Smin+$Smax;
		if ($sum>0) {
			$Hp=2*$Smin*$Smax/($sum*$sum);
			$Hpr=int($Hp*1000)/1000;
		} else {$Hp=$Hpr='Inf';}
		++$Hpc{$chr}{$Hpr};
		++$Hpa{$Hpr};
		++$N;
		$sX += $Hp;
		$sXX += $Hp*$Hp;
		print O "$chr\t$POSa\t$st\t$ed\t$Hp\t$Smin,$Smax,$sum\t$snpcount\n";
		$st += $opt_l;
		$POSa += $opt_l;
	}
	print STDERR '-';
}
close O;
print STDERR "|\b";

my $Avg=$sX/$N;
# http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
#my $Std=sqrt($sXX/$N-$Avg*$Avg);
my $Std=sqrt(($sXX-$Avg*$sX)/($N-1));

$file=$opt_o.'.dat.tmp';
open I,'<',$file or die "[x]Error opening $file: $!\n";
$file=$opt_o.'.dat';
open O,'>',$file or die "[x]Error opening $file: $!\n";
while (<I>) {
	my ($chr, $POSa, $st, $ed, $Hp, $Smms, $snpcount)=split /\t/;
	my $ZHp=($Hp-$Avg)/$Std;
	print O join("\t",$chr, $POSa, $st, $ed, $ZHp, $Hp, $Smms, $snpcount);
}
close O;
warn ";\n";
close I;
unlink $opt_o.'.dat.tmp';

print STDERR "[!]Stating Hp ";
$file=$opt_o.'.stat';
open S,'>',$file or die "[x]Error opening $file: $!\n";
print S "# Avg:\t$Avg\n# Std:\t$Std\n# N: $N, sX: $sX, sXX: $sXX\n";
sub sortInf($$) {
	if ($_[0] eq 'Inf') { return 1; }
	 elsif ($_[1] eq 'Inf') { return -1; }
	 else { return $_[0] <=> $_[1]; }
}

print S "ALL\t$_\t$Hpa{$_}\n" for sort(sortInf keys %Hpa);

for my $chr (@ChrID) {
	print S "$chr\t$_\t$Hpc{$chr}{$_}\n" for sort(sortInf keys %{$Hpc{$chr}});
}
close S;
warn "\n[!] Done !\n";

my $stop_time = [gettimeofday];
print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";

__END__
a=read.delim('v4Hp_40k_20k.dat',F,comment.char='#')

out <- function (c) {
	png(paste("E:\\BGI\\toGit\\popuSNP\\new\\",c,'.png',sep=''),2048,768)
	plot(a$V3[a$V1==c],a$V5[a$V1==c],type='p',xlab=c,ylab='Hp')
	dev.off()
}

out('seg09')

# pdf("E:\\BGI\\toGit\\popuSNP\\new\\wm_a_w40k_s20k.pdf",12,9,fonts='Arial',title="Watermelon Hp density Plot")
plot(density(a$V5),main='Watermelon Hp density Plot',sub='Window:40k bp, Slip:20k bp')
dev.off()

plot(a$V2[a$V5<=0],a$V5[a$V5<=0],type='p',xlab=c,ylab='Hp',cex=.4)
