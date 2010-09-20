#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use File::Basename;
use Galaxy::ShowHelp;
use FindBin qw($RealBin);
#use Data::Dump qw(ddx);

$main::VERSION=0.0.1;
my $SCRIPTS="$RealBin/../scripts";

our $opts='i:o:g:v:bd';
our($opt_i, $opt_o, $opt_v, $opt_g, $opt_b, $opt_d);

#our $desc='';
our $help=<<EOH;
\t-i Raw Genotype Output file (ril.rgt)
\t-g RIL Selfing Generations (10)
\t-o Filtered Genotype Output file (ril.fgt)
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

$opt_i='ril.rgt' if ! $opt_i;
$opt_o='ril.fgt' if ! $opt_o;

no warnings;
$opt_v=int $opt_v;
$opt_g=int $opt_g;
use warnings;
$opt_g=10 if $opt_g < 1;
die "[x]-i $opt_i not exists !\n" unless -f $opt_i;
my $p=0.5 ** $opt_g;
my $vp=(1-$p)/2;

print STDERR "From [$opt_i] with [$vp:$p:$vp] to [$opt_o]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
=pod
`makeEmissionFUN` <- function (errorRate = 0.01) {
	E <- log(errorRate)
	E2 <- log(1 - errorRate)
	E3 <- log(0.5)
	function(h, x, n) {
		if (h != 3)
			return(ifelse(h == x, E2, E))
		else return(n * E3)
	}
}
=cut
sub emissionFUN($$$) {
	my ($h, $x, $n)=@_;
	my $errorRate = 0.01;
	my $E=log($errorRate);
	my $E2=log(1 - $errorRate);
	my $E3=log(0.5);
	if ($h != 2) {	# 0,1,2
		return ($h == $x) ? $E2 : $E;
	} else { return $n * $E3; }
}

=pod
`._rice_phy2get_factor_` <- 244000
`phy2get.haldane.rils` <- function (a, b, l) {
	d <- l/(._rice_phy2get_factor_ * 100)
	p <- (1 - exp(-2 * d))
	p <- p/(1 + p)
	ifelse(a == b, 1 - p, p)
}
=cut
sub transitionFUN() {
	my ($a, $b, $l)=@_;
	my $rice_phy2get_factor=244000;
	my $d=$l/($rice_phy2get_factor * 100);
	my $p=(1 - exp(-2 * $d));
	$p = $p/(1 + $p);
	return ($a==$b) ? (1-$p) : $p;
}

=pod
`hmm.vitFUN.rils` <- function (geno, position, geno.probability, ...) {
	n.obs <- length(geno)
	n.state <- length(geno.probability)
	psi <- delta <- matrix(0, nrow = n.state, ncol = n.obs)
	n.con <- geno.cr <- numeric(n.obs)
	geno.dis <- abs(diff(position))
	n.con[1] <- 1
	g <- geno[1]
	for (i in 2:n.obs) {
		n.con[i] <- ifelse(geno[i] == g, n.con[i - 1] + 1, 1)
		g <- geno[i]
	}
	for (i in 1:n.state) delta[i, 1] <- log(geno.probability[i]) + emissionFUN(i, geno[1], n.con[1])
	preProb <- numeric(n.state)
	for (t in 2:n.obs) {
		for (j in 1:n.state) {
			for (i in 1:n.state) preProb[i] <- delta[i, t - 1] + log(transitionFUN(i, j, geno.dis[t - 1]))
			psi[j, t] <- which.max(preProb)
			delta[j, t] <- max(preProb) + emissionFUN(j, geno[t], n.con[t])
		}
	}
	geno.cr[n.obs] <- which.max(delta[, n.obs])
	for (t in seq(n.obs - 1, 1, by = -1)) geno.cr[t] <- psi[geno.cr[t + 1], t + 1]
	geno.cr
}
=cut
sub absdiff($) {
	my @a=@{$_[0]};
	return undef if scalar @a <= 1;
	my @r;
	my $b=shift @a;
	while (my $i=shift @a) {
		push @r,abs($i-$b);
		$b=$i;
	}
	#print "-[@a]-\n-[@r]-\n";
	return \@r;
}
sub Rmax($) {
	my $arrayr=$_[0];
	my $max=$$arrayr[0];
	my ($i,$which)=(0,0);
	for my $v (@$arrayr) {
		if ($v > $max) {
			$max = $v;
			$which = $i;
		}
		++$i;
	}
	return [$max,$which];
}
sub HMMvitFUNrils($$$) {
	my ($genoref, $position, $probability)=@_;
	my @geno=@$genoref;
	--$_ for @geno;	# Well, @$probability starts from 0
	my $nobs = $#geno;
	my $nstate = $#$probability;
	my (@psi, @delta);	# [$nobs][$nstate]. For delta[, n.obs], defined as [$nobs]->[$nstate], which is reversed to the written order in R
	my (@ncon, @genocr);	# [$nobs]
	my $genodis=absdiff($position);	# [$nobs-1]
#print "-[@$genodis]-\n";
	$ncon[0] = 1;
	my $g = $geno[0];
	for my $i (1..$nobs) {
		$ncon[$i] = ($geno[$i] == $g) ? ($ncon[$i-1]+1) : 1;
		$g = $geno[$i];
	}
	for my $i (0..$nstate) {
		$delta[0][$i] = log($$probability[$i]) + &emissionFUN($i, $geno[0], $ncon[0]);
	}
	my @preProb;	# [$nstate]
	for my $t (1..$nobs) {
		for my $j (0..$nstate) {
			for my $i (0..$nstate) {
				$preProb[$i] = $delta[$t-1][$i] + log(&transitionFUN($i, $j, $$genodis[$t-1]));
			}
			my ($max,$which)=@{&Rmax(\@preProb)};
			$psi[$j][$t] = $which;
			$delta[$t][$j] = $max + &emissionFUN($j, $geno[$t],$ncon[$t]);
		}
	}
	$genocr[$nobs] = &Rmax($delta[-1])->[1];
	for (my $t=$nobs-1;$t>=0;$t--) { $genocr[$t] = $psi[ $genocr[$t+1] ][$t+1]; }
	++$_ for @genocr;
	return \@genocr;
}

my %ReFormatGT=(0=>1, 1=>2, 0.5=>3,'NA'=>0);
my %FormatGT=(1=>0, 2=>1, 3=>0.5, 0=>'NA');
sub filter($$) {
	my @Dat=@{$_[0]};
	my @Pos=@{$_[1]};
	#my $ftype=$_[2] || 0;
	warn 'x' if $#Dat != $#Pos;
	my (@D,@P);
	my $len=$#Dat>$#Pos ? $#Pos : $#Dat;
	for my $i (0..$len) {
		#$Dat[$i]=$ReFormatGT{$Dat[$i]} if $ftype==0;
		next if $Dat[$i] !~ /[123]/;
		push @D,$Dat[$i];
		push @P,$Pos[$i];
	}
	return [\@D,\@P];
}

## main()
open G,'<',$opt_i  or die "[x]Error opening $opt_i: $!\n";
open O,'>',$opt_o  or die "[x]Error opening $opt_o: $!\n";
my $line='#';
my ($Ftype,@Values,%Count)=(0);	# $Ftype -> 0 => 0,1,0.5,NA; 1 => 1,2,3,0
while ($line =~ /^#/) {
	$line=<G>;
	push @Values,$1 if $line =~ /^[AB]=(\d+)\D/;
	print O $line;
}
if (@Values) {
	@Values=sort {$a <=> $b} @Values;
	$Ftype=1 if $Values[1]==2;
}
chomp $line;
#$line=(split "\t",$line)[1];
my @Samples=split /\t/,$line;
print STDERR '[!]Inputed [',scalar @Samples,'] Samples, ';
my (@Poses,%Dat,$i);
while (<G>) {
	next if /^#/;
	chomp;
	my ($Pos,@DatLine)=split "\t";
	push @Poses,$Pos;
	#push @{$Dat{$_}},shift @DatLine for @Samples;
	for my $Sample (@Samples) {
		my $t=shift @DatLine;
		$t=$ReFormatGT{$t} if $Ftype==0;
		push @{$Dat{$Sample}},$t;
		++$Count{$t};
	}
}
print STDERR '[',scalar @Poses,"] Poses.\n";
close G;
#ddx \%Dat;
@Values=();
push @Values,"$_:$Count{$_}" for sort {$a<=>$b} keys %Count;
$i=join ', ',@Values;
warn "[!]In: $i\n";
print O "#In: $i\n";

my (%DatHmm,%SampleH,%Diff);
%Count=();
for my $Sample (@Samples) {
	my ($D,$P)=@{&filter($Dat{$Sample},\@Poses)};
	$line = &HMMvitFUNrils($D,$P,[$vp, $vp, $p]);
	%SampleH=();
	@SampleH{@$P} = @$line;
	my ($diff,@v)=0;
	$i=0;
	for (@Poses) {
		if (exists $SampleH{$_}) {
			push @v,$SampleH{$_};
			++$Count{$SampleH{$_}};
			++$diff if $SampleH{$_} ne ${$Dat{$Sample}}[$i];
			++$i;
		} else { push @v,0;++$Count{0}; }
	}
	$DatHmm{$Sample}=\@v;
	$Diff{$Sample}=[$diff,$i];
}
#ddx \%Diff;
@Values=();
push @Values,"$_:$Count{$_}" for sort {$a<=>$b} keys %Count;
$i=join ', ',@Values;
warn "[!]Out: $i\n";
print O "#Out: $i\n";

for my $Pos (@Poses) {
	print O $Pos;
	for my $Sample (@Samples) {
		my $t=shift @{$DatHmm{$Sample}};
		$t=$FormatGT{$t};
		print O "\t",$t;
	}
	print O "\n";
}
print O '#Filtered of ',scalar @Samples,' x ',scalar @Poses,":\n";
my $SumDif;
$i=0;
for my $Sample (@Samples) {
	my ($diff,$ii)=@{$Diff{$Sample}};
	$SumDif+=$diff; $i+=$ii;
	print O "# $Sample\t$diff/$ii = ",$diff/$ii,"\n";
}
print O "# Average:\t$SumDif/$i = ",$SumDif/$i,"\n";
close O;
warn "[!]Filtered $SumDif/$i = ",$SumDif/$i,"\n";

=pod
my @O = (1,1,1,1,1,1,1,2,1,1,1,1,1,1,2,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,1,2,2,2,2,
     2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,1,1,2,2,1,2,1,2,1,2,1,2,1,2,1,1,1,1,1,1,1,2,1,1,1,1,
     1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,2,1,1,1,1,1);
my (@Opos,$out);
my $t=30e3;
for (@O) {
	push @Opos,$t;
	$t += 30e3;
};
$out=&HMMvitFUNrils(\@O,\@Opos,[0.4975, 0.4975,0.005]);
my $shouldbe='1111111111111111111111222222222222222222222222222222222222223333333333333333111111111111111111111111111111111';
my $got=join('',@$out);
print "Passed !\n" if $shouldbe eq $got;
print 'I [',join('',@O),"]\n";
print 'O [',$got,"]\n";
=cut

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
