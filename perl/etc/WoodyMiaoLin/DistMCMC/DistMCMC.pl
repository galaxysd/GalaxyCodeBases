#!/usr/bin/perl
use strict;
use warnings;

die("
Program: DistMCMC
DistMCMC does maximum likelihood extimation of pairwise distance using Markov chain Monte Carlo algorithm.
Only JC69, JC69+G, K80 and K80+G are available now.\n
Version: 1.0
Release: Aug. 3, 2014
Auther: Woody\n
Usage: $0 <in.fasta> <out.txt>
Please edit the Perl script directly if you want to make changes of parameters.
\n") if (@ARGV<2);

my $NoItera = 100000; # Number of iterations
my $Burnin = 10000; # Burn in
my $Savfreq = 100; # Frequency of saving states after burn in
my $Prifreq = 50000; # Frequency of printing states on screen
my $dSlidWin = 0.1; # Sliding window size of distance
my $kSlidWin = 10; # Sliding window size of kappa
my @dRange = (0, 1); # Range of distance
my @kRange = (0, 100); # Range of kappa
my $alpha1 = 1; # Gamma shape parameter
my $alpha2 = 0.6;
my $alpha3 = 0.2;

my $i = shift;
my $o = shift;
open my $i1, "<", "$i";
open O, ">", "$o";

my $nchar;
my %dna = %{&readfasta($i1)};
close $i1;

my %nseg; # Hash of pairwise number of segregation sites
foreach my $s1 (sort keys %dna) {
	foreach my $s2 (sort keys %dna) {
		if ($nseg{"$s1 $s2"} or $nseg{"$s2 $s1"}) {
			next;
		} elsif ($s1 eq $s2) {
			next;
		} else {
			$nseg{"$s1 $s2"} = &calnseg($dna{$s1}, $dna{$s2});
		}
	}
}

print O "#Pair\tDistance\tDistance_mean\tDistance_median\tDistance_%95_HPD_lower\tDistance_%95_HPD_upper\tKappa_mean\tKappa_median\tKappa_%95_HPD_lower\tKappa_%95_HPD_upper\n";

foreach (sort keys %nseg) {
	print "\n$_";
	print "\nJC69 MCMC simulating...\nNo. Iteration\tDistance\n";
	my @a = &JC69MCMC(@{$nseg{$_}});
	print O $_, "\tJC69\t", join("\t", @a), "\n";
	print "\nJC69+G MCMC simulating...\nNo. Iteration\tDistance\n";
	@a = &JC69gammaMCMC(@{$nseg{$_}}, $alpha1);
	print O $_, "\tJC69+G a=$alpha1\t", join("\t", @a), "\n";
	@a = &JC69gammaMCMC(@{$nseg{$_}}, $alpha2);
	print O $_, "\tJC69+G a=$alpha2\t", join("\t", @a), "\n";
	@a = &JC69gammaMCMC(@{$nseg{$_}}, $alpha3);
	print O $_, "\tJC69+G a=$alpha3\t", join("\t", @a), "\n";
	print "\nK80 MCMC simulating...\nNo. Iteration\tDistance\tKappa\n";
	@a = &K80MCMC(@{$nseg{$_}});
	print O $_, "\tK80\t", join("\t", @a), "\n";
	print "\nK80+G MCMC simulating...\nNo. Iteration\tDistance\tKappa\n";
	@a = &K80gammaMCMC(@{$nseg{$_}}, $alpha1);
	print O $_, "\tK80+G a=$alpha1\t", join("\t", @a), "\n";
	@a = &K80gammaMCMC(@{$nseg{$_}}, $alpha2);
	print O $_, "\tK80+G a=$alpha2\t", join("\t", @a), "\n";
	@a = &K80gammaMCMC(@{$nseg{$_}}, $alpha3);
	print O $_, "\tK80+G a=$alpha3\t", join("\t", @a), "\n";
}
close O;
	
sub readfasta { # Read fasta file
	my $in = $_[0];
	my %fa;
	while (<$in>) {
		$/ = ">";
		my $seq = <$in>;
		chomp $seq;
		$seq =~ s/\s//g;
		$/ = "\n";
		$nchar = length $seq unless $nchar;
		s/>//;
		/^(\S+)\s/;
		$fa{$1} = $seq;
	}
	return \%fa;
}

sub calnseg {
	my ($s1, $s2) = @_;
	my $total;
	my $ti = 0;
	my $tv = 0;
	foreach my $i (0 .. $nchar-1) {
	       next unless (substr($s1, $i, 1) =~ /[ATCG]/) and (substr($s2, $i, 1) =~ /[ATCG]/);
	       ++$total;
	       if (substr($s1, $i, 1) eq "A") {
		       if (substr($s2, $i, 1) eq "G") {
			       ++$ti;
		       } elsif (substr($s2, $i, 1) eq "A") {
			       next;
		       } else {
			       ++$tv;
		       }
	       } elsif (substr($s1, $i, 1) eq "G") {
		       if (substr($s2, $i, 1) eq "G") {
			       next;
		       } elsif (substr($s2, $i, 1) eq "A") {
			       ++$ti;
		       } else {
			       ++$tv;
		       }
	       } elsif (substr($s1, $i, 1) eq "C") {
		       if (substr($s2, $i, 1) eq "T") {
			       ++$ti;
		       } elsif (substr($s2, $i, 1) eq "C") {
			       next;
		       } else {
			       ++$tv;
		       }
	       } elsif (substr($s1, $i, 1) eq "T") {
		       if (substr($s2, $i, 1) eq "T") {
			       next;
		       } elsif (substr($s2, $i, 1) eq "C") {
			       ++$ti;
		       } else {
			       ++$tv;
		       }
	       }
       }
       return [$ti, $tv, $total];
}


sub NewSta { # Propose new state
	my ($OldStat, $SlidWin, $LoBound, $UpBound) = @_;
	my $r1 = rand;
	my $TemStat = $OldStat + ($r1-0.5)*$SlidWin; # Temp state
	my $NewStat; # New state
	if ($TemStat < $LoBound) {
		$NewStat = 2*$LoBound - $TemStat;
	} elsif ($TemStat > $UpBound) {
		$NewStat = 2*$UpBound - $TemStat;
	} else {
		$NewStat = $TemStat;
	}
	return $NewStat;
}

sub statis { # Calculate mean, median and 95% HPD of numbers
	my @d_sort = sort {$a <=> $b} @_;
	my $num_d = @_; # Number of samples in @_
	my $num_mid = int($num_d/2); # Middle number of @d_sort;
	my $d_mid = $d_sort[$num_mid];
	my $sum;
	foreach (@_) {
		$sum += $_;
	}
	my $d_mean = $sum/@_;
	my $num95_d = 0.95*$num_d; # %95 HPD number of samples
	my $intHPD95 = $d_sort[$num95_d-1] - $d_sort[0]; # %95 HPD interval
	my $intHPD95_l = $d_sort[0]; # Left boundary of %95 HPD interval
	my $intHPD95_r = $d_sort[$num95_d-1]; # Right boundary of %95 HPD interval
	for (my $i = 1; $num95_d+$i <= $num_d; ++$i) {
		my $int95 = $d_sort[$num95_d+$i-1] - $d_sort[$i]; # %95 interval
		if ($int95 < $intHPD95) {
			$intHPD95 = $int95;
			$intHPD95_l = $d_sort[$i];
			$intHPD95_r = $d_sort[$num95_d+$i-1];
		}
	}
	return ($d_mean, $d_mid, $intHPD95_l, $intHPD95_r);
}

sub JC69ML {
	my ($ti, $tv, $n) = @_;
	my $p = ($ti+$tv)/$n;
	my $d_mean = -0.75*log(1-4/3*$p);
	my $d_stde = sqrt($p*(1-$p)/(1-4/3*$p)**2/$n);
	my $ci_l = $d_mean - 1.96*$d_stde;
	my $ci_r = $d_mean + 1.96*$d_stde;
	return ($d_mean, $d_stde, $ci_l, $ci_r);
}

sub JC69MCMC {
	my ($ti, $tv, $n) = @_;
	my $d = ($dRange[1]-$dRange[0])*(rand) + $dRange[0]; # Random initial state
	my @d; # An array of saved states
	foreach (1 .. $NoItera) {
		my $d_old = $d; # Old state
		my $d_new = &NewSta($d_old, $dSlidWin, @dRange);
		my $logOld = ($ti+$tv)*log(0.75-0.75*exp(-4/3*$d_old)) + ($n-$ti-$tv)*log(0.25+0.75*exp(-4/3*$d_old)); # Log likelihood of old state
		my $logNew = ($ti+$tv)*log(0.75-0.75*exp(-4/3*$d_new)) + ($n-$ti-$tv)*log(0.25+0.75*exp(-4/3*$d_new)); # Log likelihood of new state
		my $accRat; # Acceptance ratio
		if ($logNew > $logOld) {
			$accRat = 1;
		} else {
			$accRat = exp($logNew-$logOld);
		}
		my $r2 = rand;
		if ($r2 < $accRat) {
			$d = $d_new;
		}
		my $pri = $_/$Prifreq;
		print "$_\t$d\n" if $pri == int($pri);
		if ($_ > $Burnin) {
			my $sav = $_/$Savfreq;
			push @d, $d if $sav == int($sav);
		}
	}
	my @d_stat = &statis(@d);
	return (@d_stat);
}

sub K80ML {
	my ($ti, $tv, $n) = @_;
	my $S = $ti/$n;
	my $V = $tv/$n;
	my $d_mean = -0.5*log(1-2*$S-$V) - 0.25*log(1-2*$V);
	my $a = 1/(1-2*$S-$V);
	my $b = 0.5*(1/(1-2*$S-$V)+1/(1-2*$V));
	my $d_stde = sqrt(($a**2*$S+$b**2*$V-($a*$S+$b*$V)**2)/$n);
	my $ci_l = $d_mean - 1.96*$d_stde;
	my $ci_r = $d_mean + 1.96*$d_stde;
#	my $k_mean = 2*log(1-2*$S-$V)/log(1-2*$V)-1;
	return ($d_mean, $d_stde, $ci_l, $ci_r);
}

sub K80MCMC {
	my ($ti, $tv, $n) = @_;
	my $d = ($dRange[1]-$dRange[0])*(rand) + $dRange[0]; # Random initial state of distance
	my $k = ($kRange[1]-$kRange[0])*(rand) + $kRange[0]; # Random initial state of kappa
	my @d; # An array of saved states of distance
	my @k; # An array of saved states of kappa
	foreach (1 .. $NoItera) {
		my $d_old = $d; # Old state of distance
		my $d_new = &NewSta($d_old, $dSlidWin, @dRange);
		my $k_old = $k; # Old state of kappa
		my $k_new = &NewSta($k_old, $kSlidWin, @kRange);
		my $old1 = exp(-4*$d_old/($k_old+2));
		my $old2 = exp(-2*$d_old*($k_old+1)/($k_old+2));
		my $logOld = $ti*log(0.25+0.25*$old1-0.5*$old2) + $tv*log(0.5-0.5*$old1) + ($n-$ti-$tv)*log(0.25+0.25*$old1+0.5*$old2); # Log likelihood of old state
		my $new1 = exp(-4*$d_new/($k_new+2));
		my $new2 = exp(-2*$d_new*($k_new+1)/($k_new+2));
		my $logNew = $ti*log(0.25+0.25*$new1-0.5*$new2) + $tv*log(0.5-0.5*$new1) + ($n-$ti-$tv)*log(0.25+0.25*$new1+0.5*$new2); # Log likelihood of new state
		my $accRat; # Acceptance ratio
		if ($logNew > $logOld) {
			$accRat = 1;
		} else {
			$accRat = exp($logNew-$logOld);
		}
		my $r2 = rand;
		if ($r2 < $accRat) {
			$d = $d_new;
			$k = $k_new;
		}
		my $pri = $_/$Prifreq;
		print "$_\t$d\t$k\n" if $pri == int($pri);
		if ($_ > $Burnin) {
			my $sav = $_/$Savfreq;
			if ($sav == int($sav)) {
				push @d, $d;
				push @k, $k;
			}
		}
	}
	my @d_stat = &statis(@d);
	my @k_stat = &statis(@k);
	return (@d_stat, @k_stat);
}

sub JC69gammaML {
	my ($ti, $tv, $n, $alpha) = @_;
	my $p = ($ti+$tv)/$n;
	my $d_mean = 0.75*$alpha*((1-4/3*$p)**(-1/$alpha)-1);
	my $d_stde = sqrt($p*(1-$p)/$n*(1-4/3*$p)**(-2/($alpha+1))); # There is something wrong.
	my $ci_l = $d_mean - 1.96*$d_stde;
	my $ci_r = $d_mean + 1.96*$d_stde;
	return ($d_mean, $d_stde, $ci_l, $ci_r);
}

sub JC69gammaMCMC {
	my ($ti, $tv, $n, $alpha) = @_;
	my $d = ($dRange[1]-$dRange[0])*(rand) + $dRange[0]; # Random initial state
	my @d; # An array of saved states
	foreach (1 .. $NoItera) {
		my $d_old = $d; # Old state
		my $d_new = &NewSta($d_old, $dSlidWin, @dRange);
		my $old1 = 1/(1+4*$d_old/3/$alpha)**$alpha;
		my $logOld = ($ti+$tv)*log(0.75-0.75*$old1) + ($n-$ti-$tv)*log(0.25+0.75*$old1); # Log likelihood of old state
		my $new1 = 1/(1+4*$d_new/3/$alpha)**$alpha;
		my $logNew = ($ti+$tv)*log(0.75-0.75*$new1) + ($n-$ti-$tv)*log(0.25+0.75*$new1); # Log likelihood of new state
		my $accRat; # Acceptance ratio
		if ($logNew > $logOld) {
			$accRat = 1;
		} else {
			$accRat = exp($logNew-$logOld);
		}
		my $r2 = rand;
		if ($r2 < $accRat) {
			$d = $d_new;
		}
		my $pri = $_/$Prifreq;
		print "$_\t$d\n" if $pri == int($pri);
		if ($_ > $Burnin) {
			my $sav = $_/$Savfreq;
			push @d, $d if $sav == int($sav);
		}
	}
	my @d_stat = &statis(@d);
	return (@d_stat);
}

sub K80gammaML {
	my ($ti, $tv, $n, $alpha) = @_;
	my $S = $ti/$n;
	my $V = $tv/$n;
	my $d_mean = $alpha/2*((1-2*$S-$V)**(-1/$alpha)-1) + $alpha/4*((1-2*$V)**(-1/$alpha)-1);
	my $a = (1-2*$S-$V)**(-1/(1+$alpha));
	my $b = 0.5*((1-2*$S-$V)**(-1/(1+$alpha))+(1-2*$V)**(-1/(1+$alpha)));
	my $d_stde = sqrt(($a**2*$S+$b**2*$V-($a*$S+$b*$V)**2)/$n); # There is something wrong.
	my $ci_l = $d_mean - 1.96*$d_stde;
	my $ci_r = $d_mean + 1.96*$d_stde;
	return ($d_mean, $d_stde, $ci_l, $ci_r);
}

sub K80gammaMCMC {
	my ($ti, $tv, $n, $alpha) = @_;
	my $d = ($dRange[1]-$dRange[0])*(rand) + $dRange[0]; # Random initial state of distance
	my $k = ($kRange[1]-$kRange[0])*(rand) + $kRange[0]; # Random initial state of kappa
	my @d; # An array of saved states of distance
	my @k; # An array of saved states of kappa
	foreach (1 .. $NoItera) {
		my $d_old = $d; # Old state of distance
		my $d_new = &NewSta($d_old, $dSlidWin, @dRange);
		my $k_old = $k; # Old state of kappa
		my $k_new = &NewSta($k_old, $kSlidWin, @kRange);
		my $old1 = 1/(1+4*$d_old/($k_old+2)/$alpha)**$alpha;
		my $old2 = 1/(1+2*$d_old*($k_old+1)/($k_old+2)/$alpha)**$alpha;
		my $logOld = $ti*log(0.25+0.25*$old1-0.5*$old2) + $tv*log(0.5-0.5*$old1) + ($n-$ti-$tv)*log(0.25+0.25*$old1+0.5*$old2); # Log likelihood of old state
		my $new1 = 1/(1+4*$d_new/($k_new+2)/$alpha)**$alpha;
		my $new2 = 1/(1+2*$d_new*($k_new+1)/($k_new+2)/$alpha)**$alpha;
		my $logNew = $ti*log(0.25+0.25*$new1-0.5*$new2) + $tv*log(0.5-0.5*$new1) + ($n-$ti-$tv)*log(0.25+0.25*$new1+0.5*$new2); # Log likelihood of new state
		my $accRat; # Acceptance ratio
		if ($logNew > $logOld) {
			$accRat = 1;
		} else {
			$accRat = exp($logNew-$logOld);
		}
		my $r2 = rand;
		if ($r2 < $accRat) {
			$d = $d_new;
			$k = $k_new;
		}
		my $pri = $_/$Prifreq;
		print "$_\t$d\t$k\n" if $pri == int($pri);
		if ($_ > $Burnin) {
			my $sav = $_/$Savfreq;
			if ($sav == int($sav)) {
				push @d, $d;
				push @k, $k;
			}
		}
	}
	my @d_stat = &statis(@d);
	my @k_stat = &statis(@k);
	return (@d_stat, @k_stat);
}
