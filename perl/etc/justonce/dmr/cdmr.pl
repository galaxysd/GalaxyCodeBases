#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180614
=cut
use strict;
use warnings;
use Cwd qw(abs_path cwd);
use Parse::CSV;
use List::Util qw[min max];
use Statistics::TTest;
#use lib '.';
use Data::Dump qw(ddx);

my ($bedfn,$cgLfn,$cgRfn,$outfn) = @ARGV;
if (@ARGV < 4) {
	die "Usage: $0 <bed> <Left.cg> <Right.cg> <out.tsv>\n";;
}

my $cgLin = Parse::CSV->new(file => $cgLfn, names => 1, sep_char => "\t");
my $cgRin = Parse::CSV->new(file => $cgRfn, names => 1, sep_char => "\t");

open B,'<',$bedfn or die $?;
my (%Dbed,%Rbed);
while (<B>) {
	my @L = split /\t/;
	next unless $L[2];
	push @{$Dbed{$L[0]}},[$L[1],$L[2]];
	my @rng = $L[1] .. $L[2];
	for (@rng) {
		++$Rbed{$L[0]}{$_};
	}
}
#ddx %Dbed;
close B;

sub loadCG($$) {
	my ($rngD,$csvh) = @_;
	my %Dcg;
	while ( my $value = $csvh->fetch ) {
		#ddx $value;
		if (exists $rngD->{ $value->{'#CHROM'} }{ $value->{'POS'} }) {
			my ($cgCnt,$cgC,$cgW,$flag)=(0,0,0,0);
			if ($value->{'Crick-COVERAGE'} ne '.') {
				if ($value->{'Crick-COVERAGE'} >= 5) {
					$cgC = $value->{'Crick-METH'}/$value->{'Crick-COVERAGE'};
					$cgCnt += 1;
					$flag |= 1;
				}
			}
			if ($value->{'Watson-COVERAGE'} ne '.') {
				if ($value->{'Watson-COVERAGE'} >= 5) {
					$cgW = $value->{'Watson-METH'}/$value->{'Watson-COVERAGE'};
					$cgCnt += 1;
					$flag |= 2;
				}
			}
			if ($cgCnt) {
				$Dcg{$value->{'#CHROM'}}{$value->{'POS'}} = [$cgCnt/2, ($cgC+$cgW)/$cgCnt,$flag];
			}
		} else {
			next;
		}
	}
	return \%Dcg;
}

my %DcgL = %{ loadCG(\%Rbed,$cgLin) };
#ddx \%DcgL;
my %DcgR = %{ loadCG(\%Rbed,$cgRin) };
#ddx \%DcgR;

sub getCG($$$) {
    my ($dat,$st,$ed)=@_;
    my %ret;
    for my $p ($st .. $ed) {
        if (exists $dat->{$p}) {
            $ret{$p} = $dat->{$p};
        }
    }
    return %ret;
}

open O,'>',$outfn or die $?;
print O "# L:[$cgLfn], R:[$cgRfn].\n";
open B,'<',$bedfn or die $?;
while (<B>) {
	my @L = split /\t/;
	next unless $L[2];
	print O join("\t",@L[0..2]),"\t";
	my %cgLr = getCG($DcgL{$L[0]},$L[1],$L[2]);
	my %cgRr = getCG($DcgR{$L[0]},$L[1],$L[2]);
	#ddx \%cgLr;
	my (%rPoses,%cgR);
	for (keys %cgLr,keys %cgRr) {
        ++$rPoses{$_};
	}
	for (keys %rPoses) {
        my $dL = [0,0,0];
        my $dR = [0,0,0];
        $dL = $cgLr{$_} if exists $cgLr{$_};
        $dR = $cgRr{$_} if exists $cgRr{$_};
        $cgR{$_} = [@$dL,@$dR];
	}
	#ddx \%cgR;
	(%cgLr,%cgRr)=();
	my ($shared,$left,$right,$unioned)=(0,0,0,0);
	my ($Lsum,$Lcnt,$Rsum,$Rcnt)=(0,0,0,0);
	my (@tL,@tR);
	for (keys %cgR) {
        my @d = @{ $cgR{$_} };
        $shared += min($d[0],$d[3]);
        $left += $d[0];
        $right += $d[3];
        my $flag = $d[2] | $d[5];
        if ($flag==3) {
            $unioned += 1;
            push @tL,$d[1];
            push @tR,$d[4];
        } elsif ($flag==1 or $flag==2) {
            $unioned += 0.5;
        }
        #ddx $flag;
        $Lsum += $d[1]*$d[0]; $Lcnt += $d[0];
        $Rsum += $d[4]*$d[3]; $Rcnt += $d[3];
	}
    my $p='NA';
    if (@tL>1 and @tR>1) {
        my $ttest = new Statistics::TTest;
        $ttest->load_data(\@tL,\@tR);
        $p = $ttest->{t_prob};
    }
    my ($avgL,$avgR)=qw(NA NA);
    $avgL = $Lsum/$Lcnt if $Lcnt;
    $avgR = $Rsum/$Rcnt if $Rcnt;
    print O join(',',$shared,$left,$right,$unioned),"\t",join("\t",$avgL,$avgR,$p),"\n";
}
close B;
close O;

__END__
./cdmr.pl gencode.v30.annotation.bed.h B7B.cg.h Normal.cg.h o.tsv
