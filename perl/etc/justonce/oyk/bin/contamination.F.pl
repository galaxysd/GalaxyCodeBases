#!/usr/bin/env perl
=pod
Author: HU Xuesong @ BGI <huxuesong@genomics.org.cn>, LI Bowen <libowen@genomics.cn>
Version: 1.0.0 @ 20180516
=cut
use strict;
use warnings;
use POSIX;

use FindBin qw($RealBin);
if ($FindBin::VERSION < 1.51) {
	warn "[!]Your Perl is too old, thus there can only be ONE `bsuit` file in your PATH. [FindBin Version: $FindBin::VERSION < 1.51]\n\n"
}
FindBin::again();
use lib "$RealBin/../";
require FGI::GetCPI;

my @Modes = qw(CHIP PCR);
my %Mode = map { $_ => 1 } @Modes;
my $Verbose = 0;

my $theMode = uc shift;
unless (exists $Mode{$theMode}) {
	die "[x]mode can only be:[",join(',',@Modes),"].\n";
}

our @Bases;

my $list = shift;
my $store = shift;
my $output = shift;

my @fams;
my %Fs;
open LI,"<$list" or die($!);
open OUT,">$output" or die($!);
while (my $info = <LI>){
        chomp($info);
        my ($M,$F,$C) = split /\s+/,$info;
	push @fams,"p$C";
        $Fs{$F} = "p$C.F";
}
close LI;

my @fathers = sort keys %Fs;
print OUT "\t";
print OUT join("\t",@fathers),"\n";
foreach my $family (@fams){
	my $Mfile = "$store/$family.M.tsv";
	my $Cfile = "$store/$family.C.tsv";
	print OUT "$family\t";
	my @outputs;
	foreach my $father (@fathers){
		my $total = 0;
		my $mismatch = 0;
		my $Ffile = "$store/$Fs{$father}.tsv";
		my @info = &match($Mfile,$Ffile,$Cfile);
		foreach my $line (@info){
			chomp($line);
			next if ($line =~ /^#/);
			my @data = split /\t/,$line;
			next unless (defined $data[7]);
			$total++;
			if ($data[7] == 0.0001){
				$mismatch++;
			}
		}
		my $outinfo = join "/",$mismatch,$total;
		push @outputs,$outinfo;
	}
	print OUT join("\t",@outputs),"\n";
}
close OUT;


########################################################
sub deBayes($) {
	my $p = $_[0];
	my %Dep;
	for my $i (1 .. $#$p) {
		$Dep{$i-1} = $p->[$i];
	}
	#ddx %Dep;
	my @dKeys = sort { $Dep{$b} <=> $Dep{$a} } keys %Dep;
	if ( @dKeys>1 and $Dep{$dKeys[1]} >= $Dep{$dKeys[0]} * 0.02) {	# 2%
		my @rKeys = sort {$a<=>$b} @dKeys[0,1];
		my $gt = join('/',$Bases[$rKeys[0]],$Bases[$rKeys[1]]);
		$p->[0] = $gt;
	} elsif (@dKeys>1 && ($Dep{$dKeys[1]} < $Dep{$dKeys[0]} * 0.02) && ($Dep{$dKeys[1]} > $Dep{$dKeys[0]} * 0.001)){
		$p->[0] = "NA";
	} elsif (@dKeys == 1 or ($Dep{$dKeys[1]} <= $Dep{$dKeys[0]} * 0.001)){
		my $gt = join('/',$Bases[$dKeys[0]],$Bases[$dKeys[0]]);
		$p->[0] = $gt;
	}
}
sub deBayes2($) {
	my $p = $_[0];
	my %Dep;
	for my $i (1 .. $#$p) {
		$Dep{$i-1} = $p->[$i];
	}
	#ddx %Dep;
	my @dKeys = sort { $Dep{$b} <=> $Dep{$a} } keys %Dep;
	if ( @dKeys>1 and $Dep{$dKeys[1]}  >= $Dep{$dKeys[0]} * 0.1 ) {	# 10%
		my @rKeys = sort {$a<=>$b} @dKeys[0,1];
		my $gt = join('/',$Bases[$rKeys[0]],$Bases[$rKeys[1]]);
		$p->[0] = $gt;
	} elsif (@dKeys>1 && ($Dep{$dKeys[1]}  > $Dep{$dKeys[0]} * 0.01) && ($Dep{$dKeys[1]}  < $Dep{$dKeys[0]} * 0.1)){
		$p->[0] = "NA";
	}
}

sub getBolsheviks(@) {
	my $type = shift;
	my @dat = map { [split /[;,]/,$_] } @_;
	if ($type) {
		deBayes($_) for @dat;
	} else {
		deBayes2($_) for @dat;
	}
	#ddx \@dat;

	my (%GT);
	for (@dat) {
		++$GT{$_->[0]};
	}
	if (defined $GT{NA}){
		return ["NA",0,"NA"];
	}
	my $Bolsheviks = (sort {$GT{$b} <=> $GT{$a}} keys %GT)[0];
	my @GTdep;
	for (@dat) {
		next if $_->[0] ne $Bolsheviks;
		for my $i (1 .. $#$_) {
			$GTdep[$i-1] += $_->[$i];
		}
	}
	my @GTs = split /[\/|]/,$Bolsheviks;
	my $Hom = 0;
	$Hom = 1 if $GTs[0] eq $GTs[1];
	#ddx $Bolsheviks,$Hom,\@GTdep,@dat if (keys %GT)>1;
	return [$Bolsheviks,$Hom,\@GTdep];
}

sub getequal(@) {
	my $type = shift;
	my @dat = map { [split /[;,]/,$_] } @_;
#ddx \@dat;
	if ($type) {
		deBayes($_) for @dat;
	} else {
		deBayes2($_) for @dat;
	}
	my (%GT);
	for (@dat) {
		++$GT{$_->[0]};
	}
	for (values %GT) {
		return 1 if $_ == 2;	# 两个样品GT一致
	}
	return 0;
}

sub tstat(%) {
	my %d = @_;
	unless ($d{'n'}) {
		return ('NA','NA');
	}
	my $mean1 = $d{'x'}/$d{'n'};
	my $mean2 = $d{'y'}/$d{'n'};
	my $std1 = sqrt($d{'xx'}/$d{'n'} - $mean1*$mean1);
	my $std2 = sqrt($d{'yy'}/$d{'n'} - $mean2*$mean2);
	my $srt1 = join(' ± ',$mean1,$std1);
	my $srt2 = join(' ± ',$mean2,$std2);
	return ($srt1,$srt2);
}
sub printExp($) {
	my $lnV = $_[0]/log(10);
	my $lnInt = floor($lnV);
	my $lnExt = $lnV - $lnInt;
	my $prefix = exp($lnExt*log(10));
	my $str = join('e',$prefix,$lnInt);
	return $str;
}

sub get_locus{
        my $file = shift;
        my %temp;
        open TE,"<$file" or die($!);
        while (<TE>){
                chomp;
                my @data = split /\t/,$_;
                $temp{$data[0]}{$data[1]}++;
        }
        close TE;
        return %temp;
}

sub reshape(@) {
        my $bases = shift;
        my @TBases = @$bases;
        for (@_){
                my %tempValue;
                my @dep = split /[;,]/,$_;
                my @newDep;
                for my $i(1..scalar @dep - 1){
                        $tempValue{$TBases[$i - 1]} = $dep[$i];
                }
                foreach my $allele (@Bases){
                        if (defined $tempValue{$allele}){
                                push @newDep,$tempValue{$allele};
                        }else{
                                push @newDep,0;
                        }
                }
                $_ = join(";",$dep[0],join(",",@newDep));
        }
}

sub match {
	my ($mother,$father,$child) = @_;
	my @temp_out;

	my %locusM = get_locus($mother);
	my %locusF = get_locus($father);
	my %locusC = get_locus($child);
	my %need;
	foreach my $chr (keys %locusM){
	        foreach my $pos (keys %{$locusM{$chr}}){
	                if (defined $locusF{$chr}{$pos} && $locusC{$chr}{$pos}){
	                        $need{$chr}{$pos}++;
	                }
	        }
	}

open FM,'<',$mother or die "[x]Mom: $!\n";
open FF,'<',$father or die "[x]Dad: $!\n";
open FC,'<',$child or die "[x]Child: $!\n";

my ($logcpi,$spe,$trioN,$lFC,$lFF,$lFM)=(0,0,0);

my %check_dup;
while (<FM>) {
	chomp;
	my @datM = split /\t/;
	next unless (defined $need{$datM[0]}{$datM[1]});
	$check_dup{$datM[0]}{$datM[1]}++;
	next if ($check_dup{$datM[0]}{$datM[1]} > 1);
	my (@datF,@datC);
        while ($lFF = <FF>){
                chomp($lFF);
                @datF = split /\t/,$lFF;
                if ($datF[0] eq $datM[0] && $datF[1] eq $datM[1]){
                        last;
                }
        }
        while ($lFC = <FC>){
                chomp($lFC);
                @datC = split /\t/,$lFC;
                if ($datC[0] eq $datM[0] && $datC[1] eq $datM[1]){
                        last;
                }
        }
        unless ($datF[0] eq $datM[0] && $datF[1] eq $datM[1] && $datC[0] eq $datM[0] && $datC[1] eq $datM[1]){
                last;
        }
	#my ($chr,undef,$bases,$qual,@data) = split /\t/;
	next if $datM[3] !~ /\d/ or $datM[3] < 100;
	next if $datF[3] !~ /\d/ or $datF[3] < 100;
	next if $datC[3] !~ /\d/ or $datC[3] < 100;
	die if $datM[0] ne $datC[0] or $datF[0] ne $datC[0];
	my @tM = splice @datM,4;
	my @tF = splice @datF,4;
	my @tC = splice @datC,4;
	@Bases = split /,/,$datM[2];	# $bases = ref,alt
        unless ($datM[2] eq $datF[2] && $datM[2] eq $datC[2]){
                my @MBases = split /,/,$datM[2];
                my @FBases = split /,/,$datF[2];
                my @CBases = split /,/,$datC[2];
                unless ($MBases[0] eq $FBases[0] && $MBases[0] eq $CBases[0]){
                        next;
                }
                my %alts;
                for my $i (1..scalar @MBases - 1){$alts{$MBases[$i]}++;}
                for my $i (1..scalar @FBases - 1){$alts{$FBases[$i]}++;}
                for my $i (1..scalar @CBases - 1){$alts{$CBases[$i]}++;}
                @Bases = sort keys %alts;
                unshift @Bases,$MBases[0];
                reshape(\@MBases,@tM);
                reshape(\@FBases,@tF);
                reshape(\@CBases,@tC);
        }
	next if $Bases[1] eq '.';
	next if "@tM @tF @tC" =~ /\./;

	my $retM = getBolsheviks(0,@tM);
	my $retF = getBolsheviks(0,@tF);
	#ddx $retM if $retM->[1];

	my $check_dep = 1;
	for (@tM){
		my @info = split /[;,]/,$_;
		my $sum;
		for my $i(1..scalar @info - 1){
			$sum += $info[$i];
		}
		if ($sum > 50){
			$check_dep *= 1;
		}else{
			$check_dep *= 0;
		}
	}
	for (@tF){
		my @info = split /[;,]/,$_;
		my $sum;
		for my $i(1..scalar @info - 1){
			$sum += $info[$i];
		}
		if ($sum > 50){
			$check_dep *= 1;
		}else{
			$check_dep *= 0;
		}
	}
	next if ($check_dep == 0);

	#T/T;6,2245      C/C;1698,0
	#print "> @tM , @tF , @tC\n@datM\n";
	#my $retM = getBolsheviks(0,@tM);
	next unless $retM->[1];
	#my $retF = getBolsheviks(0,@tF);
	next if ($retM->[0] eq "NA" or $retF->[0] eq "NA");
	#ddx $retM,$retF;
	my $xx = getequal(0,@tM);
	my $yy = getequal(0,@tF);
	my $zz = getequal(1,@tC);
	my $t=$xx*$yy*$zz;
	my $REP = 1;
	if ($theMode eq 'PCR') {
		next unless $t;
	} else {
		$REP = 2 if $t;
	}
	my @sdatC = map { [split /[;,]/,$_] } @tC;
	my @GTdepC;
	for (@sdatC) {
		for my $i (1 .. $#$_) {
			$GTdepC[$i-1] += $_->[$i];
		}
	}
	my %CntM;
	my @GTM = @{$retM->[2]};
	for my $i (0 .. $#GTM) {
		$CntM{$i} = $GTM[$i];
	}
	my ($x,$y) = sort { $CntM{$b} <=> $CntM{$a} } keys %CntM;
	#ddx [$x,$y,$Bases[$x],$Bases[$y]],\@sdatC,\@GTdepC;
	my ($n12,$n22);
	my $n11 = $retM->[2]->[$x];
	my $n21 = $GTdepC[$x];
	if (defined $y) {
		$n12 = $retM->[2]->[$y];
		$n22 = $GTdepC[$y];
	} else {
		$n12 = 0;
		$n22 = 0;
	}
	next unless defined $n22;
	if ($theMode eq 'PCR') {
		next if ($n21+$n22) < 200;
	} elsif ($theMode eq 'CHIP') {
		next if ($n21+$n22) < (100 * $REP);
	}
	my $GTtC;
	$GTtC = join('/',$Bases[$x],$Bases[$x]);
	my $Cdep = $n21 + $n22;

	my $retC = getBolsheviks(1,@tC);
	#ddx $retM,$retF,$retC;
	next if ($retC->[0] eq "NA");
	my @fgeno=split /\//,$retF->[0];
	my @mgeno=split /\//,$retM->[0];
	my @cgeno=split /\//,$retC->[0];

	if ($theMode eq 'PCR') {
		next if $fgeno[0] eq $fgeno[1] and $mgeno[0] eq $mgeno[1] and $mgeno[0] eq $fgeno[0] and (($retM->[2][0]>0 and $retM->[2][1]>0) or ($retF->[2][0]>0 and $retF->[2][1]>0));
	} elsif ($theMode eq 'CHIP') {
		next if $mgeno[0] eq $mgeno[1] && $cgeno[0] eq $cgeno[1] && $mgeno[0] eq $cgeno[0];
	} else {
		die;
	}

	my @mnum=@{$retM->[2]};
	my @fnum=@{$retF->[2]};	
	my $resM = join(';',$retM->[0],join(',',@{$retM->[2]}));
	my $resF = join(';',$retF->[0],join(',',@{$retF->[2]}));
	my $resC = join(';',$retC->[0],join(',',@{$retC->[2]}));
#	my $resC = join(';',$GTtC,join(',',@GTdepC),$twotailedFisher,
#						$retC->[0],join(',',@{$retC->[2]})
#					);
	my $cret = getcpi(@datM,$resM,$resF,$resC);
	#ddx $cret;
	$logcpi += log($cret->[0]);
	$spe += log(1-$cret->[1]);
	my $tempout = join("\t",@datM,$resM,$resF,$resC,@$cret,$logcpi/log(10),$spe/log(10)),"\n";
	push @temp_out,$tempout;
}

close FM; close FF; close FC;
return @temp_out;
}
