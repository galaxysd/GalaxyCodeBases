#!/usr/bin/env perl
use strict;
use warnings;

use Data::Dump qw(ddx);

my $Usage = "Usage: $0 <vcf.gz> <out prefix>\n";
die $Usage if @ARGV < 1;

my ($filename,$outp) = @ARGV;
my $cmd = "bcftools view -U -m2 -i '%QUAL>=40 & MIN(FMT/GQ)>20' -v snps $filename |bcftools view -e 'FMT/DP=\".\"'| bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%QUAL[\t%TGT:%AD:%GQ]\n'";
my @SampleIDs;
open(SID,'-|',"bcftools query -l $filename") or die "Error opening [$filename]: $!\n";
while(<SID>) {
	chomp;
	push @SampleIDs,$_;
}
close SID;
ddx \@SampleIDs;

my $theKiller = 'T3C';
my $theVictim = 'T10C';
my $Mixture = 'mixed';

sub mergeGT($) {
	my $hashref = $_[0];
	my @GTs = sort keys %{$hashref};
	s/\./0/ for @GTs;
	if (scalar @GTs ==1) {
		push @GTs,@GTs;
	}
	return join(' ',@GTs);
}
my $HetR = 0.2;
sub doCal($$) {
	my ($mix,$vit) = @_;
	my @mixGT = keys %{$mix};
	my @vitGT = keys %{$vit};
	my %mixGT = map { $_ => 1 } @mixGT;
	my ($sameGT,$vitGT);
	for (@vitGT) {
		if (exists $mixGT{$_}) {
			delete $mixGT{$_};
		}
	}
	my @extraBase = keys %mixGT;
	#ddx ($mix,$vit);
	if (scalar @extraBase > 1) {
		warn "[@extraBase]";
		return -1;
	}
	my $Extradepth = $mix->{$extraBase[0]};
	my $Totaldepth = 0;
	for (values %{$mix}) {
		$Totaldepth += $_;
	}
	my $depth0 = $Extradepth / $Totaldepth;
	my $ret1 = $depth0 / (1-0.5*$HetR);
	my $ret2 = $depth0 * (2/3) / 0.6;
	#ddx [$ret1,$ret2,$Extradepth,$Totaldepth];
	return ($ret1,$Extradepth,$Totaldepth);
}

open O,'>',"$outp.tsv" or die "[$outp.tsv]:$!\n";
my (%Calcu,%Results);
open(IN,"-|",$cmd) or die "Error opening [$filename]: $!\n";
while (<IN>) {
	chomp;
	my ($Chrom,$Pos,$RefAlt,$Qual,@GTs) = split /\t/,$_;
	next if $Chrom =~ /_/;
	#print "$_\n";
	my (%GTs,%GTstr,%GTped);
	@GTstr{@SampleIDs} = @GTs;
	#ddx \%GTstr;
	for my $k (keys %GTstr) {
		my ($sGT,$sAD,$sGQ) = split /:/,$GTstr{$k};
		my @aGT = split /[|\/]/,$sGT;
		my @aAD = split /\,/,$sAD;
		my %tmp;
		@tmp{@aGT} = @aAD;
		$GTs{$k} = \%tmp;
	}
	#ddx \%GTs;
	my %Killer = %{$GTs{$theKiller}};
	my %Victim = %{$GTs{$theVictim}};
	my %Mixture = %{$GTs{$Mixture}};
	next if scalar(keys %Mixture) != 2;
	next if scalar(keys %Victim) != 1;
	my $gt1 = mergeGT(\%Mixture);
	next if $gt1 =~ /0/;
	next if length($gt1)>3;
	my $gt2 = mergeGT(\%Victim);
	next if $gt2 =~ /0/;
	next if length($gt2)>3;
	#ddx [\%Killer,\%Victim,\%Mixture,$gt1,$gt2];
	my ($ret,$Extradepth,$Totaldepth) = doCal(\%Mixture,\%Victim);
	#my $ret = doCal(\%Mixture,\%Killer);
	next if $ret == -1;
	++$Results{int($ret*100)/100};
	$Calcu{'S'} += $ret;
	$Calcu{'SS'} += $ret*$ret;
	++$Calcu{'N'};
	#ddx \%Results,\%Calcu;
	print O join("\t",$Chrom,$Pos,$Qual,$gt1,$gt2,$ret,$Extradepth,$Totaldepth),"\n";
}

my $mean = $Calcu{'S'} / $Calcu{'N'};
my $std = sqrt($Calcu{'SS'}/$Calcu{'N'} - $mean*$mean);
my $var = $std/$mean;

print O "# $mean ± $std, $var\n";
close O;
print "$mean ± $std, $var\n";

__END__
my $theKiller = 'T3C';
my $theVictim = 'T10C';

T3C: 85,982,846,136 bp
T10C:18,277,246,581 bp

K% = 82%

#   { N => 3859, S => 1342.27179667765, SS => 676.187587538486 },
# )
0.347828918548237 ± 0.232891755122676, 0.669558345219586

#   { N => 370586, S => 175277.892782687, SS => 98608.7057863698 },
# )
0.472974944500566 ± 0.205872025106518, 0.43527046728428
