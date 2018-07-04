#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my $SampleList = 'samples.lst';
my $SNPtsv = 'micro.snp.tsv';
my $Out1 = 'micro.consisted.tsv';

my (%Samples,@Samples,@UIDs,%ID2Sample,%SampleCnts);
open I,'<',$SampleList or die $?;
while (<I>) {
	chomp;
	my ($id,$rep) = /([A-Za-z]+)(\w+)/;
	push @Samples,$_;
	push @{$Samples{$1}},$_;
	$ID2Sample{$_} = $1;
	++$SampleCnts{$1};
}
close I;
@UIDs = sort keys %Samples;

my (%SameCntSum);
open O,'>',$Out1 or die $?;
print O join("\t",'Chr','Pos','Consisted',@UIDs),"\n";
open I,'<',$SNPtsv or die $?;
while (<I>) {
	chomp;
	my ($Chr,$Pos,$GTstr,$Qual,@Dat) = split /\t/;
	my $Loc = "$Chr\t$Pos";
	my @GTs = split /,/,$GTstr;
	my @DatA = map { (split /;/,$_)[0] } @Dat;
	#ddx [$Chr,$Pos,\@GTs,$Qual,\@DatA];
	my %SmpGTp;
	for my $i (0 .. $#Samples) {
		my $sid = $Samples[$i];
		my $uid = $ID2Sample{$sid};
		my $GT = $DatA[$i];
		++$SmpGTp{$uid}{$GT};
	}
	#ddx \%SmpGTp;
	my %SameCntLoc;
	for my $uid (@UIDs) {
		my @SmpGTpKeys = keys %{$SmpGTp{$uid}};
		$SameCntLoc{$uid} = 0;
		if (@SmpGTpKeys == 1) {
			$SameCntLoc{$uid}  = 1;
			$SameCntLoc{'\t'} += 1;
		}
	}
	#ddx \%SameCnt;
	++$SameCntSum{$SameCntLoc{'\t'}};
	print O join("\t",$Chr,$Pos,$SameCntLoc{'\t'},(map {$SameCntLoc{$_}} @UIDs)),"\n";
}
close I;

print O "\n# Total ",scalar @UIDs," Individuals.\n\n# Summary\n";
for my $i (sort {$a <=> $b} keys %SameCntSum) {
	print O "# $i:",$SameCntSum{$i},"\n";
}

close O;


__END__

ddx \%Samples;
ddx \%ID2Sample;
ddx \%SampleCnts;
print scalar keys %SampleCnts,"\n";

