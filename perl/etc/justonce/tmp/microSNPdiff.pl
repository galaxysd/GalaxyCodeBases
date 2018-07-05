#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my $SampleList = 'samples.lst';
my $SNPtsv = 'micro.snp.tsv';
my $Out1 = 'micro.consisted.tsv';
my $Out2 = 'micro.same.tsv';
my $Out3 = 'micro.GT.tsv';
my $Out4 = 'micro.diff.tsv';
my $Out5 = 'micro.Tsame.tsv';
my $Out6 = 'micro.Tdiff.tsv';

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

my (%SameCntSum,%SameGTdat,%SameCnt,%TopLoc);
open O,'>',$Out1 or die $?;
open P,'>',$Out3 or die $?;
print O join("\t",'Chr','Pos','Consisted',@UIDs),"\n";
print P join("\t",'Chr','Pos','Consisted',@UIDs),"\n";
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
	#my %SameCntLoc;
	for my $uid (@UIDs) {
		my @SmpGTpKeys = sort keys %{$SmpGTp{$uid}};
		$SameCnt{$Loc}{$uid} = 0;
		if (@SmpGTpKeys == 1) {
			$SameCnt{$Loc}{$uid}  = 1;
			$SameCnt{$Loc}{'\t'} += 1;
		}
		$SameGTdat{$Loc}{$uid} = join(',',sort @SmpGTpKeys);
	}
	#ddx \%SameCnt;
	++$SameCntSum{$SameCnt{$Loc}{'\t'}};
	print O join("\t",$Chr,$Pos,$SameCnt{$Loc}{'\t'},(map {$SameCnt{$Loc}{$_}} @UIDs)),"\n";
	print P join("\t",$Chr,$Pos,$SameCnt{$Loc}{'\t'},(map {$SameGTdat{$Loc}{$_}} @UIDs)),"\n";
	if ($SameCnt{$Loc}{'\t'} == scalar @UIDs) {
		++$TopLoc{$Loc};
	}
}
close I;
close P;

print O "\n# Total ",scalar @UIDs," Individuals.\n\n# Summary\n";
for my $i (sort {$a <=> $b} keys %SameCntSum) {
	print O "# $i:",$SameCntSum{$i},"\n";
}
close O;

#my %SameDiffSum;
open O,'>',$Out2 or die $?;
open E,'>',$Out4 or die $?;
open X,'>',$Out5 or die $?;
open Y,'>',$Out6 or die $?;
print O join("\t",'=',@UIDs),"\n";
print E join("\t",'=',@UIDs),"\n";
print X join("\t",'=',@UIDs),"\n";
print Y join("\t",'=',@UIDs),"\n";
for my $i (0 .. $#UIDs) {
	my $uidi = $UIDs[$i];
	print O $uidi;
	print E $uidi;
	print X $uidi;
	print Y $uidi;
	for my $j (0 .. $i) {
		my $uidj = $UIDs[$j];
		my ($cntij,$cntijE,$Tcntij,$TcntijE) = (0,0,0,0);
		for my $Loc (keys %SameGTdat) {
			my $GTi = $SameGTdat{$Loc}{$uidi};
			my $GTj = $SameGTdat{$Loc}{$uidj};
			next unless $SameCnt{$Loc}{$uidi};
			next unless $SameCnt{$Loc}{$uidj};
			if ($GTi eq $GTj) {
				#next if $GTi =~ /,/;
				++$cntij;
			} else {
				++$cntijE;
			}
			if ($TopLoc{$Loc}) {
				if ($GTi eq $GTj) {
					++$Tcntij;
				} else {
					++$TcntijE;
				}
			}
		}
		#++$SameDiffSum{$cntij};
		print O "\t$cntij";
		print E "\t$cntijE";
		print X "\t$Tcntij";
		print Y "\t$TcntijE";
	}
	print O "\n";
	print E "\n";
	print X "\n";
	print Y "\n";
}

close O;
close E;
close X;
close Y;

#ddx \%SameDiffSum;
__END__

ddx \%Samples;
ddx \%ID2Sample;
ddx \%SampleCnts;
print scalar keys %SampleCnts,"\n";

