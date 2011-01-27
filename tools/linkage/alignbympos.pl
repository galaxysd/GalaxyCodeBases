#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
#use GalaxyXS::ChromByte;
use Data::Dump qw(ddx);

my $maxkbRchg=0.05;
my $scaffmainR=0.7;

unless (@ARGV > 0) {
    print "perl $0 <markerpos scaffold> <markerpos chr> <out>\n";
    exit 0;
}

my ($scaff,$chrf,$outf)=@ARGV;
my $opt_v=1;
my (%MarkerDat,%ScaffAlign);

open S,'<',$scaff or die "Error:[$scaff] $!\n";
#Markerid       MarkercM        Sid     pos     strand  Pidentity       E       BTOP
while (<S>) {
	next if /^#/;
	chomp;
	my ($Qid,$mcM,$Sid,$pos,$strand,$Pidentity,$E,$BTOP)=split /\t/;
	$MarkerDat{$Qid}=[[$Sid,$pos,$strand]];
}
close S;
open C,'<',$chrf or die "Error:[$chrf] $!\n";
while (<C>) {
	next if /^#/;
	chomp;
	my ($Qid,$mcM,$Sid,$pos,$strand,$Pidentity,$E,$BTOP)=split /\t/;
	next unless exists $MarkerDat{$Qid};
	push @{$MarkerDat{$Qid}},[$Sid,$pos,$strand];
}
close C;
#ddx \%MarkerDat;
for my $markid (keys %MarkerDat) {
	my ($Scaff,$Spos,$Ss)=@{$MarkerDat{$markid}->[0]};
	next if @{$MarkerDat{$markid}} != 2;
	my ($Chr,$Cpos,$Cs)=@{$MarkerDat{$markid}->[1]};
	push @{$ScaffAlign{$Scaff}{$Chr}{$Ss*$Cs}},[$Spos,$Cpos];
}
%MarkerDat=();
#ddx \%ScaffAlign;
#   scaffold9788  => {
#                      Chr10 => {
#                        1 => [
#                               [3465, 10475529],
#                               [2231, 10474295],
#                               [2636, 10474700],
#                               [3397, 10475461],
#                               [2074, 10474138],
#                               [2211, 10474275],
#                               [2241, 10474305],
#                               [2692, 10474756],
#                               [3325, 10475389],
#                               [2072, 10474136],
#                               [3516, 10475580],
#                             ],
#                      },
#                    },

my $valuecount=1;
sub addavg ($$) {
	my ($base,$add)=@_;
	$base=($base*$valuecount+$add)/($valuecount+1);
	return $base;
}
my $inScaff = keys %ScaffAlign;
my $outScaff=$inScaff;
for my $Scaff (keys %ScaffAlign) {
	my ($posSum,$posMax,$datAref,@datMPath)=(0,0);
	for my $Chr (keys %{$ScaffAlign{$Scaff}}) {
		for my $Strand (keys %{$ScaffAlign{$Scaff}{$Chr}}) {
			$datAref=$ScaffAlign{$Scaff}{$Chr}{$Strand};
			$posSum += @$datAref;
			if ($posMax < @$datAref) {
				$posMax = @$datAref;
				@datMPath=($Chr,$Strand);
			}
		}
	}
	$datAref=$ScaffAlign{$Scaff}{$datMPath[0]}{$datMPath[1]};
	if ($opt_v && $posMax != $posSum) {
		print '-'x75,"\n[!]$posMax\t$posSum\n";
		ddx $ScaffAlign{$Scaff};
	}
	if ($posMax < $posSum * $scaffmainR) {
		delete $ScaffAlign{$Scaff};
		--$outScaff;
		print "Deleted !\n" if $opt_v;
		next;
	} elsif ($posMax>1) {
		my @dat=@{$datAref};
		my ($S0pos,$C0pos)=@{shift @dat};
		my ($S1pos,$C1pos)=@{shift @dat};
	warn $Scaff,"\n";
	ddx $ScaffAlign{$Scaff} if $S1pos==$S0pos;
		my $k=($C1pos-$C0pos)/($S1pos-$S0pos);
		my $b=($S1pos*$C0pos-$S0pos*$C1pos)/($S1pos-$S0pos);
		my ($i,%newDat)=(1);
		$valuecount=1;
		$newDat{$i}=[ [$S0pos,$C0pos],[$S1pos,$C1pos] ];
		for (@dat) {
			my ($Spos,$Cpos)=@$_;
			my $k1=($C1pos-$Cpos)/($S1pos-$Spos);
			my $b1=($S1pos*$Cpos-$Spos*$C1pos)/($S1pos-$Spos);
			if ( ($k1-$k)/$k > $maxkbRchg or ($b1-$b)/$b > $maxkbRchg ) {
				++$i;
				$valuecount=0;
				$k=$k1; $b=$b1;
				$C1pos=$Cpos; $S1pos=$Spos;
			} else {
				$k=&addavg($k,$k1); $b=&addavg($b,$b1);
				$C1pos=&addavg($C1pos,$Cpos); $S1pos=&addavg($S1pos,$Spos);
			}
			push @{$newDat{$i}},[$Spos,$Cpos];
			++$valuecount;
		}
		;
	} else {
		;
	}
	$ScaffAlign{$Scaff}=[@datMPath,$datAref];
}
warn "[!]Scaffold Count: $inScaff -> $outScaff\n";
#ddx \%ScaffAlign;
#   scaffold94519 => [
#                      "Chr10",
#                      1,
#                      [
#                        [1322, 16088080],
#                        [2610, 16089368],
#                        [1271, 16088029],
#                        [1738, 16088496],
#                        [1597, 16088355],
#                        [1777, 16088535],
#                        [1855, 16088613],
#                      ],
#                    ],

__END__
./alignbympos.pl markerpos/m2sChr10.pos.f markerpos/m2cChr10.pos.f t.out
