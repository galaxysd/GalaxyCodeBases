#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
use GalaxyXS::ChromByte;
use Data::Dump qw(ddx);

my $scaffKnot1R=0.1;
my $maxKdiffto1=0.1;
my $maxIndelR=0.002;
my $scaffmainR=0.7;

unless (@ARGV > 0) {
    print "perl $0 <markerpos scaffold> <markerpos chr> <chr Nzone> <chr nfo> <out>\n";
    exit 0;
}

my ($scaff,$chrf,$NZchrf,$ChrNFO,$outf)=@ARGV;
my $opt_v=1;
my (%MarkerDat,%ScaffAlign,%NZone);

open S,'<',$ChrNFO or die "Error:[$ChrNFO] $!\n";
while (<S>) {
	next if /^#/;
	#chomp;
	my ($chrid,$len)=split /\t/;
	$NZone{$chrid}=initchr($len);
}
close S;
open N,'<',$NZchrf or die "Error:[$NZchrf] $!\n";
while (<N>) {
	next if /^#/;
	#chomp;
	my ($chrid,$start,$end)=split /\t/;
	die "[x]Wrong ChrID as [$chrid] !\n" unless exists $NZone{$chrid};
	setbases($NZone{$chrid},$start,$end,1);
}
close N;
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
	# There can be markers with the same sequence >_<
	#push @{$ScaffAlign{$Scaff}{$Chr}{$Ss*$Cs}},[$Spos,$Cpos];
	unless (exists $ScaffAlign{$Scaff}{$Chr}{$Ss*$Cs}{$Spos}) {
		$ScaffAlign{$Scaff}{$Chr}{$Ss*$Cs}{$Spos}=$Cpos;
	} else {
		$ScaffAlign{$Scaff}{$Chr}{$Ss*$Cs}{$Spos}=undef;
	}
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
=pod
my $valuecount=1;
sub addavg ($$) {
	my ($base,$add)=@_;
	$base=($base*$valuecount+$add)/($valuecount+1);
	return $base;
}
=cut
my $inScaff = keys %ScaffAlign;
my $outScaff=$inScaff;
for my $Scaff (keys %ScaffAlign) {
	my ($posSum,$posMax,$datAref,@datMPath)=(0,0);
	for my $Chr (keys %{$ScaffAlign{$Scaff}}) {
		for my $Strand (keys %{$ScaffAlign{$Scaff}{$Chr}}) {
			$datAref=[];
			for my $Spos (keys %{$ScaffAlign{$Scaff}{$Chr}{$Strand}}) {
				my $Cpos=$ScaffAlign{$Scaff}{$Chr}{$Strand}{$Spos};
				next unless defined $Cpos;
				push @$datAref,[$Spos,$Cpos];
			}
			#$datAref=$ScaffAlign{$Scaff}{$Chr}{$Strand};
			$ScaffAlign{$Scaff}{$Chr}{$Strand}=$datAref;
			$posSum += @$datAref;
			if ($posMax < @$datAref) {
				$posMax = @$datAref;
				@datMPath=($Chr,$Strand);
			}
		}
	}
	$datAref=$ScaffAlign{$Scaff}{$datMPath[0]}{$datMPath[1]};
	if ($opt_v && $posMax != $posSum) {
		print '-'x75,"\n[!][DiffStrand]$Scaff:\t$posSum -> $posMax\n";
		ddx $ScaffAlign{$Scaff};
	}
	if ($posMax < $posSum * $scaffmainR) {
		delete $ScaffAlign{$Scaff};
		--$outScaff;
		print "Deleted !\n" if $opt_v;
		next;
	} elsif ($posMax>1) {
		#$datAref = [ sort {$a->[0] <=> $b->[0] } @{$datAref} ];
		#my @dat=@{$datAref};
		my @dat=sort {$a->[0] <=> $b->[0] } @{$datAref};
		my (%IndelsDat,@IndelsGroup);
		for (@dat) {
			my ($Spos,$Cpos)=@$_;
			push @{$IndelsDat{$Cpos-$Spos}},$_;
		}
		#ddx \%IndelsDat;
		my $majorIndelGroup2R=0.6;
		my $majorIndelGroup3R=0.85;
		@dat=sort {$#{$IndelsDat{$b}} <=> $#{$IndelsDat{$a}} } keys %IndelsDat;
		my ($Sum,@CountArr,$t)=(0);
		for (@dat) {
			$t=scalar @{$IndelsDat{$_}};
			push @CountArr,$t;
			$Sum += $t;
		}
		$t=0;
		my ($flag,$ss)=(0,0);
		for (@CountArr) {
			++$t;
			$ss += $_;
			if ($t == 2) {
				$flag=1 if $ss/$Sum >= $majorIndelGroup2R;
			} elsif ($t == 3) {
				$flag=1 if $ss/$Sum >= $majorIndelGroup3R;
			}
		}
		$flag=1 if $t==1;
		if ($flag==0) {
			delete $ScaffAlign{$Scaff};
			--$outScaff;
			print "Deleted,too !\n" if $opt_v;
			next;
		}
		#ddx $ScaffAlign{$Scaff};
		@dat=sort { ${$IndelsDat{$a}}[0][0] <=> ${$IndelsDat{$b}}[0][0] } keys %IndelsDat;
		$datAref=[$#dat];
		for (@dat) {
			push @$datAref,$_,$IndelsDat{$_};
		}
	} else {
		$datAref=[0,$$datAref[0][1]-$$datAref[0][0],$datAref];
	}
	$ScaffAlign{$Scaff}=[@datMPath,$datAref];
}
warn "[!]Scaffold Count: $inScaff -> $outScaff\n";
#ddx \%ScaffAlign;
#   scaffold95088 => [
#                      "Chr10",
#                      1,
#                      [
#                        2,
#                        12418282,
#                        [[1476, 12419758], [1690, 12419972]],
#                        12418248,
#                        [
#                          [2154, 12420402],
#                          [2175, 12420423],
#                          [3093, 12421341],
#                          [5973, 12424221],
#                          [6197, 12424445],
#                        ],
#                        12418246,
#                        [
#                          [6235, 12424481],
#                          [6255, 12424501],
#                          [7279, 12425525],
#                        ],
#                      ],
#                    ],
#   scaffold9564  => [
#                      "Chr10",
#                      -1,
#                      [1, 8711880, [[445, 8712325]], 8710850, [[997, 8711847]]]
#                    ],
#   scaffold92619 => ["Chr10", 1, [0, 17878303, [[2999, 17881302]] ] ],




__END__
./alignbympos.pl markerpos/m2sChr10.pos.f markerpos/m2cChr10.pos.f ../9311.Nzone ../chr.nfo t.out

1. change to contig.
2. indel only
