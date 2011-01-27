#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
#use GalaxyXS::ChromByte;
use Data::Dump qw(ddx);

my $scaffKnot1R=0.1;
my $maxKdiffto1=0.05;
my $maxbRchg=0.001;
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
		my $maxKnot1cnt=int($scaffKnot1R * @dat);
		++$maxKnot1cnt unless $maxKnot1cnt;	# at least 1.
		my ($item,$Knot1cnt,%newDat)=(1,0);
		my ($Spos,$Cpos)=@{$dat[0]};
		push @{$newDat{$item}},[$Spos,$Cpos,$datMPath[1],$Cpos-$Spos];
		for (my $i=1;$i<=$#dat;$i++) {
			my ($S0pos,$C0pos)=@{$dat[$i-1]};
			my ($S1pos,$C1pos)=@{$dat[$i]};
			my $kp=($C1pos-$C0pos)/($S1pos-$S0pos);
			#my $bp=($S1pos*$C0pos-$S0pos*$C1pos)/($S1pos-$S0pos);
			my $bp=$C0pos-$S0pos;
			my $b=$C1pos-$S1pos; #my $k=$kp;
			if ( abs(abs($kp)-1) > $maxKdiffto1 ) {
				++$Knot1cnt;
push @{$newDat{$item}},[-$S1pos,-$C1pos,$kp,$bp,$b];
				next;	# skip those with too much indels
			} elsif ( abs(($b-$bp)/$bp) < $maxbRchg ) {
				push @{$newDat{$item}},[$S1pos,$C1pos,$kp-1*$datMPath[1],$bp,$b];
			} else {
				++$item;
				push @{$newDat{$item}},[$S1pos,$C1pos,$kp-1*$datMPath[1],$bp,$b];
			}
		}
=pod
		my ($S0pos,$C0pos)=@{shift @dat};
		my ($S1pos,$C1pos)=@{shift @dat};
	#warn $Scaff,"\n";
	#ddx $ScaffAlign{$Scaff} if $S1pos==$S0pos;
		my $k=($C1pos-$C0pos)/($S1pos-$S0pos);
		my $b=($S1pos*$C0pos-$S0pos*$C1pos)/($S1pos-$S0pos);
		my ($i,$Knot1cnt,%newDat)=(1,0);
		$valuecount=1;
		my $datCount=2;
		$newDat{$i}=[ [$S0pos,$C0pos],[$S1pos,$C1pos,$k-1,$b] ];
		for (@dat) {
			my ($Spos,$Cpos)=@$_;
			my $k1=($C1pos-$Cpos)/($S1pos-$Spos);
			my $b1=($S1pos*$Cpos-$Spos*$C1pos)/($S1pos-$Spos);
			++$datCount;
			if ( abs($k1-1) > $maxKdiffto1 ) {
				++$Knot1cnt;
				next;	# skip those with too much indels
			}
			if ( abs(($b1-$b)/$b) > $maxbRchg ) {
				++$i;
				$valuecount=0;
				$k=$k1; $b=$b1;
				$C1pos=$Cpos; $S1pos=$Spos;
			} else {
				$k=&addavg($k,$k1); $b=&addavg($b,$b1);
				$C1pos=&addavg($C1pos,$Cpos); $S1pos=&addavg($S1pos,$Spos);
			}
			push @{$newDat{$i}},[$Spos,$Cpos,$k-1,$b,$k1-1,$b1];
			++$valuecount;
		}
=cut
		if ($Knot1cnt > $maxKnot1cnt) {
			if ($opt_v) {
				print '-'x75,"\n[!][MoreBadK]Deleted: [$Scaff]<@datMPath> ",scalar @dat," $Knot1cnt > $maxKnot1cnt\t";
				#ddx $datAref;
				ddx \@dat;
				ddx \%newDat;
				print '-'x75,"\n";
			}
			delete $ScaffAlign{$Scaff};
			--$outScaff;
			next;	# skip this
		}
		print "[!]<@datMPath> ";
		ddx \%newDat;
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
