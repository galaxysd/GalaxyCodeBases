#!/bin/env perl
use strict;
use warnings;

my $path='E:\BGI\gaosj\dat\\';

unless (@ARGV){
	print "perl $0 <begin date(\"10/22/2009\")> [days for UD (2)]\n";
	exit;
}

my $date = shift;
my $dayUD = shift;

my $Days=5;
my $UpRatio=.6;
my $UpRate=.02;
my $BorderRatioU=1/4;	# 1/4
my $BorderRatioD=1/5;	# 1/5

$dayUD=2 unless $dayUD;
$dayUD=2 unless $dayUD == 1;
$date='10/22/2009' unless $date=~/^\d+\/\d+\/\d{4}$/;
print "Caltulating from [$date] to $Days days later with [$dayUD]...\n";
my (%Names,%Max,%Min,%Files,%Values,%P);

opendir INDIR,$path or warn "[!]Error opening [$path]: $!\n";
for my $file (readdir INDIR) {
	if ($file =~ /([^.]+)\.txt$/i) { $Files{$1}=$file;}
}
closedir INDIR;

sub Cal($$$) {
	my ($inref,$max,$min)=@_;
	my ($cur,$p,$ud)=@$inref;
	my ($av,$ap,$aud,$bv,$bp,$bud);
	if ($cur*(1+$UpRate) > $max) {
		$aud=$bud=-1;
		$ap=$bp=$p*0.5;
	} elsif ($cur*(1-$UpRate) < $min) {
		$aud=$bud=1;
		$ap=$bp=$p*0.5;
	} else {
		$aud=$ud;
		$bud=-$ud;
		$ap=$p*$UpRatio;
		$bp=$p*(1-$UpRatio);
	}

	$av=$cur*(1+$aud*$UpRate);
	$bv=$cur*(1+$bud*$UpRate);
	#my @t=([$av,$ap,$aud],[$bv,$bp,$bud]);
	return [[$av,$ap,$aud],[$bv,$bp,$bud]];
}

my $skip=0;
MAIN: for my $id (keys %Files) {
	open IN,'<',$path.$Files{$id} or die "[!]Error opening [$path$Files{$id}]: $!\n";
	my $line=<IN>;
	$line =~ /(\d+) (.+) \S+$/;
	my $ID=$1;
	$Names{$ID}=$2;
#print "$id\t$1\t",$Names{$ID},"\n";
	<IN>;
	my ($max0,$min0,$end0,$endpre,$endpt)=(-999999,999999,-1,-1);
	while (<IN>) {
		$endpt=$end0;
#      日期	    开盘	    最高	    最低	    收盘	    成交量	    成交额
		my ($dateitem,$begin,$max,$min,$end)=split /\t/;
		$max0=$max if $max0<$max;
		$min0=$min if $min0>$min;
		if ($date eq $dateitem) {
			$end0=$end;
			$endpre=($dayUD==2)?$endpt:$begin;
		}
	}
	close IN;
	$Max{$id}=$max0;
	$Min{$id}=$min0;
#print "$Max{$id},$Min{$id}\t$end0\n";
	unless ($end0) {
		#warn "[!]$date not found in $id !\n";
		++$skip;
		next MAIN;
		#$end0=0;
	}
#print "$Max{$id},$Min{$id}\t$end0\t";
	#my $Result=&Cal($max0,$min0,$end0);
	my $UDbegin=($end0>=$endpre)?1:-1;
	my @Days=([[$end0,1,$UDbegin]]);	# ($cur,$p,$ud)
	my $max=($max0-$end0)*$BorderRatioU+$end0;
	my $min=($min0-$end0)*$BorderRatioU+$end0;
	for (my $i=1;$i<$Days;$i++) {	# 0..4
		$Days[$i]=[];
		push @{$Days[$i]},@{&Cal($_,$max,$min)} for @{$Days[$i-1]};
	}
	my ($i,$resp)=(0);
	for my $day (@Days) {
		++$i;
		my ($ups,$upp,$downs,$downp)=(0,0,0,0);
		for (@$day) {
			my ($cur,$p,$ud)=@$_;
			if ($cur >= $end0) {
				$ups += $cur*$p;
				$upp += $p;
			} else {
				$downs += $cur*$p;
				$downp += $p;
			}
		}
		$resp += $upp;
	}
	#$Values{$ID}=$resv;
	$P{$ID}=$resp/$i;
#print "$resv\t$resp\n";# if $resv > 1;
}

warn "[!]$skip skipped !\n" if $skip;

my @KEYs=sort {$P{$b} <=> $P{$a}} keys %P;

for my $i (0..9) {
	my $key=$KEYs[$i];
	print "$key\t$Names{$key}:\t$P{$key}\n";
}
print "\n";
for my $i (-5..-1) {
	my $key=$KEYs[$i];
	print "$key\t$Names{$key}:\t$P{$key}\n";
}

my %t;
for (values %P) {
	++$t{$_};
}
for (sort {$b <=> $a} keys %t) {
	print "$_\t$t{$_}\n";
}
