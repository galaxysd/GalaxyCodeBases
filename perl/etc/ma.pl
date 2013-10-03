#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 boss_atk boss_hp [my_hp] [my_atk]\n" if @ARGV < 2;
my ($batk,$bhp,$myhp,$myatk)=@ARGV;
my $tmp='';
$tmp .= ", MyHP:$myhp" if $myhp;
$tmp .= ", MyATK:$myatk" if $myatk;
print "Input: BossATK:$batk, BossHP:$bhp$tmp\n";

my $MaxHP = 200000;
my $MaxATK = 300000;

my $Times = int(1+ $MaxHP / $batk);
for my $times ( 1 .. $Times) {
	my $needATK = int(1+ $bhp / $times);
	my $minHP = $batk * $times;
	next if $needATK > $MaxATK;
	last if $myhp and 3*$myhp < $minHP;
	my $tmp='';
	if ($times>1 and $myhp and $myhp > $minHP) {
		my $maxBossATK = int($myhp/$times);
		$tmp .= "\tMaxBossATK: $maxBossATK";
	}
	if ($myatk) {
		my $maxBossHP = 2*$myatk*$times;
		$tmp .= "\tMaxBossHP: $maxBossHP (2x)";
	}
	print "Times: $times\tHP: $minHP\tAtk: ".int(1+$needATK/3),"(x3) ",int(1+$needATK/2),"(x2) $needATK(x1)$tmp\n";
}

print "$Times\n";
__END__
MaxHP=154060
MaxATK=167033
