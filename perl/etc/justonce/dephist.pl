#!/bin/env perl
use strict;
use warnings;

my %dat=(
'Drone1' => 'd1.depth.xz',
'Scouts1' => 's1.depth.xz',
'Scouts2' => 's2.depth.xz',
'Recruit2' => 'r2.depth.xz',
'Recruit3' => 'r3.depth.xz',
);
my @id=sort keys %dat;

my $EffLen=199723619;
my $Ncnt= 219629612 - $EffLen;

sub openfile($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.xz$/) {
	    open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
    return $infile;
}

my (%Hist,%Stack);
sub dohist($$) {
	my ($id,$value)=@_;
	$value=0 if $value == 65535;
	if ($value <= 150) {
		++$Hist{$value}{$id};	# 0 .. 100
	} elsif ($value <= 1500) {
		++$Hist{200 + int(($value-100)/10)}{$id};	# 150~159->205, 1490~1499->339, 1500->340
	} else { ++$Hist{999}{$id}; }	# >1000
}
sub getV($){
	my $v=$_[0];
	if ($v<=150) {
		return $v;
	} elsif ($v<=340) {
		return '~'.(109 + ($v-200)*10);
	} else {
		return '>1500';
	}
}

for my $f (@id) {
	my $fh=openfile($dat{$f});
	$Hist{0}{$f} = -$Ncnt;
	while (<$fh>) {
		next if /^>/;
		chomp;
		my @d = split /\s+/;
		dohist($f,$_) for @d;
	}
}

my %Sum;
for my $v (sort {$b<=>$a} keys %Hist) {
	for my $id (@id) {
		unless (exists $Hist{$v}{$id}) {
			$Hist{$v}{$id}=0;
		}
		$Sum{$id} += $Hist{$v}{$id};
		$Stack{$v}{$id} = $Sum{$id};
	}
}
open O,'>',"depth.hist";
print O join("\t",'Depth',@id),"\n";
for my $v (sort {$a<=>$b} keys %Hist) {
	print O getV($v),"\t";
	for my $id (@id) {
		if (exists $Hist{$v}{$id}) {
			print O $Hist{$v}{$id},', ',$Stack{$v}{$id},',',int(10000*$Stack{$v}{$id}/$EffLen)/10000,"\t";
		}
	}
	print O "\n";
}
close O;
