#!/usr/bin/perl
use strict;
use warnings;

open I1, "<", "../1.find_need_to_submit/MT/MT_needto.fasta";
open I2, "<", "../1.find_need_to_submit/Y/Y_needto.fasta";
open I3, "<", "../1.find_need_to_submit/PLP/PLP_needto.fasta";
#open O, ">", "BankIt_submit.fasta";

my %STS;
while (<I1>) {
	/^>(\w{3})(MT\d+)/;
	my $s = <I1>;
	$STS{"16S"}{$1}{$2} = substr $s, 0, 367, "";
	$STS{ATP8}{$1}{$2} = substr $s, 0, 176, "";
	$STS{CytB}{$1}{$2} = substr $s, 0, 1249, "";
	warn $s unless $s eq "\r\n";
	delete $STS{ATP8};
}
close I1;
while (<I2>) {
	/^>(\w{3})(Y)/;
	my $s = <I2>;
	$STS{SMCY3}{$1}{$2} = substr $s, 0, 841, "";
	$STS{SMCY7_STR_upstream}{$1}{$2} = substr $s, 0, 224, "";
	$STS{SMCY7_STR_downstream}{$1}{$2} = substr $s, 0, 325, "";
	$STS{DBY7}{$1}{$2} = substr $s, 0, 280, "";
	$STS{UTY11}{$1}{$2} = substr $s, 0, 484, "";
	warn $s unless $s eq "\r\n";
}
close I2;
while (<I3>) {
	$/ = ">";
	my $s = <I3>;
	chomp $s;
	$s =~ s/\s//g;
	$/ = "\n";
	/^>?(\w{3})([\--Z]+)/;
	$STS{PLP}{$1}{$2} = $s;
}
close I3;

my %SOURCE = (
	"Pbe"	=>	"Prionailurus bengalensis",
	"Pvi"	=>	"Prionailurus viverrinus",
	"Ipl"	=>	"Prionailurus planiceps",
	"Pte"	=>	"Pardofelis temminckii",
	"Pma"	=>	"Pardofelis marmorata",
	"Pti"	=>	"Panthera tigris",
	"Ppa"	=>	"Panthera pardus",
);

my %desc = (
	PLP => 'X chromosome, PLP1, partial sequence',
	CytB => 'mitochondrial, cytochrome b complete CDS, tRNA-Thr complete sequence, tRNA-Pro partial sequence',
	'16S' => 'mitochondrial, 16S ribosomal RNA, partial sequence',
	UTY11 => 'Y chromosome, UTY, partial sequence',
	DBY7 => 'Y chromosome, DBY, partial sequence',
	SMCY3 => 'Y chromosome, SMCY, partial sequence',
	SMCY7_STR_upstream => 'Y chromosome, SMCY, partial sequence',
	SMCY7_STR_downstream => 'Y chromosome, SMCY, partial sequence'
);

foreach my $a (sort keys %STS) {
	open O, ">", "submit_$a.fa";
	foreach my $b (sort keys %{$STS{$a}}) {
		foreach my $c (sort keys %{$STS{$a}{$b}}) {
			$STS{$a}{$b}{$c} =~ s/-|\?//g;
			my $l = length $STS{$a}{$b}{$c};
			next if $l < 200;
			next unless $SOURCE{$b};
			my $str = "$SOURCE{$b}, $desc{$a}";
			if ($a eq "PLP") {
				print O ">$b$c [organism=$SOURCE{$b}] [chromosome=X] [gcode=1] $str\n$STS{$a}{$b}{$c}\n";
			} elsif ($a eq "CytB") {
				if ($l > 1000) {
					print O ">$a-$b$c [organism=$SOURCE{$b}] [location=mitochondrion] [mgcode=2] $str\n$STS{$a}{$b}{$c}\n";
				} else {
					print O ">$a-$b$c [organism=$SOURCE{$b}] [location=mitochondrion] [mgcode=2]";
					print O "$SOURCE{$b}, mitochondrial, cytochrome b, partical CDS\n$STS{$a}{$b}{$c}\n";
				}
			} elsif ($a eq "16S"){
				print O ">$a-$b$c [organism=$SOURCE{$b}] [location=mitochondrion] $str\n$STS{$a}{$b}{$c}\n";
			} elsif ($a =~ /^(SMCY|DBY7|UTY11)/) {
				print O ">$a-$b$c [organism=$SOURCE{$b}] [chromosome=Y] [gcode=1] $str\n$STS{$a}{$b}{$c}\n";
			} else {
				print O ">$a-$b$c [organism=$SOURCE{$b}] [gcode=1]\n$STS{$a}{$b}{$c}\n";
				warn "Bad gene name!";
			}
		}
	}
	close O;
}
#close O;
