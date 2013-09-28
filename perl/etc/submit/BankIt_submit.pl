#!/usr/bin/perl
use strict;
use warnings;

open I1, "<", "Mt.fasta";
open I2, "<", "Y.fasta";
open I3, "<", "PLP.fasta";
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
	/^>(\w{3})(Y\w)/;
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
	"Ppl"	=>	"Prionailurus planiceps",
	"Pte"	=>	"Pardofelis temminckii",
	"Pma"	=>	"Pardofelis marmorata",
	"Pti"	=>	"Panthera tigris",
	"Ppa"	=>	"Panthera pardus",
);

=pod
print O "TYPE: Pub
TITLE:
$CITATION
AUTHORS:
$AUTHORS
YEAR: 2014
STATUS: 1
||\n\n";

print O "TYPE: Source
NAME: $SOURCE{$_}
ORGANISM: $SOURCE{$_}
||\n\n" foreach sort keys %SOURCE;

print O "TYPE: Cont
NAME: $CONT_NAME
TEL: +86-10-62752307
EMAIL: Luo.shujin\@pku.edu.cn
LAB: Peking-Tsinghua Center for Life Sciences, College of Life Sciences
INST: Peking University
ADDR: No. 5 Yiheyuan Road, Haidian District, Beijing 100871, China
||\n\n";
=cut

my %desc = (
	PLP => 'X chromosome, PLP1, partial CDs',
	CytB => 'mitochondrial haplotype Cytochrome B, partial CDs',
	'16S' => 'mitochondrial haplotype 16s Ribosomal RNA, partial sequence',
	UTY11 => 'Y-chromosome haplotype UTY, partial sequence',
	DBY7 => 'Y-chromosome haplotype DBY7, partial sequence',
	SMCY3 => 'Y-chromosome haplotype SMCY3, partial sequence',
	SMCY7_STR_upstream => 'Y-chromosome haplotype SMCY7, partial sequence',
	SMCY7_STR_downstream => 'Y-chromosome haplotype SMCY7, partial sequence'
);

foreach my $a (sort keys %STS) {
	open O, ">", "submit_$a.fa";
	foreach my $b (sort keys %{$STS{$a}}) {
		foreach my $c (sort keys %{$STS{$a}{$b}}) {
			$STS{$a}{$b}{$c} =~ s/-|\?//g;
			next unless $STS{$a}{$b}{$c};
			next unless $SOURCE{$b};
			my $str = "$SOURCE{$b} $desc{$a}";
			if ($a eq "PLP") {
				print O ">$b$c [organism=$SOURCE{$b}] [chromosome=X] [gcode=1] $str\n$STS{$a}{$b}{$c}\n";
			} elsif ( $a eq "CytB" or $a eq "16S" ) {
			print O ">$a-$b$c [organism=$SOURCE{$b}] [location=mitochondrion] [mgcode=2] $str\n$STS{$a}{$b}{$c}\n";
			} elsif ( $a =~ /^(SMCY|DBY7|UTY11)/ ) {
				print O ">$a-$b$c [organism=$SOURCE{$b}] [chromosome=Y] [gcode=1] $str\n$STS{$a}{$b}{$c}\n";
			} else {
				print O ">$a-$b$c [organism=$SOURCE{$b}] [gcode=1]\n$STS{$a}{$b}{$c}\n";
				die;
			}
		}
	}
	close O;
}
#close O;

system("cat submit_*.fa > submitall.fa");
system('find submit*.fa|sed \'s/\.fa$//\'|xargs -n1 perl BankIt_Feature.pl');
print '-' x 75,"\n";
system('ls -l *.val');

print "To clean, run:[ rm submit* ]\n";
