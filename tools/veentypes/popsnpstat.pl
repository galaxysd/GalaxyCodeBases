#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:\BGI\toGit\perlib\etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
#use GalaxyXS::ChromByte 1.02;
#use DBI;
#
$main::VERSION=0.0.1;

our $opts='i:o:d:bv';
our ($opt_i, $opt_o, $opt_v, $opt_b, $opt_d);

our $help=<<EOH;
\t-i Group list (./veen.lst)
\t-o Output Stat (stat.txt)
\t-d Details dump to (details.lst)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='./veen.lst' if ! defined $opt_i;
$opt_o='stat.txt' if ! $opt_o;
$opt_d='details.lst' if ! $opt_d;

print STDERR "From [$opt_i] to [$opt_o] [$opt_d]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

#http://doc.bioperl.org/bioperl-live/Bio/Tools/IUPAC.html#BEGIN1
our %IUB = ( A => [qw(A)],
	     C => [qw(C)],
	     G => [qw(G)],
	     T => [qw(T)],
	     U => [qw(U)],
	     M => [qw(A C)],
	     R => [qw(A G)],
	     W => [qw(A T)],
	     S => [qw(C G)],
	     Y => [qw(C T)],
	     K => [qw(G T)],
	     V => [qw(A C G)],
	     H => [qw(A C T)],
	     D => [qw(A G T)],
	     B => [qw(C G T)],
	     X => [qw(G A T C)],
	     N => [qw(G A T C)]
	     );
our %REV_IUB = (A	=> 'A',
		T	=> 'T',
		C	=> 'C',
		G 	=> 'G',
		AC	=> 'M',
		AG	=> 'R',
		AT	=> 'W',
		CG	=> 'S',
		CT	=> 'Y',
		'GT'	=> 'K',
		ACG	=> 'V',
		ACT	=> 'H',
		AGT	=> 'D',
		CGT	=> 'B',
		ACGT=> 'N',
		N	=> 'N'
		);

my ($bit,$i,%Log2,%files,%bit2id)=(1,0);
open( L,'<',$opt_i) or die "[x]Error: $!\n";
while (<L>) {
	use integer;
	chomp;
	my ($id,$file)=split /\t/;
	$Log2{$bit}=$i;
	$bit2id{$bit}=$id;
	$files{$bit}=$file;
	warn "[!]$i: $bit,$id -> $file\n";
	++$i;
	$bit *= 2;
}

my (%ResultS,%ResultD,%Tables,%Bits);
my @BitsA=sort {$a <=> $b} keys %Log2;
--$bit;
for my $v (1..$bit) {
	use integer;
	for my $b (@BitsA) {
		my $n=$Log2{$b};
		my $t=$v & $b;
		next if $t == 0;
		next if $v >> $n == 0;	# if $v >> $n+1 == 0, then $n self will be ommited.
		my $w=($v << $i-$n) & $bit;
		next if $w != 0;
		$ResultS{$v}=0;
	}
	my @tt=split //,(sprintf "%b",$v);	# "%.${i}b"
	my ($ii,$tt,@Tab)=1;
	for (reverse @tt) {
		if ($_ == 1) {
			++$tt;
			push @Tab,$bit2id{$ii};
		}
		$ii *= 2;
	}
	$Bits{$v}=$tt;
	$Tables{$v}=join '',@Tab;
	warn "$v\t@tt\t$tt $Tables{$v}\n";
}

my %SNP;
for my $id (keys %files) {
	open IN,'<',$files{$id} or die "[x]Error: $!\n";
	++$i;
	while (<IN>) {
		chomp;
		my ($fh,%Bases,$Base);
		print STDERR ">> $_ .\b";
		open $fh,'<',$_ or die "[x]Error: $!\n";
		while (<$fh>) {
			my ($chr,$pos,$ref,$tail)=split /\t/;
			%Bases=();
			for (split / /,$tail) {
				next unless /[ACGTRYMKSWHBVDNX]/;
				for (@{$IUB{$_}}) {
					++$Bases{$_};
				}
			}
			$Base=$REV_IUB{join('',sort keys %Bases)};
			$SNP{$chr}{$pos}{$id}=$Base;
			#warn "$Base\t$tail";
		}
		close $fh;
		warn "<\n";
	}
}

#my (%ResultS,%ResultD);
for my $chr (keys %SNP) {
	for my $pos (keys %{$SNP{$chr}}) {
		my %Values;
		for my $bit (keys %{$SNP{$chr}{$pos}}) {
			$Values{$SNP{$chr}{$pos}{$bit}} += $bit;
		}
		for (keys %Values) {
			++$ResultS{$Values{$_}};
			push @{$ResultD{$Values{$_}}},[$chr,$pos,$_];
		}
		if ($opt_v or scalar (keys %Values) > 1) {
			print "[v]";
			print "$_,$SNP{$chr}{$pos}{$_} " for keys %{$SNP{$chr}{$pos}};
			print "-> ";
			print "$Values{$_},$_ " for keys %Values;
			print "\n";
		}
		delete $SNP{$chr}{$pos};
	}
}


open O,'>',$opt_o or die "[x]Error: $!\n";
for (sort { $Bits{$b} <=> $Bits{$a} } keys %ResultS) {
	print O "$Tables{$_}\t$ResultS{$_}\n";
}
close O;

open S,'>',$opt_d or die "[x]Error: $!\n";
for (sort { $Bits{$a} <=> $Bits{$b} } keys %ResultD) {
	print S "\n[$Tables{$_}]\n";
	my $array=$ResultD{$_};
	for (@$array) {
		print S join("\t",@$_),"\n";
	}
}
close S;
