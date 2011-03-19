#!/bin/env perl
use lib '/share/raid010/resequencing/resequencing/tmp/bin/annotation/glfsqlite';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.2;

our $opts='i:o:s:m:x:q:l:bv';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_m, $opt_x, $opt_q, $opt_l);

our $help=<<EOH;
\t-i PSNP list (./psnp.lst) for chrid.add_ref
\t-l Indel list (./indel.lst) in [SampleID\\tpath/to/indel-result.list]
\t-s GLF list (undef), will use \$1 of (([^/]+)/[^/]+$) for sample names
\t-o Output Prefix (./indel_f/f_)
\t-m minimal depth of the indel (3)
\t-x maximal depth of the indel (20)
\t-q the minimal quality of the indel (20)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='./indel_f/f_' if ! defined $opt_o;
$opt_i='./psnp.lst' if ! $opt_i;
$opt_l='./indel.lst' if ! $opt_l;
#$opt_s='./glf.lst' if ! $opt_s;
$opt_s='_UseAdd_ref_' if ! $opt_s;
$opt_m=3 if ! $opt_m;
$opt_x=20 if ! $opt_x;
$opt_q=20 if ! $opt_q;

my $dis = 4;

$opt_i =~ s/\/$//;
print STDERR "From [$opt_i][$opt_l] to [$opt_o], with [$opt_s]($opt_m,$opt_x)[$opt_q]/\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];

system('mkdir','-p',$opt_o);
system('rmdir',$opt_o) if $opt_o =~ m#/[\s.]*[^/\s.]+[^/]*$#;

my @Samples;
if ($opt_s eq '_UseAdd_ref_') {
	my $line;
	open L,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";
	chomp($line=<L>);
	close L;
	open L,'<',$line or die "[x]Error opening $line: $!\n";
	chomp($line=<L>);
	$line=(split /\t/,$line)[3];
	@Samples=split / /,$line;
	print STDERR '[!]Sample Order: ',(scalar @Samples),':[',join('] [',@Samples),']';

} else {
	open L,'<',$opt_s or die "[x]Error opening $opt_s: $!\n";
	print STDERR '[!]Sample Order: ';
	while (<L>) {
		m#([^/]+)/[^/]+$#;
		push @Samples,$1;
		print STDERR (scalar @Samples),':[',$1,"] ";
	}
}
close L;
print STDERR "\n";

my ($C_PSNP,$C_iSNP,%SNP)=(0,0);
print STDERR "[!]Parsing SNP ";
open P,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";
while (my $file=<P>) {
	chomp $file;
	open SNP,'<',$file or (warn "\n[!]Error opening $file: $!\n" and next);
	print STDERR ".\b";
	my ($chr,$pos,$ref,$tail,$i);
	while (<SNP>) {
		($chr,$pos,$ref,$tail)=split /\t/;
		++$C_PSNP;
		my @indSNP=split / /,$tail;	# /[ACGTRYMKSWHBVDNX-]/
		for my $s (@Samples) {
			$tail=shift @indSNP;	# ReCycle ...
			if ($tail ne '-') {
				#$SNP{$chr}{$_}{$s}='.' for ($pos-$dis .. $pos+$dis);
				$SNP{$chr}{$pos}{$s}='.';	# the above requires about 9 times memory ...
				++$C_iSNP;
			} else {
				next;
			}
		}
	}
	close SNP;
	print STDERR '-';
}
close P;
warn "\n[!]PSNP: $C_PSNP <= $C_iSNP iSNP in ",scalar @Samples," samples\n";

my $OutNameP;;
print STDERR "[!]Filtering indel:\n";
open P,'<',$opt_l or die "[x]Error opening $opt_l: $!\n";
while (<P>) {
	chomp;
	my ($s,$file)=split /\t/;
	open FH,'<',$file or (warn "\n[!]Error opening $file: $!\n" and next);
	$OutNameP=$opt_o.$s.'.indel-result.';
	open OUT,'>',$OutNameP.'filter' or warn "\n[!]Error creating ${OutNameP}filter: $!\n";
	open OUTC,'>',$OutNameP.'summary' or warn "\n[!]Error creating ${OutNameP}summary: $!\n";
	#print STDERR ".\b";
	my $snp_indel = 0;
	my $homo = 0;
	my $hete = 0;
	my $insertion = 0;
	my $deletion = 0;
	my $total = 0;
	my $sum = 0;

	INDEL: while (<FH>) {
		++$sum;
		my @indel=split /\t/;
		next if ($indel[6] < $opt_q);
		next if (($indel[7] < $opt_m) || ($indel[7] > $opt_x));
		for ($indel[1]-$dis .. $indel[1]+$dis) {
			if (exists $SNP{$indel[0]}{$_}{$s}){
				$snp_indel ++;
				next INDEL;
			}
		}
		$total ++;
		if ($indel[5] eq "homo"){
			$homo ++;
		}else{
			$hete ++;
		}
		if ($indel[2] =~ /^I/){
			$insertion ++;
		}else{
			$deletion ++;
		}
		print OUT $_;
	}
	warn "-> $s: $sum -> $total\n";
	next if ($sum == 0 or $total==0);
	print OUTC "total\thomo\thomo%\thete\thete%\tinsertion\tins%\tdeletion\tdel%\tINDEL_SNP\tINDEL_SNP%\n";
	printf OUTC ("%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\n", $total, $homo, $homo/$total*100, $hete, $hete/$total*100, $insertion, $insertion/$total*100, $deletion, $deletion/$total*100, $snp_indel, $snp_indel/$sum*100);
	print OUTC "
Raw Indel:\t$sum
SNP nearby:\t$snp_indel\t",100*$snp_indel/$sum," %

Filtered Indel:\t$total\t",100*$total/$sum," %
Home Indel:\t$homo\t",100*$homo/$total," %
Hete Indel:\t$hete\t",100*$hete/$total," %
Insertion:\t$insertion\t",100*$insertion/$total," %
Deletion:\t$deletion\t",100*$deletion/$total," %
";
	close OUTC;
	close OUT;
	close FH;
	#print STDERR '-';
}
warn "\n[!]Done !\n";






my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
68629431 iSNP => 8.6GB

find ./Indel/ -name indel-result.list|grep -P '\d'|perl -lane '$_ =~ /\/(\w+)\/indel-result\.list$/;print "$1\t$_"' > indel.lst.n

find ./indel_f/ -name *indel-result.filter|grep -P '\d'|perl -lane '$_ =~ /\/f_([^.]+)\.[^\/]+$/;print "$1\t$_"' > indel.lst.n
