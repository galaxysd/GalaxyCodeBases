#!/bin/env perl
use lib '/nas/RD_09C/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.3;

our $opts='i:o:s:f:d:bv';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_f, $opt_d);

our $help=<<EOH;
\t-i PSNP list (./psnp.lst) for chrid.individual.finalSNPs
\t-d Indel list (./indel.lst), in format [^SampleID\\tindel-result.filter\\n\$]
\t-f fabyChr path (./faByChr) for [chrid].fa
\t-s GLF list (./glf.lst), will use \$1 of (([^/]+)/[^/]+\$) for sample names
\t-o Output Prefix (./indGenomes/ig_)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='./indGenomes/ig_' if ! defined $opt_o;
$opt_i='./psnp.lst' if ! $opt_i;
$opt_s='./glf.lst' if ! $opt_s;
$opt_f='./faByChr' if ! $opt_f;
$opt_d='./indel.lst' if ! $opt_d;

$opt_i =~ s/\/$//;
$opt_f =~ s/\/$//;
print STDERR "From [$opt_i] to [$opt_o], with [$opt_s][$opt_d][$opt_f]/\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];

system('mkdir','-p',$opt_o);
system('rmdir',$opt_o) if $opt_o =~ m#/[\s.]*[^/\s.]+[^/]*$#;

my (@Samples,%SampleToIndelf);
open L,'<',$opt_s or die "[x]Error opening $opt_s: $!\n";
print STDERR "[!]Sample Order: ";
while (<L>) {
	m#([^/]+)/[^/]+$#;
	push @Samples,$1;
	print STDERR (scalar @Samples),':[',$1,"] ";
}
print STDERR "\n";

open L,'<',$opt_d or die "[x]Error opening $opt_d: $!\n";
print STDERR "[!]Indel files:\n";
while (<L>) {
	chomp;
	my ($sample,$file)=split /\t/;
	$SampleToIndelf{$sample}=$file;
	print STDERR "\t$sample -> $file\n";
}
print STDERR "\n";

sub GetNextInDel($$) {
	my ($Aref,$theChr)=@_;
	my ($FH,$curChr,$curSt,$curEd,$curInDe,$curIDseq)=@$Aref;
	my $filePos = tell $FH;	# maybe useful in future ?
	$curChr=undef;	# flag
	$curSt=$curEd=-1;
	while(<$FH>) {
		my ($chr,$pos,$indel,$bases,$strand,$homhet,undef,$depth,$t) = split /\t/;
		unless ($t) {   # if file not completed.
			print STDERR '\'';
			next;
		}
		next if $chr ne $theChr;
		$bases = lc $bases if $homhet eq 'hete';
		if ($indel =~ /^([ID])(\d+)$/) {
			$curSt=$pos;
			$curIDseq=$bases;
			if ($1 eq 'I') {
				$curInDe = $2;
				$curEd=$pos;
			} else {
				$curInDe = -$2;
				$curEd = $pos-1 + $2;
			}
		} else {print STDERR '|';next;}
		$curChr=$chr;
		@$Aref=($FH,$curChr,$curSt,$curEd,$curInDe,$curIDseq);
		if (defined $curChr) { return $filePos; }
		 else { return -1; }
	}
	return -2;	# Empty file
}

sub QueryInDel($$$$) {
	my ($Aref,$theChr,$pos,$theBase)=@_;
	my ($FH,$curChr,$curSt,$curEd,$curInDe,$curIDseq)=@$Aref;
	return $theBase if $pos < $curSt;	# Since we Recursion below, these two if is a must.
	if ($pos >= $curSt and $pos <= $curEd) {
#return '<D>' if $curInDe<0;
		return '' if $curInDe<0;
		# Now, must be Insertion:
#print STDERR 'X' if $curSt ne $curEd;	# just to recheck, safe to remove this line.
#return $theBase.'<I>'.$curIDseq.'</I>';
		return $theBase.$curIDseq;
	} else {
		my $r=&GetNextInDel($Aref,$theChr);
		if ($r < 0) {
			return $theBase;
		} else { return &QueryInDel($Aref,$theChr,$pos,$theBase); }	# should be safe to Recursion on normal indel results.
	}
	print STDERR 'x';	# should never be here.
	return '.';
}

print STDERR "[!]Parsing SNP ";
open P,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";
while (my $file=<P>) {
	chomp $file;
	open SNP,'<',$file or (warn "\n[!]Error opening $file: $!\n" and next);
	print STDERR ".\b";
	my (%SNP,$chr,$pos,$ref,$tail,$i);
	while (<SNP>) {
		chomp;
		($chr,$pos,$ref,$tail)=split /\t/;
		my @indSNP;
		for (split / /,$tail) {
			next unless /[ACGTRYMKSWHBVDNX-]/;
			s/-/n/;
			push @indSNP,$_;
		}
		$SNP{$pos}=[$ref,\@indSNP];
	}
	print STDERR "+\b";

	my (@FH,@IndelH);
	for (@Samples) {
		$file=$opt_o.$_.'.'.$chr.'.fa';
		my ($fh,$fhIDl);
		open $fh,'>',$file or die "[x]Error opening $file: $!\n";
		print $fh ">${_}---$chr\n";
		push @FH,$fh;
		open $fhIDl,'<',$SampleToIndelf{$_} or die "[x]Error opening $SampleToIndelf{$_}: $!\n";
		push @IndelH,[$fhIDl];
		&GetNextInDel($IndelH[-1],$chr) >=0 || die;
	}
	warn '[!]PSNP:[',1+$#{${$SNP{$pos}}[1]},'] != File:[',(scalar @FH),"].\n" if $#FH != $#{${$SNP{$pos}}[1]};

	$file=$opt_f.'/'.$chr.'.fa';
	open FA,'<',$file or (warn "\n[!]Error opening $file: $!\n" and next);
	$ref=<FA>;
	$ref =~ />(\S+)/;
	warn "[!]Different ChrID, [>$1] in [$file] !\n" if $1 ne $chr;
	$i=1;
	while ($ref=getc(FA)) {
		unless ($ref =~ /[ACGTRYMKSWHBVDNX]/i) {
			last if $ref eq '>';
			next;
		}
		unless ($i%80) {
			print $_ "\n" for @FH;
		}
		unless ($SNP{$i}) {	# Normal or InDel, which cannot start at '-' as not covered or SNP as filtering
			my $t=0;
			for my $fh (@FH) {
				my $str;
				#my ($FH,$curChr,$curSt,$curEd,$curInDe,$curIDseq)=@$Aref;
				if ($i < $IndelH[$t]->[2]) {
					$str = $ref;
				} elsif ($i <= $IndelH[$t]->[3]) {
					if ($IndelH[$t]->[4]<0) { $str = ''; }
					 else { $str = $ref . $IndelH[$t]->[5]; }
				} else { $str = QueryInDel($IndelH[$t],$chr,$i,$ref); }	# Well, the above is for speed
				print $fh $str;

				++$t;
			}
			#print $_ $ref for @FH;
		} else {	# SNP
			my ($refbase,$indSNPr)=@{$SNP{$i}};
			warn "[!]RefBase differ, SNP:[$refbase] ne FASTA:[$ref].\n" if $refbase ne uc($ref);
			my $t=0;
			for (@$indSNPr) {
				$tail=$FH[$t];
				print $tail $_;
				++$t;
			}
		}
		++$i;
	}
	print STDERR '-';
}

my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
