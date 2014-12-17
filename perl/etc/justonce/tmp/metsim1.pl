#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my $eachDepth = 10;
my $ReadLen = 50;
my $QUAL = 'e';
my ($minRefLen,$maxRefLen) = (40,220);

my $outCnt = 1;

open I,'<','hg19chr18.bed.frag.trim' or die;
open Z,'<','zone.lst' or die;
open S,'<','snp.lst' or die;

sub revcom($) {
	my $str = $_[0];
	$str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	my $rev = reverse $str;
	$rev =~ tr/[](){}<>/][)(}{></;
	return $rev;
}

my %MetUnchgRate;
while (<Z>) {
	chomp;
	my (undef,$s,$e,$ra,$rb) = split /\t/;
	$MetUnchgRate{$s}=[$e,$ra,$rb];
}
close Z;
ddx \%MetUnchgRate;

my %SNPList;
while (<S>) {
	chomp;
	my ($HomHet,undef,$s,undef,undef,undef,undef,undef,undef,undef,$newseq) = split /\t/;
	$SNPList{"$HomHet\t$s"} = $newseq;
}
close S;

my $SNPremained = 0;
my %usedPos;
sub addSNP($$$$) {
	my ($chr,$s,$seq,$type) = @_;
	my $key = "$type\t$s";
	my $newseq = $seq;
	$newseq = $SNPList{$key} if exists $SNPList{$key};
	return $newseq;
}

sub getCpG($) {
	my ($seq) = @_;
	my @ret;
	while ($seq =~ /CG/gi) {
#print join("\t",'CG',$-[0], $+[0],"$` ($&) $'"),"\n";
		push @ret, $-[0];
	}
	return @ret;
}
sub getC($) {	# C but not CpG
	my ($seq) = @_;
	my @ret;
	while ($seq =~ /C(?!G)/gi) {
#print join("\t",'C',$-[0], $+[0],"$` ($&) $'"),"\n";
		push @ret, $-[0];
	}
	return @ret;
}

sub doCoin($) {
	my $rate = $_[0];
	my $coin = rand(1);
	if ($coin > $rate) {
		return 0;
	} else {
		return 1;
	}
}

sub simPlusMinus($$) {	#甲基化只处理fq1的单链
	my ($seq,$unchgRate) = @_;
	my @CGpos = getCpG($seq);
	my @Cpos = getC($seq);
	my $newseq = $seq;
	for (@Cpos) {
		substr $newseq,$_,1,'T';
	}
	my $flag = 'NoCpG';
	$flag = 'CpGisC' if @CGpos;
	unless (doCoin($unchgRate)) {
		for (@CGpos) {
			substr $newseq,$_,1,'T';
		}
		$flag = 'CpGtoT';
	}
	return ($flag,$newseq);
}
sub getPE ($$) {
	my ($seq,$readlen) = @_;
	my $seqlen = length $seq;
	#die if $seqlen < $readlen;
	if ($seqlen < $readlen) {
		$readlen = $seqlen;
	}
	my $r1 = substr $seq,0,$readlen;
	my $r2s = substr $seq,-$readlen,$readlen;
	my $r2 = revcom($r2s);
#print "$seqlen $r1,$r2s\t$seq\t$r2\n";
	return ($r1,$r2,$readlen);
}

sub realdosim($$$$$$$$$) {
	my ($homhet,$theseq,$thedepth,$fhref,$unchgRate,$chr,$s,$e,$MetStatRef) = @_;
	my @FH = @$fhref;
	my ($flag,$newseq,$str,$seqname,$r1,$r2,$realReadlen);
	for ( 1 .. $thedepth ) {
		my $outputhomhet;
		if (defined $usedPos{$chr}{$s}) {
			my @t = values %{$usedPos{$chr}{$s}};
			my %t; @t{@t} = 0 .. $#t;
			if (exists $t{'Het'}) {
				$outputhomhet = 'Het';
			} else {
				$outputhomhet = 'Hom';
			}
		} else {
			$outputhomhet = 'NoSNP';
		}
#print STDERR "$outputhomhet $unchgRate\t$chr,$s,$e\n";

		# Hom -> seq1
		# fq1 is Plus
		$str = "$chr:$s $e Waston $outputhomhet";
		($flag,$newseq) = simPlusMinus($theseq,$unchgRate);
		++$$MetStatRef{$flag};
#print "$flag $unchgRate $str, $seq1, $newseq\n";
		($r1,$r2,$realReadlen) = getPE($newseq,$ReadLen);
		$seqname =  "\@$outCnt $flag $unchgRate $str";
		$seqname =~ s/ /#/g;
#print "$_ $seqname\n";
		print {$FH[0]} join("\n",$seqname.'/1',$r1,'+',$QUAL x $realReadlen),"\n";
		print {$FH[1]} join("\n",$seqname.'/2',$r2,'+',$QUAL x $realReadlen),"\n";
		++$outCnt;
		
		# fq1 in Minus
		$str = "$chr:$s $e Crick $outputhomhet";
		my $revtheseq = revcom($theseq);
		($flag,$newseq) = simPlusMinus($revtheseq,$unchgRate);
		++$$MetStatRef{$flag};
#print "$flag $unchgRate $str, $seq1, $newseq\n";
		($r1,$r2,$realReadlen) = getPE($newseq,$ReadLen);
		$seqname =  "\@$outCnt $flag $unchgRate $str";
		$seqname =~ s/ /#/g;
		print {$FH[2]} join("\n",$seqname.'/1',$r1,'+',$QUAL x $realReadlen),"\n";
		print {$FH[3]} join("\n",$seqname.'/2',$r2,'+',$QUAL x $realReadlen),"\n";
		++$outCnt;
	}
}
sub dosim($$$$$$$) {
	my ($fhref,$unchgRate,$chr,$s,$e,$seq1,$seq2) = @_;
	my $depth1 = int(0.9 + $eachDepth/2);
	my $depth2 = $eachDepth - $depth1;
	my %MetStat = (
		'NoCpG' => 0, 'CpGisC' => 0, 'CpGtoT' => 0
	);
	unless (defined $unchgRate) {
		$unchgRate = 1;
		print STDERR '+';
	}
	realdosim('Hom',$seq1,$depth1,$fhref,$unchgRate,$chr,$s,$e,\%MetStat);
	realdosim('Het',$seq2,$depth2,$fhref,$unchgRate,$chr,$s,$e,\%MetStat) if $depth2>0;
	#print M join("\t",$chr,$s,$e,$e-$s+1,$MetStat{'CpGtoT'}),"\n";
}

my (@fhC,@fhN);
open $fhC[0],'|-',"gzip -9c > Cwaston.1.fq.gz" or die;
open $fhC[1],'|-',"gzip -9c > Cwaston.2.fq.gz" or die;
open $fhC[2],'|-',"gzip -9c > Ccrick.1.fq.gz" or die;
open $fhC[3],'|-',"gzip -9c > Ccrick.2.fq.gz" or die;
open $fhN[0],'|-',"gzip -9c > Nwaston.1.fq.gz" or die;
open $fhN[1],'|-',"gzip -9c > Nwaston.2.fq.gz" or die;
open $fhN[2],'|-',"gzip -9c > Ncrick.1.fq.gz" or die;
open $fhN[3],'|-',"gzip -9c > Ncrick.2.fq.gz" or die;
my @Paras;
print STDERR "[!]Run `metsim0.pl` to refresh zone.lst if see any cross below:\n";
while(<I>) {
	chomp;
	my ($chr,$s,$e,$len,$seq) = split /\t/;
	next if $len < $minRefLen or $len > $maxRefLen;
#print "$s\n";
	if (exists $MetUnchgRate{$s}) {
		@Paras = @{ $MetUnchgRate{$s} };
		ddx \@Paras;
	}
	my $seq1 = addSNP($chr,$s,$seq,'Hom');
	my $seq2 = addSNP($chr,$s,$seq1,'Het');
#ddx \@Paras; die;
my ($a,$b)=($seq1,$seq2);
	dosim(\@fhN,$Paras[1],$chr,$s,$e,$seq1,$seq2);
die if $a ne $seq1;
die if $b ne $seq2;
	dosim(\@fhC,$Paras[2],$chr,$s,$e,$seq1,$seq2);
die if $a ne $seq1;
die if $b ne $seq2;
}
close S;
#close M;
close $_ for (@fhC,@fhN);
warn "\n[!]Done !\n";
system("gzip -fd *.fq.gz");

__END__

