#!/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my $eachDepth = 10;
my $ReadLen = 50;
my $minRefLen = 50;
my $maxRefLen = 100;
my $HomSNPrate = 0.0005;
my $HetSNPrate = 0.0005;
my $TranStoV = 4;	# 转换比颠换，transitions，transversions
my $QUAL = 'e';

my $outCnt = 1;

open I,'<','hg19chr17.bed.frag' or die;
open Z,'<','zone.lst' or die;

sub getSNP($) {
	my $inbase = $_[0];
	my $outBase = $inbase;
	my $SorV = rand(1);
	if ($SorV > $TranStoV/(1+$TranStoV)) {	# transversions
		my $AB = rand(1);
		if ($AB > 0.5) {
			$outBase =~ tr/AGCTagct/CCAAccaa/;
		} else {
			$outBase =~ tr/AGCTagct/TTGGttgg/;
		}
	} else {	# transitions
		$outBase =~ tr/AGCTagct/GATCgatc/;
	}
	return $outBase;
}

sub revcom($) {
    my $str = $_[0];
    $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    my $rev = reverse $str;
    $rev    =~ tr/[](){}<>/][)(}{></;
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

open S,'>','snp.lst';
#open M,'>','Met.lst';

my $SNPremained = 0;
my %usedPos;
sub addSNP($$$$$) {
	my ($chr,$s,$seq,$rate,$type) = @_;
	my $len = length $seq;
	my $SNPcount = $len * $rate + $SNPremained;
	$SNPremained = $SNPcount - int($SNPcount);
	my $newseq = $seq;
	for ( 1 .. int($SNPcount) ) {
		my $Pos;
		do {
			$Pos = int(rand($len));
		} until ( ! exists $usedPos{$chr}{$s}{$Pos} );
		$usedPos{$chr}{$s}{$Pos} = $type;
		my $oldbase = substr $newseq,$Pos,1;
		my $newbase = getSNP($oldbase);
		$newseq = $seq;
		substr $newseq,$Pos,1,$newbase;
		print S join("\t",$type,$chr,$s,$s+$len-1,$s+$Pos,$oldbase,$newbase ,$Pos+1,$SNPcount,$seq,$newseq),"\n";
	}
	return $seq;
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
	die if $seqlen < $readlen;
	my $r1 = substr $seq,0,$readlen;
	my $r2s = substr $seq,-$readlen,$readlen;
	my $r2 = revcom($r2s);
#print "$seqlen $r1,$r2s\t$seq\t$r2\n";
	return ($r1,$r2);
}

sub realdosim($$$$$$$$$) {
	my ($homhet,$theseq,$thedepth,$fhref,$unchgRate,$chr,$s,$e,$MetStatRef) = @_;
	my @FH = @$fhref;
	my ($flag,$newseq,$str,$seqname,$r1,$r2);
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

		# Hom -> seq1
		# fq1 is Plus
		$str = "$chr:$s $e fq1_Waston $outputhomhet";
		($flag,$newseq) = simPlusMinus($theseq,$unchgRate);
		++$$MetStatRef{$flag};
#print "$flag $unchgRate $str, $seq1, $newseq\n";
		($r1,$r2) = getPE($newseq,$ReadLen);
		$seqname =  "\@$outCnt $flag $unchgRate $str";
		$seqname =~ s/ /#/g;
#print "$_ $seqname\n";
		print {$FH[0]} join("\n",$seqname.'/1',$r1,'+',$QUAL x $ReadLen),"\n";
		print {$FH[1]} join("\n",$seqname.'/2',$r2,'+',$QUAL x $ReadLen),"\n";
		++$outCnt;
		
		# fq1 in Minus
		$str = "$chr:$s $e fq1_Crick $outputhomhet";
		my $revtheseq = revcom($theseq);
		($flag,$newseq) = simPlusMinus($revtheseq,$unchgRate);
		++$$MetStatRef{$flag};
#print "$flag $unchgRate $str, $seq1, $newseq\n";
		($r1,$r2) = getPE($newseq,$ReadLen);
		$seqname =  "\@$outCnt $flag $unchgRate $str";
		$seqname =~ s/ /#/g;
		print {$FH[2]} join("\n",$seqname.'/1',$r1,'+',$QUAL x $ReadLen),"\n";
		print {$FH[3]} join("\n",$seqname.'/2',$r2,'+',$QUAL x $ReadLen),"\n";
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
while(<I>) {
	chomp;
	my ($chr,$s,$e,$len,$seq) = split /\t/;
	next if $len < $minRefLen or $len > $maxRefLen;
#print "$s\n";
	if (exists $MetUnchgRate{$s}) {
		@Paras = @{ $MetUnchgRate{$s} };
	}
	my $seq1 = addSNP($chr,$s,$seq,$HomSNPrate,'Hom');
	my $seq2 = addSNP($chr,$s,$seq1,$HetSNPrate,'Het');
#ddx \@Paras; die;
	dosim(\@fhC,$Paras[1],$chr,$s,$e,$seq1,$seq2);
	dosim(\@fhN,$Paras[2],$chr,$s,$e,$seq1,$seq2);
}
close S;
#close M;
close $_ for (@fhC,@fhN);

system("gzip -fd *.fq.gz");

__END__

