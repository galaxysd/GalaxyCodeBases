#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Vcf;
#use GTF;
use Data::Dump qw(ddx);
use Galaxy::IO::FASTAQ;
use Galaxy::SeqTools qw(translate revcom);

die "Usage: $0 <vcf.gz with tabix indexed> <out> <regions> (eg.: scaffold75:924209-5441687)\n" if @ARGV<2;
my $fafs='/share/users/huxs/work/tiger/paper/parents75n1458.fa';
my $gtfs = '/share/users/huxs/work/tiger/paper/parents75n1458.gtf';
my $vcfs = shift;
my $outfs = shift;
my $regions = shift;

warn "From [$vcfs][$gtfs][$fafs] to [$outfs]\n";
warn "Regions:[$regions]\n" if $regions;

my $code=<<Ecode;
    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ---M---------------M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
Ecode
my @t=split /\n/,$code;
my (%t,%gen_code,%start_codes);
map {s/\s//g;@_=split /=/;$t{$_[0]}=[split //,$_[1]]} @t;
#ddx \%t;
for (@{$t{AAs}}) {
#	my $aa=$_;
#	my @bases;
#	map {my $t=shift @{$t{$_}};push @bases,$t} qw/Base1 Base2 Base3/;
#	my $base=join('',@bases);
	my $base=(shift @{$t{Base1}}).(shift @{$t{Base2}}).(shift @{$t{Base3}});
#print "$base";
	my $start=shift @{$t{Starts}};
	++$start_codes{$base} if $start eq 'M';
	$gen_code{$base}=$_;
#print " -> [$_]\n";
}
#ddx [\@t,\%gen_code,\%start_codes];
print "Codon Table:\n$code\n";

my (%GeneDat,%Annoted);
open GTF,'<',$gtfs or dir $!;
while(my $in_line = <GTF>){
	next if ($in_line =~ /^\s*\#/);
	chomp $in_line;
	#remove leading whitespace and get comments
	my $comments = "";	  
	while($in_line =~ /^(.*)\#(.*)$/){
	$in_line = $1;
	$comments = $2.$comments;
	}
	if($in_line =~ /^\s+(.*)$/){
	$in_line = $1;
	}
	my @data = split /\s+/,$in_line;
	#verify line is correct length
	next if ($#data < 8);
	my %feature = map {s/\s*;$//; s/\"//g; $_;} splice @data,8;
#ddx [\@data,\%feature];
	if ($data[2] eq 'transcript') {
		die "[$feature{gene_type}]" if $feature{gene_type} ne 'protein_coding';
		$GeneDat{$data[0]}{$feature{gene_id}}=[$feature{gene_name},$data[6],[],$data[3]];
	} elsif ($data[2] eq 'CDS' or $data[2] eq 'stop_codon') {
		push @{$GeneDat{$data[0]}{$feature{gene_id}}->[2]},[$data[3],$data[4]];
		push @{$Annoted{$data[0]}},[$feature{gene_id},$data[3],$data[4]];
	}
}
close GTF;
for my $chr (keys %GeneDat) {
	$Annoted{$chr} = [ sort { $a->[1] <=> $b->[1] } @{$Annoted{$chr}} ];
	for my $id (keys %{$GeneDat{$chr}}) {
		my ($gene,$strand,$cdsA) = @{$GeneDat{$chr}{$id}};
		if ($strand eq '+') {
			$cdsA = [ sort {$a->[0] <=> $b->[0]} @$cdsA ];
		} elsif ($strand eq '-') {
			$cdsA = [ sort {$b->[0] <=> $a->[0]} @$cdsA ];
		} else { die; }
		$GeneDat{$chr}{$id}->[2] = $cdsA;
	}
}
warn "GTF loaded.\n";
#ddx \%GeneDat;
#ddx \%Annoted;

my $fafh;
my %RefSeq;
open $fafh,'<',$fafs or die "$!";
my @aux = undef;
my ($name, $comment, $seq, $ret);
my ($n, $len) = (0, 0);
while ( $ret = &readfq($fafh, \@aux) ) {
	($name, $comment, $seq) = @$ret;
	++$n;
	$len = length($seq);
	print "$n: ",join("\t", $name, $comment, $len), "\n";
	$RefSeq{$name} = $seq;
}
close $fafh;

my (%cDNA,%Protein);
open OUTCDNA,'>',"$outfs.cDNA.fa" or die $!;
open OUTPT,'>',"$outfs.protein.fa" or die $!;
for my $chr (sort keys %GeneDat) {
	for my $id (sort { $GeneDat{$chr}{$a}->[3] <=> $GeneDat{$chr}{$b}->[3] } keys %{$GeneDat{$chr}}) {
		my ($gene,$strand,$cdsA,$st) = @{$GeneDat{$chr}{$id}};
		print OUTCDNA ">${gene}_$id $chr:$strand:$st";
		print OUTPT ">${gene}_$id $chr:$strand:$st";
		my ($len,$seq,$tmpseq,$AA)=(0);
		for (@$cdsA) {
			my ($s,$e) = @$_;
			print OUTCDNA ",$s-$e";
			$len += $e-$s+1;
			$tmpseq = substr $RefSeq{$chr},$s-1,$e-$s+1;
			$tmpseq = revcom($tmpseq) if $strand eq '-';
			$seq .= $tmpseq."\n";
		}
		print OUTCDNA " =$len ",int($len*10/3)/10,"\n$seq\n";
		$seq =~ s/\s//g;
		$AA = translate($seq,\%gen_code);
		print OUTPT " $len ",int($len*10/3)/10,"\n$AA\n";
		$cDNA{$id} = $seq;
		$Protein{$id} = $AA;
	}
}
close OUTCDNA;
close OUTPT;

my $vcf = Vcf->new(file=>$vcfs,region=>$regions);
$vcf->parse_header();
my (@samples) = $vcf->get_samples();
ddx \@samples;

my (%mutGenes,%mutSNPs);
while (my $x=$vcf->next_data_hash()) {
	die unless exists $Annoted{$$x{CHROM}};
	#ddx $x if exists $$x{INFO}{INDEL};
	next if $$x{QUAL} < 20;
	my (%GTok,%GTcnt,%GTp1,%GTp2);
	my %GTs = %{$$x{gtypes}};
	my $sampleCNT=0;
	for my $sample (keys %GTs) {
		next if $GTs{$sample}{DP}<=0 or $GTs{$sample}{GQ}<20;
		++$sampleCNT;
		$GTok{$sample} = $GTs{$sample}{GT};
		++$GTcnt{ $GTs{$sample}{GT} };
		++$GTp1{ $GTs{$sample}{GT} } if $sample =~ /^B(HX|XH)01/;	# bb
		++$GTp2{ $GTs{$sample}{GT} } if $sample =~ /^JHH001/;	# Bb
	}
	next unless keys %GTok;
	next if (keys %GTcnt) != 2;
	next unless keys(%GTp1)==1 and keys(%GTp2)==1;
	my $GT1 = (keys %GTp1)[0];
	my $GT2 = (keys %GTp2)[0];
	my ($a1,$a2,$a3) = $vcf->split_gt($GT1);
	if ($a1 eq $a2) {
		$GT1 = $a1;	# b
	} else { next; }
	($a1,$a2,$a3) = $vcf->split_gt($GT2);
	if ($a1 ne $a2) {
		if ($a1 eq $GT1) {
			$GT2 = $a2;	# B
		} elsif ($a2 eq $GT1) {
			$GT2 = $a1;
		} else { next; }
	} else { next; }
#print "B:$GT2, b:$GT1 [$sampleCNT]\n";
	my @Bases = ($$x{REF},@{$$x{ALT}});
	my $GT2Base = $Bases[$GT2];
	my $GT1Base = $Bases[$GT1];
	my ($theGID,$theS,$theE) = @{ChechRange( $$x{CHROM},$$x{POS} )};
	next unless $theGID;
	die if exists $$x{INFO}{INDEL};	# No need to do INDEL as there is none.

	unless (exists $mutGenes{$$x{CHROM}}{$theGID}) {
#		if ($GeneDat{$$x{CHROM}}{$theGID}->[1] eq '+') {
#			$mutGenes{$$x{CHROM}}{$theGID} = [$cDNA{$theGID},$cDNA{$theGID}];
#		} elsif ($GeneDat{$$x{CHROM}}{$theGID}->[1] eq '-') {
#			my $tmpstr = revcom( $cDNA{$theGID} );
#			$mutGenes{$$x{CHROM}}{$theGID} = [$tmpstr,$tmpstr];
#		} else {die;}
		$mutGenes{$$x{CHROM}}{$theGID} = [$cDNA{$theGID},$cDNA{$theGID}];
	}
	push @{ $mutSNPs{$$x{CHROM}}{$theGID} },[$$x{POS}, $$x{REF}, $GT2Base, $GT1Base];

	print join(",",$$x{CHROM},$$x{POS},$GeneDat{$$x{CHROM}}{$theGID}->[0]." \tB:$GT2,$GT2Base b:$GT1,$GT1Base [$sampleCNT] ",
		$theGID,$GeneDat{$$x{CHROM}}{$theGID}->[1],$theS,$theE),"\n";
#	for my $gt (keys %GTs) {
#		my ($a1,$a2,$a3) = $vcf->split_gt($gt);
#		if ($a3 or ($a1 != $a2)) {
#			$flag = 1;
#		}
#	}
	#ddx \%GTs;
	#ddx $x;
# vcf2cds.pl:189: {
#   ALT    => ["A"],
#   CHROM  => "scaffold75",
#   FILTER => ["."],
#   FORMAT => ["GT", "PL", "DP", "SP", "GQ"],
#   gtypes => {
#               "BHX011.bam" => { DP => 19, GQ => 61, GT => "0/0", PL => "0,57,255", SP => 0 },
#               "BHX019.bam" => { DP => 26, GQ => 82, GT => "0/0", PL => "0,78,255", SP => 0 },
#               "JHH001.bam" => { DP => 28, GQ => 99, GT => "0/1", PL => "244,0,255", SP => 9 },
#             },
#   ID     => ".",
#   INFO   => {
#               AC1 => 1,
#               AF1 => 0.1667,
#               DP  => 74,
#               DP4 => "34,26,4,9",
#               FQ  => 999,
#               MQ  => 59,
#               PV4 => "0.13,1,0.12,1",
#               VDB => 0.0365,
#             },
#   POS    => 3572768,
#   QUAL   => 999,
#   REF    => "G",
# }
	#print $vcf->format_line($x);
}
#ddx $vcf;
$vcf->close;

open OUTDNA,'>',"$outfs.candidateDNA.fa" or die $!;
open OUTAA,'>',"$outfs.candidateAA.fa" or die $!;
#ddx \%mutSNPs;
for my $chr (keys %mutSNPs) {
	for my $gid (sort keys %{$mutSNPs{$chr}}) {
		for ( @{$mutSNPs{$chr}{$gid}} ) {
			my ($pos, $ref, $GT2Base, $GT1base) = @$_;
			$mutGenes{$chr}{$gid}->[0] = mutpoint($chr,$pos,$ref,$GT2Base,$gid,$mutGenes{$chr}{$gid}->[0]);
			$mutGenes{$chr}{$gid}->[1] = mutpoint($chr,$pos,$ref,$GT1base,$gid,$mutGenes{$chr}{$gid}->[1]);
		}
		my $id = $gid.$GeneDat{$chr}{$gid}->[1].$GeneDat{$chr}{$gid}->[0];
		my @diff = cmpstr($mutGenes{$chr}{$gid}->[0],$mutGenes{$chr}{$gid}->[1]);
#ddx \@diff;
		my $AA0 = translate($mutGenes{$chr}{$gid}->[0],\%gen_code);
		my $AA1 = translate($mutGenes{$chr}{$gid}->[1],\%gen_code);
		my @diffAA = cmpstr($AA0,$AA1);

		for (@diff,@diffAA) {
			$_ = join('',@$_[1,0,2]);
		}
		print OUTDNA ">${id}_B $chr\n",$mutGenes{$chr}{$gid}->[0],"\n",
			">${id}_b $chr ",join(',',scalar(@diff),@diff),"\n",$mutGenes{$chr}{$gid}->[1],"\n\n";
		print OUTAA ">${id}_B $chr\n",$AA0,"\n",
			">${id}_b $chr ",join(',',scalar(@diffAA),@diffAA),"\n",$AA1,"\n\n";
	}
}
close OUTDNA;
close OUTAA;

sub ChechRange($$) {
	my ($chr,$pos) = @_;
	my (%gidcnt,%gidat);
	if (exists $Annoted{$chr}) {
		my $datA = $Annoted{$chr};
		if ($pos > $$datA[-1][2] or $pos < $$datA[0][1]) {
			return [];
		}
		for (@{$datA}) {
			my ($gid,$s,$t) = @$_;
			next if $pos > $t;
			last if $pos < $s;
			++$gidcnt{$gid};
			push @{$gidat{$gid}},[$s,$t];
		}
		if ( (keys %gidcnt) == 1 ) {
			my $id=(keys %gidcnt)[0];
			return [$id,@{$gidat{$id}->[0]}];
		}
		die "$chr,$pos" if (keys %gidcnt) > 1;
	}
	return [];
	#push @{$Annoted{$data[0]}},[$feature{gene_id},$data[3],$data[4]];
}

sub mutpoint($$$$$$) {
	my ($chr,$pos,$ref,$gt,$gid,$seq) = @_;
	# $GeneDat{$chr}{$gene_id}=[$gene_name,$strand,[[s1,e1],[s2,e2]],$start];
	# $cDNA{$gid} = $seq;
	# $Protein{$gid} = $AA;
	my ($gname,$strand,$cdsA) = @{$GeneDat{$chr}{$gid}};
	my ($met,$len,$tmp) = (0,0);

	my $REFseq = $cDNA{$gid};

	for (@{$cdsA}) {
		last if $met;
		my ($s,$e) = @$_;
		if ($pos >= $s and $pos <= $e) {
			$met = 1;
			#$len += $pos - $s +1;
			#$seq = $cDNA{$gid};
			if ($strand eq '+') {
				$len += $pos - $s +1 -1;
			} elsif ($strand eq '-') {
				$len += $e - $pos +1 -1;
			} else {die;}
			substr($seq,$len,1) = $gt;
			#$seq = substr($cDNA{$gid},0,$len) . $gt . substr($cDNA{$gid},$len+1)
			my $test = substr($REFseq,$len,1);
			$test =~ tr/ATCG/TAGC/ if $strand eq '-';
			warn "$gid($gname),$strand: $test ne $ref" if $test ne $ref;
		} else {
			$len += $e -$s + 1;
		}
	}
	return $seq;
}
sub cmpstr {
	my ($a, $b) = @_;
	my $c = $a ^ $b;
	my @ret;
	while ($c =~ /[^\0]/g) {
		my $p = pos($c);
		push @ret,[$p,substr($a,$p-1,1),substr($b,$p-1,1)];
	}
	@ret;
}
