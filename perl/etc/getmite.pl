#!/usr/bin/env perl
use strict;
use warnings;
use Term::ANSIColor qw(:constants);
use Data::Dump qw(ddx);

=pod
>theseq
ATTACATGGACTTCTCTGGCTACCAAT

blastn to bacteria (taxid:2)

http://www.ncbi.nlm.nih.gov/nucleotide/373881580?report=genbank&log$=nuclalign&blast_rank=1&RID=927901S401N
LOCUS	  JQ031552			76221 bp    DNA	circular BCT 29-JAN-2012
DEFINITION  Aliivibrio fischeri strain KB1A-97 plasmid pKB1A97-67, complete
		  sequence.
ACCESSION   JQ031552
VERSION	JQ031552.1  GI:373881580
KEYWORDS    .
SOURCE	 Aliivibrio fischeri (Vibrio fischeri)
  ORGANISM  Aliivibrio fischeri
		  Bacteria; Proteobacteria; Gammaproteobacteria; Vibrionales;
		  Vibrionaceae; Aliivibrio.
REFERENCE   1  (bases 1 to 76221)
  AUTHORS   Summers,A.O., Wireman,J. and Williams,L.E.
  TITLE	Direct Submission
  JOURNAL   Submitted (14-NOV-2011) Department of Microbiology, University of
		  Georgia, 527 Biological Sciences, Athens, GA 30602, USA
FEATURES		   Location/Qualifiers
	source		1..76221
				 /organism="Aliivibrio fischeri"
				 /mol_type="genomic DNA"
				 /strain="KB1A-97"
				 /host="Euprymna scolopes"
				 /db_xref="taxon:668"
				 /plasmid="pKB1A97-67"
				 /country="USA: Hawaii"
				 /collection_date="2005"
				 /note="strain provided by Eric Stabb, University of
				 Georgia, Athens, GA, USA"
	CDS		   complement(53120..53287)
				 /codon_start=1
				 /transl_table=11
				 /product="hypothetical protein"
				 /protein_id="AEY78253.1"
				 /db_xref="GI:373881650"
				 /translation="MIFVCNGLQFYCLVNRHESVINFQACLWVILVLFIYEKWHYVDLRRFNGCFEWLS"
>>>  CDS		   53565..54167
				 /codon_start=1
				 /transl_table=11
				 /product="Transposase-like protein"
				 /protein_id="AEY78254.1"
				 /db_xref="GI:373881651"
				 /translation="MDFSGYQYPSDIILQAVRYYVSYKLSTRDIEEIFTERGSAIDHS
				 TINRWVITFAPMLEQNARQLKRKVSSSWRMDETYIKIKGEWWYYYRAVDKYGDIVDFY
				 LSKERDEKAAKAFLRKAIHTNGLPDKVVIDKSGANALALHNLNVKLWLSVVFMLNLIE
				 IVDVKYLNNIVEQSYRPIKQKMVQALGWKSVEGATVTMSG"
	CDS		   54975..55841
				 /codon_start=1
				 /transl_table=11
				 /product="hypothetical protein"
				 /protein_id="AEY78255.1"
				 /db_xref="GI:373881652"
				 /translation="MQHDVGSNSDVYISCYSLNIKESGLSGAHYTISDIRKSIETVKVTSSYRHLHIEMKNSLCVCFSTLDITWVNGLEHGELLAGQFTVFDSECPISYKVTRVGSLCFVFIPKYFYEGVLQKQMVRCGMFEFVYVDAIRFILTRVNSKEDGEQLLISELLALGYLLSVLERKEEAVGKKVAFEDKVHEVIKDNMLNPSLYLDDIALILGCSKRKIQHCLSLQGVSFTKLVTKYRIEYLAEQLIRKKHSRIDVLCYESGFNSPGYASNSFKVIMGMSPKEYRCRYLAKSSVF"
=cut

# http://www.perlmonks.org/?node_id=1085446
#my %invert; @invert{ qw[ A C G T ] } = qw[ T G C A ];

my @CDS = (53565, 54167);
my $irMinLen = 3;
my $getLen = 40099;

$getLen = 150;
$irMinLen = 5;

# http://www.perlmonks.org/?node_id=197793
sub revdnacomp {
  my $dna = shift; # or   my $dna = shift @_;
  my $revcomp = reverse($dna);
  $revcomp =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
  return $revcomp;
}

sub analyseIR($$) {
	my ($seq1,$seq2) = @_;
	my $seqlen1 = length $seq1;
	my $seqlen2 = length $seq2;
	my %Dat;
	for my $irLeftPos2 (0 .. ($seqlen2-1)) {
		#for my $irLen ($irMinLen .. int(($seqlen-$irLeftPos)/2)) {
		for my $irLen ($irMinLen .. ($seqlen2-1-$irLeftPos2)) {
			my $irSeq = substr $seq2,-$irLen-$irLeftPos2,$irLen;
			my $revcmp = revdnacomp($irSeq);
			my @thePoses;
			my $p = 0;
			while ($p >= 0) {
				$p = index $seq1,$revcmp,$p+1;
				push @thePoses,$p;
			}
			last if $thePoses[0] == -1;
			pop @thePoses;
			$Dat{"$_|$irLeftPos2"} = $irLen for @thePoses;
=pod
			print "$irLen\[$irSeq]\t";
			for (@thePoses) {
				my $found = substr $seq1,$_,$irLen;
				my $extra = substr $seq1,$_+$irLen,6;
				print ".$found",GREEN,BOLD,"$extra $_",RESET;
			}
			print "\t$irLeftPos2\n";
=cut
		}
	}
	#print "1:$seq1,$seq2\n";
	return \%Dat;
}
sub printRes($$$$) {
	my ($dat,$getLen,$seq1,$seq2) = @_;
	my %res;
	for my $poses (keys %{$dat}) {
		my $irlen = $dat->{$poses};
		++$res{$irlen}{$poses};
	}
	for my $irlen (sort {$a<=>$b} keys %res) {
		my %tmp;
		for my $poses (keys %{$res{$irlen}}) {
			my @Poses = split /\|/,$poses;
			my $dis = $Poses[1] + $Poses[0];
			my $upperPoses = join('|',$Poses[0]-1,$Poses[1]-1);
			if (exists $res{$irlen+1}) {
				delete $res{$irlen}{$poses} if exists $res{$irlen+1}{$upperPoses};
			}
			$tmp{$poses} = $dis;
		}
		my @Poses = keys %{$res{$irlen}};
		if (@Poses > 0) {
			my @tmp = sort { $tmp{$a} <=> $tmp{$b} } @Poses;
			$res{$irlen} = \@tmp;
		} else {
			delete $res{$irlen};
		}
	}
	for my $irlen (sort {$a<=>$b} keys %res) {
		print "--- $irlen ---\n";
		#ddx $res{$irlen};
		for my $poses ( @{$res{$irlen}} ) {
			my @Poses = split /\|/,$poses;
			my $irL = substr $seq1,$Poses[0],$irlen;
			my $irR = substr $seq2,-$irlen-$Poses[1],$irlen;
			#my $revcmp = revdnacomp($irR);
			print "$poses $irL $irR\n";
		}
	}
	#print "2:$seq1,$seq2\n";
}

open I,'<','pKB1A97-67.fa' or die $?;
while (<I>) {
	chomp;
	my $Head;
	my ($id,$desc)=split / /,$_,2;
	if ($desc && $desc !~ /^\s*$/) {
		$desc=~s/\t/_/g;
		$Head="$id $desc";
	} else { $desc='.';$Head=$id; }
	$/=">";
	my $seq=<I>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/="\n";
	my $len=length($seq);
	print STDERR "$id:$len,[$desc]\n";
	$seq = uc $seq;
	my $left = substr $seq,0,($CDS[0]-1);
	my $right = substr $seq,($CDS[1]);
	print "0 Len:(",length $left,',',length $right,")\n";
	#$left = substr $left,-$getLen;
	#$right = substr $right,0,$getLen;
	print "1 Len:(",length $left,',',length $right,")\n";
	my $res = analyseIR($left,$right);
	#ddx $res;
	printRes($res,$getLen,$left,$right);
=pod
	#print "[$left]\n\n[$right]\n";
	[TTTATTGACCGGCAGCGGTATGCGAGTGCGAACACCAACTTCTTATGGTGAATCACAAATGTACGACCAGCCTTTGCCTGTCACGCAGATGGACATGGGGCAAATTCCT
...	ATACACTTACCAACGCATCTTAATAGGGGTTCTGTCGCATGACAGGGAACTAAATCATCCATGAGATACTCAAAGTTCTGTTTTTTACATTAACTCGAGTATTAC]

	[GAGGTCTGGACCCAAATTAAACGAAGGCAAGTTGGGGACGTAAATTTACCCGTATGGGAGCGCTTTCACGCGCTCGCTGCATAATTGTGTCTAGAATTCGGCACGTCTA
...	TGGATGCACGTATGCCTGCAATGGCAAATCTGACCGTAAAAGATGGGTTGGATTT]
=cut
}
close I;

__END__
