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
my %invert; @invert{ qw[ A C G T ] } = qw[ T G C A ];

sub getTE($) {
	my $in = $_[0];
	my %ret;
	#print $in,"\n";
	for my $p1 ( 1 .. length( $in ) -2 ) {
		next unless substr( $in, $p1, 1 ) eq $invert{ substr $in, $p1+1, 1 };
		my $pals = 0;
		for my $p2 ( 1 .. $p1 -1 ) {
			last unless substr( $in, $p1-$p2, 1 ) eq $invert{ substr $in, $p1+$p2+1, 1 };
			++$pals;
		}
		#if( $pals ) {
		if( $pals > 1 ) {
			#printf "%s at %d\n",substr( $in, $p1-$pals, ($pals+1)*2 ), $p1-$pals;
			my $palindromicDNA = substr( $in, $p1-$pals, ($pals+1)*2 );
			my $leftPad = $p1-$pals;
			$leftPad = 0 if $leftPad > 272;
			printf "%s%s @ %d [%d]", '.'x$leftPad,
				$palindromicDNA, $p1-$pals, ($pals+1)*2;
			push @{$ret{$palindromicDNA}},$p1-$pals;
			if (scalar @{$ret{$palindromicDNA}} > 1) {
				print BOLD,RED,' *',RESET;
			}
			print "\n";
		}
	}
	return \%ret;
}

my @CDS = (53565, 54167);
my $getIN = 5;
my $getLen = 6000;

($getIN,$getLen) = (0,2000);

my $wholeLen = $getIN + $getLen;
my ($leftTEs,$rightTEs,%TEs);
sub analyseTE($);

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
	my $left = substr $seq,($CDS[0]-1-$getLen),$wholeLen;
	my $right = substr $seq,($CDS[1]-$getIN),$wholeLen;
	#local $Term::ANSIColor::AUTORESET = 1;
	print substr($left,0,$getLen),BOLD,GREEN,substr($left,$getLen),RESET,"\n";
	$leftTEs = getTE($left);
	print BOLD,GREEN,substr($right,0,$getIN),RESET,substr($right,$getIN,$getLen),"\n";
	$rightTEs = getTE($right);
	analyseTE($seq);
}
close I;

sub analyseTE($) {
	my $seq = $_[0];
	for my $k (keys %{$leftTEs}) {
		$TEs{$k} = [ ${$leftTEs}{$k} ];
	}
	for my $k (keys %{$rightTEs}) {
		push @{$TEs{$k}},${$rightTEs}{$k};
	}
	for my $k (keys %TEs) {
		next unless scalar @{$TEs{$k}} == 2;
		print '>',$k,"\t";
		ddx $TEs{$k};
	}

	my @Patterns = sort { length($a)<=>length($b) || $a cmp $b } keys %TEs;
	for my $k (@Patterns) {
		my @itsPoses;
		my $itsLen = length($k);
		print "\n$k\[$itsLen]: ";
		while ($seq =~ /(?=$k)/g) {	# http://www.perlmonks.org/?node_id=1090633
			my $p = pos $seq;
			print "$p\t";
		}
		print "\n";
	}
}
__END__

#http://repo.hackerzvoice.net/depot_madchat/esprit/texture/hallucinati/finding%20DNA%20palindroms.htm

#!/usr/bin/perl
$filename = "668plasmids.fa";
open (TEXT, $filename)||die"Cannot";
$line = " ";
$count = 0;
for $n (5..20)
   {
   $re = qr /[CAGT]{$n}/;
   $regexes[$n-5] = $re;
   }
NEXTLINE: while ($count < 1000)
   {
   $line = <TEXT> ;
   $count++;
   foreach my $value (@regexes)
      {
      $start = 0;
      while ($line =~ /$value/g)
         {
         $endline = $';
         $match = $&;
         $revmatch = reverse($match);
         $revmatch =~ tr/CAGT/GTCA/;
         if ($endline =~ /^([CAGT]{0,15})($revmatch)/)
            {
            $start = 1;
            $palindrome = $match . "*" . $1 . "*" . $2;
            $palhash{$palindrome}++;
            }
         }
      if ($start == 0)
         {
         goto NEXTLINE;
         }
      }
   }
close TEXT;
while(($key, $value) = each (%palhash))
   {
   print "$key => $value\n";
   }
exit;
