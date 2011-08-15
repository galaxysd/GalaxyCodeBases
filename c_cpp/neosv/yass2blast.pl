#!/usr/bin/perl -w
#
# a quick&dirty script to convert "yass -d 1 " output files into
# - axt format
# - fasta format
# - blast format
# usage: yass2blast.pl [-blast|-fasta|-axt] {yassoutputfiles}\n";

use strict;

#
# Print an alignment according to $format = {"fasta",blast","blast header","axt"}
#

sub PrintAlign() {
    # 1) read parameters
    my ($nbAlignments,$nbSSequences,$pos_q_str,$pos_q_end,$pos_s_str,$pos_s_end,$name_q,$name_s,$size_q,$size_s,$reverse,$evalue,$score,$bitscore,$bias,$ts,$tv,$proba,$entropy,$fasta1,$middle,$fasta2,$format,$newchunk_q,$newchunk_s) = @_ ;

    # 2) reverse and complement sequences to keep the "query order"
    if ($reverse eq "-") {
	# change first file positions (increasing order)
	my $tmp = $pos_q_str;
	$pos_q_str = $pos_q_end;
	$pos_q_end = $tmp;
	# sequence complement
	$fasta1 =~ tr /ATUGCRYKMSWVBDHatugcrykmswvbdh/TAACGYRMKWSBVHDtaacgyrmkwsbvhd/;
	$fasta1 = reverse($fasta1);
	$fasta2 =~ tr /ATUGCRYKMSWVBDHatugcrykmswvbdh/TAACGYRMKWSBVHDtaacgyrmkwsbvhd/;
	$fasta2 = reverse($fasta2);
	$middle = reverse($middle);
    }
    


    # 3) format selection
    if ($format eq "fasta") { 
        # 3.1) fasta but with "-" to keep alignments
	print ">".$nbAlignments."a||".$name_q.": [".$pos_q_str."-".$pos_q_end ."]".$reverse."\n";
	print $fasta1."\n";
	print ">".$nbAlignments."b||".$name_s.": [".$pos_s_str."-".$pos_s_end ."]".$reverse."\n";
	print $fasta2."\n";

    } elsif ($format eq "axt") { 
        # 3.2) axt (blastz) format
	print $nbAlignments."\t".$name_q."\t".$pos_q_str."\t".$pos_q_end ."\t".$name_s."\t".$pos_s_str."\t".$pos_s_end ."\t".$reverse."\t".$score."\n";
	print $fasta1."\n";
	print $fasta2."\n";
	
    } elsif ($format eq "blast") { 
        # 3.3) blast like format
	$middle =~ tr /:\./  /;
	my $nbMatches = ($middle =~ tr/\|//);
	# a) print Sequence Header
	if ($newchunk_s){
	    print ">".$nbSSequences."_0 ".$name_s."\n";
	    print "          Length = ".$size_s."\n";
	    print "\n";
	}
	# b) print Alignment Header
	print " Score = ".$bitscore." bits (".$score."), Expect = ".$evalue."\n";
	print " Identities = ".$nbMatches."/".(length($middle))." (".(int(100*$nbMatches/length($middle)))."%)\n";
	print " Strand = Plus / ";
	if ($reverse eq "-"){
	    print "Minus"."\n";
	}else{
	    print "Plus"."\n";
	}
	print "\n\n";
	

	# c) print Alignment
	my $c = 0;
	my $i_fasta1 = $pos_q_str;
	my $i_fasta2 = 0;
	if ($reverse eq "-") {
	    $i_fasta2 = $pos_s_end;
	} else {
	    $i_fasta2 = $pos_s_str;
	}
	while ($c < length($middle)){
	    my $subfasta1 = substr($fasta1,$c,60);
	    my $submiddle = substr($middle,$c,60);
	    my $subfasta2 = substr($fasta2,$c,60);
	    my $letters_subfasta1 = $subfasta1; $letters_subfasta1 =~ s/-//g;
	    my $letters_subfasta2 = $subfasta2; $letters_subfasta2 =~ s/-//g;
	    my $nbletters_subfasta1 = length($letters_subfasta1);
	    my $nbletters_subfasta2 = length($letters_subfasta2);
	    printf ("Query: %-9d ",$i_fasta1); $i_fasta1 += $nbletters_subfasta1-1;
	    print($subfasta1);
	    printf (" %d\n",       $i_fasta1); $i_fasta1 += 1;

	    printf ("                 ");  
	    print($submiddle);
	    printf ("\n");

	    if ($reverse eq "-") {
		printf ("Sbjct: %-9d ",$i_fasta2); $i_fasta2 -= $nbletters_subfasta2-1;
		print($subfasta2);
		printf (" %d\n",       $i_fasta2); $i_fasta2 -= 1;
	    } else {
		printf ("Sbjct: %-9d ",$i_fasta2); $i_fasta2 += $nbletters_subfasta2-1;
		print($subfasta2); 
		printf (" %d\n",       $i_fasta2); $i_fasta2 += 1;
	    }
	    $c += 60;	  
	    print "\n\n";
	}
    } elsif ($format eq "blastheader") {
	# 3.4) just the header of blast
	if ($newchunk_s){
	    if (length ($name_s) > 60) {
		print "".$nbSSequences."_0 ".(substr($name_s,0,60))."...   ".$bitscore."   ".$evalue."\n";
	    } else {
		print "".$nbSSequences."_0 ".$name_s;
		for(my $e = length ($name_s); $e <= 63 ;$e ++) {
		    print " ";
		}
		print "   ".$bitscore."   ".$evalue."\n";
	    }
	}#1_0 embl|AF417609|AF417609 Colpidium campylum telomerase RNA gen...    32   0.063
    }
}


sub ScanFile($$) {

    my ($yassOutputFile,$format) = @_ ;

    open(FIC,$yassOutputFile) or die print "cant open ".$yassOutputFile."\n";

    my $nbAlignments = 0;
    my $nbSSequences = 0;
    my $selector     = 0;
    my $previous_name_s = "";
    my $previous_name_q = "";

    #A) alignments
    my ($pos_q_str,
	$pos_q_end,
	$pos_s_str,
	$pos_s_end,
	$name_q,
	$name_s,
	$size_q,
	$size_s,
	$reverse,
	$evalue,
	$score,
	$bitscore,
	$bias,
	$ts,
	$tv,
	$proba,
	$entropy,
	$fasta1,
	$middle,
	$fasta2,)
	=
       (
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"");
    
    while (my $line = <FIC>){

	# potential header (for alignemnts)
	if ($format eq "blastheader" && !($name_q eq $previous_name_q)) {
	    print "Query= ".$name_q."\n";
	    print "     (".$size_q." letters)\n\n\n";
	    print "                                                                 Score    E\n";
	    print "Sequences producing significant alignments:                      (bits) Value\n\n";
	    $previous_name_q = $name_q;
	}
	
        # lines begining with "*"
	if($line  =~ /^\*/ ){
	    
	    # lines begining with "*("
	    if($line  =~ /^\*\(/  ) {
	
		# print previous alignments
		if (!($middle eq "")){
		    &PrintAlign($nbAlignments,$nbSSequences,
				$pos_q_str,$pos_q_end,$pos_s_str,$pos_s_end,
				$name_q,$name_s,$size_q,$size_s,
				$reverse,$evalue,$score,$bitscore,
				$bias,$ts,$tv,$proba,$entropy,
				$fasta1,$middle,$fasta2,$format,
				($previous_name_q eq "") || !($previous_name_q eq $name_q),
				($previous_name_s eq "") || !($previous_name_s eq $name_s)
			);
		    $fasta1 = "";
		    $middle = "";
		    $fasta2 = "";
		}
		
		
		# count number of new alignments (and new sequences)
		$nbAlignments++;
		if (($previous_name_s eq "") || 
		    !($previous_name_s eq $name_s)){
		    $nbSSequences++;
		}


		#first line :
                #*(546486-566813)(515659-536310) Ev: 0 s: 20328/20652 f
		$line  =~ m/\(([0-9]+)-([0-9]+)\)\(([0-9]+)-([0-9]+)\) Ev: (.+) s: ([0-9]+)\/([0-9]+) ([rf])$/;
                # get positions
		$pos_q_str  = $1;
		$pos_q_end  = $2;
		$pos_s_str  = $3;
		$pos_s_end  = $4;
		$evalue     = $5;
		$size_q     = $6;
		$size_s     = $7;
		$reverse    = $8;
		if ($reverse =~ /^f/) {
		    $reverse = "+";
		} else {
		    $reverse = "-";
		}
	   } elsif ($line  =~ /^\* score/ ){
	       # third line :
	       $line  =~  m/([0-9]+).*bitscore = (.*)/;
	       $score    = $1 + 0;
	       $bitscore = $2 + 0;
	   } elsif ($line =~ /^\* mutations/){
                # fourth line :
	        #* mutations per triplet 129, 103, 124 (4.54e-04) | ts : 190 tv : 166 | entropy : 5.88743	       
 		$line  =~ m/^\* mutations per triplet ([0-9]+, [0-9]+, [0-9]+) \((.*)\) \| ts : ([0-9]+) tv : ([0-9]+) \| entropy : (.*)/;
 		$bias    = $1;
		$proba   = $2;
		$ts      = $3;
		$tv      = $4;
		$entropy = $5;
	    }else{
		# second line :
		$line  =~ m/\*[ ]*"(.+)" \(([0-9]+) bp\)[ ]*\/[ ]*"(.+)" \(([0-9]+) bp\)[ ]*$/;
		$previous_name_q = $name_q;
		$name_q = $1;
		$size_q = $2;
		$previous_name_s = $name_s;
		$name_s = $3;
		$size_s = $4;
	    }
	} else {
	    # alignment :
	    if($line =~ /^[A-Za-z-]+$/){
		$line =~ s/\n//;
                $selector++;
                if ($selector % 2){
                    $fasta1 .= "".$line;
                } else {
                    $fasta2 .= "".$line;
                }
            } elsif ($line =~ /^[| :.]+$/ && (length($fasta1) > length($fasta2)) ){  
		$line =~ s/\n//;
		$middle .= $line;
	    }
	}
    } # while <FIC>

    if (!($middle eq "")){
	&PrintAlign($nbAlignments,$nbSSequences,
		    $pos_q_str,$pos_q_end,$pos_s_str,$pos_s_end,
		    $name_q,$name_s,$size_q,$size_s,
		    $reverse,$evalue,$score,$bitscore,
		    $bias,$ts,$tv,$proba,$entropy,
		    $fasta1,$middle,$fasta2,$format,
		    ($previous_name_q eq "") &&!($previous_name_q eq $name_q),
		    ($previous_name_s eq "") &&!($previous_name_s eq $name_s)
	    );
	$fasta1 = "";
	$middle = "";
	$fasta2 = "";
    }
    close FIC;
    return $nbAlignments;
}





#------#
# Main #
#------#

my $outputformat = "blast";

($#ARGV >= 0) or die "yass2blast:\n".
                     "  a quick&dirty script to convert \"yass -d 1 \" output files into\n".
                     "  - axt format\n".
                     "  - fasta format\n".
                     "  - blast format\n".
                     "\n".
                     "  usage: yass2blast.pl [-blast|-fasta|-axt] {yassoutputfiles}\n".
                     "\n";


for (my $i = 0 ; $i <= $#ARGV ; $i++) {
    #a) select on the fly parameters ...
    if (($ARGV[$i]) eq "-blast") {
	$outputformat = "blast";
    } elsif (($ARGV[$i]) eq "-fasta") {
	$outputformat = "fasta";
    } elsif (($ARGV[$i]) eq "-axt") {
	$outputformat = "axt";
    } else {
	#b) or scan files
	if ($outputformat eq "blast") {
	    &ScanFile($ARGV[$i],"blastheader");
	    print "\n";
	}
	&ScanFile($ARGV[$i],$outputformat);
    }
}
