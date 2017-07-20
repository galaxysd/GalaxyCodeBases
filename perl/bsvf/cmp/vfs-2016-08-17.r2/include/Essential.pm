###########################################################################
##  ViralFusionSeq
##  Software for discovering and annotating viral integration event and
##    fusion transcript
##
##  Version 2.0
##
##  Copyright (C) 2012 by Jing-Woei Li & Raymond Wan, All rights reserved.
##  Contact:  marcoli@cuhk.edu.hk, rwan@cuhk.edu.hk
##  Organization:  Hong Kong Bioinformatics Centre, School of Life Sciences, The
##                 Chinese University of Hong Kong, Shatin, NT,
##                 Hong Kong SAR
##
##  This file is part of ViralFusionSeq.
##
##  ViralFusionSeq is free software; you can redistribute it and/or
##  modify it under the terms of the GNU General Public License
##  as published by the Free Software Foundation; either version
##  3 of the License, or (at your option) any later version.
##
##  ViralFusionSeq is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public
##  License along with ViralFusionSeq; if not, see
##  <http://www.gnu.org/licenses/>.
###########################################################################


##  $LastChangedDate: 2013-12-13 17:45:55 +0800 (Fri, 13 Dec 2013) $
##  $LastChangedRevision: 1282 $
package Essential;

use diagnostics;
use strict;
use warnings;

$| = 1; #disable buffer

use Cwd;
use Exporter;
use Bio::DB::Sam;
use Bio::SeqIO;
use Bio::SearchIO;
use File::Copy;

use FindBin;
use lib "$FindBin::Bin/include";
use General; # Module for general subroutines
use BWAauto; #BWA sampe/samse pipeline
use PerlIO::gzip;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = qw(bwa_index bwa_sw rtn_CS extract_fa_by_len cap3_run blast_run blast_formatdb ssake_reads_prep ssake_run map_2_genome create_simplified_readsID_fa rtn_reads_cap3_contigs file_2_array file_2_array_skip_header blast_parse CS_read_level_analysis tx_targeted_assembly Assemble_CS clean_tempfiles rtn_fastq_phred);
our @EXPORT_OK   = qw(bwa_index bwa_sw rtn_CS extract_fa_by_len cap3_run blast_run blast_formatdb ssake_reads_prep ssake_run map_2_genome create_simplified_readsID_fa rtn_reads_cap3_contigs file_2_array file_2_array_skip_header blast_parse CS_read_level_analysis  tx_targeted_assembly Assemble_CS clean_tempfiles rtn_fastq_phred);
our %EXPORT_TAGS = ( DEFAULT => [qw(&bwa_index &bwa_sw &split_n_take_first &rtn_CS &extract_fa_by_len &cap3_run &blast_run &blast_formatdb &ssake_reads_prep &ssake_run &map_2_genome &create_simplified_readsID_fa &rtn_reads_cap3_contigs &file_2_array &file_2_array_skip_header &blast_parse &CS_read_level_analysis &tx_targeted_assembly &Assemble_CS &clean_tempfiles &rtn_fastq_phred)]);

##  stderr of the indexing phase is sent to /dev/null
sub bwa_index {
	my $curr_process =(caller(0))[3];
	print STDERR "Entering $curr_process.\n";
	my $mapper = $_[0];
	my $genome = $_[1];
	my $small_genome = $_[2];
	my $logdir = $_[3];

	my @bwa_index = split (/ /, "$mapper index -a bwtsw $genome 2>/dev/null");
	if ($small_genome) {
		@bwa_index = split (/ /, "$mapper index -a is $genome 2>/dev/null");
	}
	ExecCmd (\@bwa_index, $DEBUG_MODE);

	return;
}

sub bwa_sw {
	my $curr_process =(caller(0))[3];
	print STDERR "Entering $curr_process.\n";
	#my $dir = getcwd; #main dir
	my $runID = $_[0];
	my $mapper = $_[1];
	my $F_readfile_f = $_[2];
	my $R_readfile_f = $_[3];
	my $bwa_index = $_[4];
	my $thread = $_[5];
	my $logdir = $_[6];
	#my $mapper = "$dir/thirdparty/bwa-0.6.1/bwa bwasw";
	#my $bwasw = "$mapper -t $thread $bwa_index $F_readfile_f $R_readfile_f -f $runID.sam 2>$runID.bwasw.log";
	my @bwasw = ();
	if ($R_readfile_f eq ""){
		@bwasw = split (/ /, "$mapper bwasw -t $thread $bwa_index $F_readfile_f 2>>$logdir/$runID.bwasw.log | samtools view -bSu -T $bwa_index - | samtools sort -m 2000000000 - $runID");
	} else {
		@bwasw = split (/ /, "$mapper bwasw -t $thread $bwa_index $F_readfile_f $R_readfile_f 2>>$logdir/$runID.bwasw.log | samtools view -bSu -T $bwa_index - | samtools sort -m 2000000000 - $runID");
	}
	ExecCmd (\@bwasw, $DEBUG_MODE);
	#my $cmd = "samtools view -bSu $runID.sam | samtools sort -m 2000000000 - $runID"; #$runID.bam
	#exec_cmd($cmd, 1);
	#my $cmd2 = "rm $runID.sam";
	#exec_cmd($cmd2, 1);
	#exec_cmd("rm -rf $runID.bwasw.log", 1);
	my $fileout = "$runID.bam";
	return $fileout;
}

sub rtn_CS {
	#May be good to check another end of clipped read can be mapped to another genome
	my $curr_process =(caller(0))[3];
	print STDERR "Entering $curr_process.\n";
	my $bam_file = $_[0];
	my $fa_file = $_[1];
	my $CS_file = $_[2];
	my $MS_file = $_[3];
	my $CS_MS_read_file = $_[4];
	my $F_read_file = $_[5];
	my $R_read_file = $_[6];

	if (-z $bam_file){
		print STDERR "BAM file: $bam_file is 0 byte. VFS will terminate\n";
		exit(1);
	}
	print STDERR "Processing alignment file: $bam_file\n";
	open my $CS_MS_fullread_relation_wFH, ">", $CS_MS_read_file or die "Can't open $CS_MS_read_file for writing : $!";
	my $CS_MS_mate_fullread_relation_wFH;
	if (defined $R_read_file){
		my $file_out = "$CS_MS_read_file.mates.fa";
		open $CS_MS_mate_fullread_relation_wFH, ">", $file_out or die "Can't open $CS_MS_read_file for writing : $!"; #Mates of CS-read
	}
	open my $CS_wFH, ">", $CS_file or die "Can't open $CS_file for writing : $!";
	open my $MS_wFH, ">", $MS_file or die "Can't open $MS_file for writing : $!";

	#1st pass to record the readID need to be hashed
	print STDERR "First pass on $bam_file\n";
	my %readID_needed = ();
	my $line_count = 0;
	my $dot_count = 0;
	my $remainder = 5000000;
	my $line_dot_count = 60;
	my $bam = Bio::DB::Sam->new(-bam  =>"$bam_file",
								-fasta=>"$fa_file",
								);
	my $aln  = $bam->features(-iterator=>1);
	while (my $aln_read = $aln->next_seq) {
		my $ifSC = 0;
		my @AUX_tags = $aln_read->get_all_tags;
		my @CIGAR = @{$aln_read->cigar_array};
		#test if soft clipped
		for my $ref (@CIGAR){
			my @cigarA = @{$ref};
			my $SCq = "S";
			if (grep {$_ eq $SCq} @cigarA) { #assertain if a read is soft-clipped
				#print "Clipped length: $cigarA[1]\n";
				$ifSC = 1;
			}
		}
		if ($ifSC == 1){ #extract Clipped sequence and its associated mapped sequence
			my $read_id = $aln_read->query->name;
			#print "need $read_id\n";
			$readID_needed{$read_id} = "";
		}
		$line_count++;
		if (($line_count % $remainder) == 0){
			print ".";
			$dot_count++;
			if (($dot_count % $line_dot_count) == 0){
				print "\n";
				$line_count = 0;
				$dot_count = 0;
			}
		}
	}
	
	#reset count for Second pass
	$line_count = 0;
	$dot_count = 0;

	my %F_readID_to_seq_h = %{hash_fastq_readID_to_sequence_only_needed($F_read_file, \%readID_needed)};
	my %R_readID_to_seq_h = ();
	if (defined $R_read_file){
		%R_readID_to_seq_h = %{hash_fastq_readID_to_sequence_only_needed($R_read_file, \%readID_needed)};
	}

	print STDERR "Extracting Clipped sequences\n";
	my $bam2 = Bio::DB::Sam->new(-bam  =>"$bam_file",
								-fasta=>"$fa_file",
								);
	my $aln2  = $bam2->features(-iterator=>1);
	while (my $aln_read = $aln2->next_seq) {
		my $ifSC = 0;
		my @AUX_tags = $aln_read->get_all_tags;
		my @CIGAR = @{$aln_read->cigar_array};
		#test if soft clipped
		for my $ref (@CIGAR){
			my @cigarA = @{$ref};
			my $SCq = "S";
			if (grep {$_ eq $SCq} @cigarA) { #assertain if a read is soft-clipped
				#print "Clipped length: $cigarA[1]\n";
				$ifSC = 1;
			}
		}
		if ($ifSC == 1){ #extract Clipped sequence and its associated mapped sequence
			my $read_id = $aln_read->query->name;
			my $ref_id = $aln_read->seq_id;
			my $ref_dna = $aln_read->dna;
			my $q_len = $aln_read->query->length; #this is alignment length (i.e. how many nt aligned from the read)
			my $q_st = $aln_read->query->start; $q_st--; #1-based coordinate
			my $q_ed = $aln_read->query->end; #1-based coordinate
			my $read_nt = $aln_read->query->dna;
			my $bitwise_flag = $aln_read->flag;
			#print "$read_id\t$bitwise_flag\n";
			#print "q start: $q_st\tq end: $q_ed\tq_len: $q_len\n";
			#my $F_or_R_read;
			if (defined $R_read_file){ #Paired-end data
				if (($bitwise_flag & 0x40) == 0x40){ #Forward read
					#$F_or_R_read = "F";
					#print "$read_id\t-$F_readID_to_seq_h{$read_id}-\n";
					print $CS_MS_fullread_relation_wFH '>'.$read_id."\n".$F_readID_to_seq_h{$read_id}."\n";
					print $CS_MS_mate_fullread_relation_wFH '>'.$read_id."\n".$R_readID_to_seq_h{$read_id}."\n";
				}
				else { #Reverse read
					#$F_or_R_read = "R"
					#print "$read_id\t-$R_readID_to_seq_h{$read_id}-\n";
					print $CS_MS_fullread_relation_wFH '>'.$read_id."\n".$R_readID_to_seq_h{$read_id}."\n";
					print $CS_MS_mate_fullread_relation_wFH '>'.$read_id."\n".$F_readID_to_seq_h{$read_id}."\n";
				}
			} else { #Single-end data
				unless (($bitwise_flag & 0x4) == 0x4){ # 0x4 is unmapped. For SE data, if it's mapped, extract them
					print $CS_MS_fullread_relation_wFH '>'.$read_id."\n".$F_readID_to_seq_h{$read_id}."\n";
				}
			}

			my $aligned_nt = substr $read_nt, $q_st, $q_len;
			#get clipped seq
			my $c_seq = "";
			my $m_seq = "";
			#S(\d+)M(\d+)
			if ($q_st > 0){
				$c_seq = substr $read_nt, 0, $q_st;
			}
			else { #M(\d+)S(\d+)
				#get number of matches in front
				my @cigar_array = @{$aln_read->cigar_array};
				my $ext = 0;
				for (@cigar_array){
					my @a = @{$_};
					if ($a[0] eq "M"){
						$ext = $a[1];
					}
				}
				$c_seq = substr $read_nt, $q_ed, $ext;
			}
			#get mapped sequence
			$m_seq = substr $read_nt, $q_st, $q_len;

			#print "Read ID: $read_id\n";
			#print "Read nt: $read_nt\n";
			#print "Query length: $q_len\n";
			#print "Query start (0-based): $q_st\n";
			#print "Query end (0-based): $q_ed\n";
			#print "Aligned nt: $aligned_nt\n";
			#print "Clipped nt: $c_seq\n";
			my $CS_data_out = ">$read_id\_$ref_id\n$c_seq"; #readID + underscore + ref ID + sequence
			#print "$data_in\n";
			print $CS_wFH "$CS_data_out\n";
			#output mapped seq file
			my $MS_data_out = ">$read_id\_$ref_id\n$m_seq"; #readID + underscore + ref ID + sequence
			print $MS_wFH "$MS_data_out\n";
		}
		$line_count++;
		if (($line_count % $remainder) == 0){
			print ".";
			$dot_count++;
			if (($dot_count % $line_dot_count) == 0){
				print "\n";
				$line_count = 0;
				$dot_count = 0;
			}
		}
	}
	close ($CS_wFH);
	close ($MS_wFH);
}

##  Input file can be either gzipped or not
sub hash_fastq_readID_to_sequence_only_needed {
	my $file = $_[0];
	my %readID_needed = %{$_[1]};
	print STDERR "Hashing $file\n";
	my %h = ();
	my $is_gzip_input = DetermineFileType ($file);

	my $rFQ_seq_FH;
	if ($is_gzip_input) {
		open ($rFQ_seq_FH, "<:gzip", $file) or die "Can't open $file (gzipped) for reading : $!";
	}
	else {
		open $rFQ_seq_FH, '<', $file or die "Can't open $file for reading : $!";
	}

	while (1){
		my $id = "";
		my $seq = "";
		my $aux = "";
		my $qual = "";

		$id = <$rFQ_seq_FH>;

		if (!defined $id) {
			last;  ##  End of file reached
		}

		chomp ($id);
		#my @s = split (/ /, $_);
		my @s = split (/ /, $id);
		my $full_readID = $s[0];
		my $str_readID = substr $full_readID, 1;

		$seq = <$rFQ_seq_FH>;
		$aux = <$rFQ_seq_FH>;
		$qual = <$rFQ_seq_FH>;

		chomp $seq;
		chomp $aux;
		chomp $qual;

		if (defined $readID_needed{$str_readID}){
			$h{$str_readID}=$seq;
		}
	}

	close ($rFQ_seq_FH);

	return \%h;
}

sub hash_fastq_readID_to_sequence {
	my $curr_process =(caller(0))[3];
	print STDERR "Entering $curr_process.\n";
	my $file = $_[0];
	print STDERR "hashing $file\n";
	my %h = ();
	open my $rFQ_seq_FH, "<", $file or die "Can't open $file for reading : $!";
	while (<$rFQ_seq_FH>){
		chomp;
		my @s = split (/ /, $_);
		my $full_readID = $s[0];
		my $str_readID = substr $full_readID, 1;
		my $seq = <$rFQ_seq_FH>; chomp $seq;
		my $aux = <$rFQ_seq_FH>; chomp $aux;
		my $qual = <$rFQ_seq_FH>; chomp $qual;
		$h{$str_readID}=$seq;
	}
	close ($rFQ_seq_FH);
	return \%h;
}

sub extract_fa_by_len {
	my $curr_process =(caller(0))[3];
	print STDERR "Entering $curr_process.\n";
	my $usage = qq{Usage: $0 <fasta in> <fasta.out> <min length>
	};
	die($usage) if (@_!=5);
	#print @_;
	my $file_in = $_[0];
	my $file_out = $_[1];
	my $min_len = $_[2];
	my $MS_file = $_[3];
	my $MS_file_out = $_[4];
	
	my %MS_retain_h = ();
	my $seqin  = Bio::SeqIO->new(-file => "$file_in", -format => "fasta");
	my $seqout = Bio::SeqIO->new(-file => ">$file_out", -format => "fasta");
	while(my $seq = $seqin->next_seq) {
	  if($seq->length >= $min_len) {
		$seqout->write_seq($seq);
		my $header = $seq->primary_id; #fasta header without ">"
		$MS_retain_h{$header}="";
	  }
	}
	my $seqin2  = Bio::SeqIO->new(-file => "$MS_file", -format => "fasta");
	my $seqout2 = Bio::SeqIO->new(-file => ">$MS_file_out", -format => "fasta");
	while(my $seq = $seqin2->next_seq) {
	  if(defined($MS_retain_h{$seq->primary_id})){
		$seqout2->write_seq($seq);
	  }
	}
}

sub create_simplified_readsID_fa {
	my $curr_process =(caller(0))[3];
	print STDERR "Entering $curr_process.\n";
	#hard-coded for Illumina 1.8 fastq
	my $fa_in = $_[0];
	my $fa_out = $_[1];
	my $tracking_file = $_[2];
	my %h = ();
	open my $FH, "<", $fa_in or die "Can't open $fa_in for reading : $!";
	open my $wFH, ">", $fa_out or die "Can't open $fa_out for writing : $!";
	open my $trackingFH, ">", $tracking_file or die "Can't open $tracking_file for writing : $!";
	my $c = 1;
	while (<$FH>){
		chomp;
		if (/^>/){
			#my @s = split (/:/, $_); #retain only unique info
			#my $str = join("", @s[3..($#s)]);
			#
			$_ =~ s/>//g;
			print $trackingFH $c."\t".$_."\n";
			$h{$c} = $_;
			print $wFH ">$c\n";
			$c++;
		}
		else {
			print $wFH "$_\n";
		}
	}
	close ($FH);
	close ($wFH);
	close ($trackingFH);
	return \%h;
}

sub cap3_run {
	my $curr_process =(caller(0))[3];
	print STDERR "Entering $curr_process.\n";
	#my $dir = getcwd; #main dir
	my $seq_flie_in = $_[0];
	my $cap3 = $_[1];
	my $logdir = $_[2];

	#my $cap3 = "$dir/thirdparty/CAP3/cap3";
	#my @cmd = split (/ /, "$cap3 $seq_flie_in -x cap3 -p 90 >$seq_flie_in.cap3.log 2>>$logdir/cap3.log");
	my @cmd = split (/ /, "$cap3 $seq_flie_in -x cap3 >$seq_flie_in.cap3.log 2>>$logdir/cap3.log");
	ExecCmd (\@cmd, $DEBUG_MODE);

	#parse CAP3 log file (reads' identity used in assembly. Must use simplified fa to do CAP3, otherwise the reads'ID are truncated, cannot trace back)
	my %contig_to_reads = ();
	my @a = @{file_2_array("$seq_flie_in.cap3.log")};
	my $cont_parse = 1;
	my $curr_contig;
	#my $st_new_parse = 0;
	for (@a){
		if ($cont_parse == 1){
			if (/^\*/){
				#print "$_\n";
				if (/\*+ Contig (\d*) \*+/){
					#print "Contig $1\n";
					$curr_contig = $1;
				}
			}
			if ((/^(\d*)[+|-]/)||(/(\d*)[+|-] is/)){ #fasta used need to be already simplifed to integer
				push @{ $contig_to_reads{"Contig $curr_contig"} }, "$1";
				#print $curr_contig."\t".$1."\n";
			}
			#skip when see "DETAILED DISPLAY OF CONTIGS" - may be used later
			if (/DETAILED DISPLAY OF CONTIGS/i){
				$cont_parse = 0;
				next;
			}
		}
	}
	return \%contig_to_reads;
}

sub map_2_genome {
	my $mapper = $_[0];
	my $seq_in = $_[1];
	my $reference = $_[2];
	my $thread = $_[3];
	my $out_sam = $_[4];
	my $logdir = $_[5];

	#my $dir = getcwd; #main dir
	#my $mapper = "$dir/thirdparty/bwa-0.6.1/bwa bwasw";
	my @cmd = split (/ /, "$mapper bwasw -t $thread $reference $seq_in -f $out_sam 2>>$logdir/map2genome.$seq_in.log");
	ExecCmd (\@cmd, $DEBUG_MODE);
}

sub blast_run {
	my $curr_process =(caller(0))[3];
	#print $curr_process."\n";
	my $usage = qq{Usage: $curr_process <fa in> <blast_db> <out file>
	};
	die($usage) if (@_!=7);

	#my $dir = getcwd; #main dir
	my $blast = $_[0];
	my $seq_file_in = $_[1];
	my $blast_db = $_[2];
	my $out = $_[3];
	my $evalue = $_[4];
	my $thread = $_[5];
	my $logdir = $_[6];

	#my $blast = "$dir/thirdparty/blast-2.2.26/bin/blastall";
	my @cmd = split (/ /, "$blast -p blastn -d $blast_db -m 7 -i $seq_file_in -a $thread -e $evalue -o $out 2>>$logdir/blastn.$seq_file_in.log"); #XML output
	ExecCmd (\@cmd, $DEBUG_MODE);
}

sub blast_formatdb { #This function is not needed because nt database are pre-formatted
	my $curr_process =(caller(0))[3];
	my $usage = qq{Usage: $curr_process <file in>
	};
	die($usage) if (@_!=1);
	my $dir = getcwd;
	my $file_in = $_[0];
	#my $formatdb = "$dir/thirdparty/blast-2.2.26/bin/formatdb";
	my $formatdb = `which formatdb`;
	my @cmd = split (/ /, "$formatdb -p F $file_in");
	ExecCmd (\@cmd, $DEBUG_MODE);
}

sub blast_parse {
	my $file_in = $_[0];
	open my $wFH, ">", "$file_in.parsed" or die "Can't open $file_in.parsed for writing : $!";
	print $wFH "query_name\thit_name\thit_description\tQuery_start\tQuery_end\tquery_aln_seq\tevalue\thsp_frac_identical\thsp_gap_count\tHit_start\tHit_end\tquery_length\thsp_length\thsp_over_query_len\tquery_strand\tmismatches pos on query(0-based)\n";
	my $in = new Bio::SearchIO(-format => 'blastxml',
							   -file   => $file_in);
	while( my $result = $in->next_result ) {
		while( my $hit = $result->next_hit ) {
			while( my $hsp = $hit->next_hsp ) {
				my $query_name = $result->query_name;
				my $query_length = $result->query_length;
				my $query_align_seq = $hsp->query_string;
				my $number_of_hits = $result->num_hits;
				my $hit_name = $hit->name;
				my $hit_length = $hit->length;
				my $evalue = $hit->significance;
				my $hit_description = $hit->description;
				my $hsp_length = $hsp->length('hit');
				my $hsp_over_query_len = $hsp_length/$query_length;
				my @seq_inds = $hsp->seq_inds('query','identical');
				my $query_strand = $hsp->strand('query');
				my $query_start_QUERY = $hsp->start('query');
				my $query_end_QUERY = $hsp->end('query');
				my $query_start_HIT = $hsp->start('hit');
				my $query_end_HIT = $hsp->end('hit');
				my $hsp_gap_count = $hsp->gaps;
				my $hsp_frac_identical = $hsp->frac_identical;
				#obtain mismatch pos.
				my @hsp_array = ($query_start_QUERY..$query_end_QUERY);
				my %hsp_pos = ();
				for (@seq_inds){
					$hsp_pos{$_} = "";
				}
				my @hsp_non_identical = ();
				for (@hsp_array){
					if (!defined $hsp_pos{$_}){
						my $act_pos = ($_ - $query_start_QUERY);
						#print "$act_pos\n";
						push (@hsp_non_identical, $act_pos);
					}
				}
				#mismatch pos to str
				#still need to decide if the position should be reverse if the query strand in blast is reverse (i.e. !1)
				my $MMapos = "";
				my $hspc = @hsp_non_identical;
				unless ($hspc == 0){
					$MMapos = join (";",@hsp_non_identical); # [single element case: 1] [more than 1: 1;2]
				}
				else {
					$MMapos = "NA";
				}
				print $wFH "$query_name\t$hit_name\t$hit_description\t$query_start_QUERY\t$query_end_QUERY\t$query_align_seq\t$evalue\t$hsp_frac_identical\t$hsp_gap_count\t$query_start_HIT\t$query_end_HIT\t$query_length\t$hsp_length\t$hsp_over_query_len\t$query_strand\t$MMapos\n";
			}
		}
	}
	close($wFH);
}

sub CS_read_level_analysis {
	my $curr_process =(caller(0))[3];
	print STDERR "Entering $curr_process.\n";
	my $usage = qq{Usage: $curr_process <runID> <MS BLAST result> <CS BLAST result> <MS fasta file> <CS fasta file> <CS_MS_ori_read_file> <Reference to clipped sequence keywords> <Reference to mapped segments keywords>
	};
	die($usage) if (@_<9);

	my ($runID,
	    $CS_RL_result_file,
	    $MS_BLAST_file,
	    $CS_BLAST_file,
	    $MS_fasta_file,
	    $CS_fasta_file,
	    $CS_MS_ori_read_file,
	    $clippedSeqKeywordsRef,
	    $mappedSeqKeywordsRef) = @_;

	############################
	##  Pre-defined variables ##
	############################
	my @CS_keywords;
	for my $keyword (@$clippedSeqKeywordsRef) {
		print STDERR "CS keys: $keyword\n";
		push (@CS_keywords, $keyword);
	}
	my @MS_keywords;
	for my $keyword (@$mappedSeqKeywordsRef) {
		print STDERR "MS keys: $keyword\n";
		push (@MS_keywords, $keyword);
	}

	#Internal defined threshold for specificity
	my $first_pass_cov_threshold = 0.80;
	my $CS_reporting_threshold = 0.80;

	#Define CSmethod output file
	open my $CSoutFH, '>' , $CS_RL_result_file or die "Can't open $CS_RL_result_file for writing : $!";
	my @MS = file_2_array_skip_header($MS_BLAST_file);
	my @CS = file_2_array_skip_header($CS_BLAST_file);
	my %MSs = %{blast_to_data_str(\@MS)};
	my %CSs = %{blast_to_data_str(\@CS)};
	my @MSreadsID = (); my %unify =();
	for (@MS){
		my @s = split(/\t/,$_);
		$unify{$s[0]} = "";
	}
	for my $k (keys %unify){ #reduce readIDs to non-redundant form
		push (@MSreadsID, $k);
	}
	# May have to check if readsID in MS and CS are identical (to make sure every previoius steps finsihed correctly). Should run find without checking, but is good to check. NOT PRIORITY TASK

	############################
	#MS parsing
	############################
	my $total_analyzed_read_c = 0;
	my @MS_optimal_outcome = ();
	my $MS_optimal_outcome_v;
	for my $readID (@MSreadsID){
		my $MS_analyzed_read_c = 0; #the sum of blast result of each sequence read
		$total_analyzed_read_c++; #the sum of blast results in MS file
		#print "analyzing $readID\n";
		my @MSa = @{$MSs{$readID}};
		my ($MS_total_qualified_q, $MS_q_status, $MS_q_positive, $MS_q_negative) = (0,0,0,0);
		my %MMEF_MS = ();
		##MS specificty check
		for (@MSa){
			$MS_q_status = 0;
			my @s = split(/\t/, $_);
			if ($s[-3]>$first_pass_cov_threshold){
				for my $q_word (@MS_keywords){
					if ($s[1] =~ /$q_word/gi){
						$MS_q_status = 1;
					}
				}
				if ($MS_q_status == 1){
					$MS_q_positive++;
				}
				else {
					$MS_q_negative++;
				}
				$MS_total_qualified_q++;
			}
			$MS_analyzed_read_c++;
			##Prepare for MMEF calculation
			#s[3] is query alignment length
			#s[-1] contians the MS positions
			#get number of mismatches
			my $MM_c = 0;
			#infer MM from ";"
			unless ($s[-1] =~/NA/gi){
				if ($s[-1] =~/,/){
					$MM_c = $s[-1] =~ tr/;//;
				}
				else {
					$MM_c = 1;
				}
			}
			my $MMEF_MS_score = abs($s[-4] - $MM_c);
			push (@{$MMEF_MS{$MMEF_MS_score}}, $s[1]);
			#push (@{$MMEF_MS{$MMEF_MS_score}}, $_);
		}
		my ($MS_specific_percent_positive, $MS_negative_match_percent);
		if ($MS_total_qualified_q > 0){
			$MS_specific_percent_positive = ($MS_q_positive / $MS_total_qualified_q)*100;
			$MS_negative_match_percent = ($MS_q_negative / $MS_total_qualified_q)*100;
		}
		else {
			$MS_specific_percent_positive = "NA";
			$MS_negative_match_percent = "NA";
		}
		if ($MS_specific_percent_positive =~/^\d+$/){
			if ($MS_specific_percent_positive < 95){
				#high quality non target match
			}
			else {
				#High quality specific match to target
			}
		}
		else {
			#the sequence has no good match (defined by coverage) to any target in NT
		}
		#print "$readID\t$MS_q_positive\t$MS_analyzed_read_c\t$MS_total_qualified_q\t$MS_specific_percent_positive\t$MS_negative_match_percent\n";
		#the more equally good matches to MS target is good, it means the sequence originated from a conserved part of target species
		my @sorted_MMEF_MS_k = sort {$b <=> $a} keys %MMEF_MS; #max key value
		#print "@{$MMEF_MS{$sorted_MMEF_MS_k[0]}}\n"; exit;
		@MS_optimal_outcome = @{$MMEF_MS{$sorted_MMEF_MS_k[0]}};
		$MS_optimal_outcome_v = $sorted_MMEF_MS_k[0];

		############################
		#Parse the MS associated CS
		############################
		if ($CSs{$readID}){
			my @CSa = @{$CSs{$readID}};
			my @CS_optimal_outcome = ();
			my $CS_optimal_outcome_v;
			#get best CS part, highest coverage, longest alignment length, minimal mismatch
			my ($CS_total_qualified_q, $CS_q_status, $CS_q_positive, $CS_q_negative) = (0,0,0,0);
			my %MMEF_CS = ();
			##CS specificity check
			#exclude those with key words in @MS, record #non-target matches (not desired)
			for (@CSa){
				$CS_q_status = 0;
				my @s = split (/\t/, $_);
				if ($s[-3]>$first_pass_cov_threshold){
					for my $q_word (@MS_keywords){
						if ($s[1] =~ /$q_word/gi){
							$CS_q_status = 1;
						}
					}
					if ($CS_q_status == 1){
						$CS_q_negative++;
					}
					else {
						for my $q_word (@CS_keywords){
							if ($s[1] =~ /$q_word/gi){
								$CS_q_status = 1;
							}
						}
						if ($CS_q_status == 1){
							$CS_q_positive++;
						}
					}
					$CS_total_qualified_q++;
				}
				##Prepare for MMEF calculation
				#s[3] is query alignment length
				#s[-1] contians the MS positions
				#get number of mismatches
				my $MM_c = 0;
				#infer MM from ";"
				unless ($s[-1] =~/NA/gi){
					if ($s[-1] =~/,/){
						$MM_c = $s[-1] =~ tr/;//;
					}
					else {
						$MM_c = 1;
					}
				}
				my $MMEF_CS_score = abs($s[-4] - $MM_c);
				push (@{$MMEF_CS{$MMEF_CS_score}}, $s[1]);
				#push (@{$MMEF_CS{$MMEF_CS_score}}, $_);
			}
			##Obtain best CS
			my @sorted_MMEF_CS_k = sort {$b <=> $a} keys %MMEF_CS; #max key value
			#print "@{$MMEF_MS{$sorted_MMEF_MS_k[0]}}\n"; exit;
			@CS_optimal_outcome = @{$MMEF_CS{$sorted_MMEF_CS_k[0]}};
			$CS_optimal_outcome_v = $sorted_MMEF_CS_k[0];
			##Claculate combined MMEF
			my $MMEF = min($MS_optimal_outcome_v,$CS_optimal_outcome_v);
			##Hash MS&CS sequences previously used for BLAST
			my %MS_readID_seq = %{hash_fasta_ID_2_seq($MS_fasta_file)};
			my %CS_readID_seq = %{hash_fasta_ID_2_seq($CS_fasta_file)};
			my %MS_CS_ori_readID_seq = %{hash_fasta_ID_2_seq($CS_MS_ori_read_file)};
			##Report result
			if ($CS_total_qualified_q > 0){
				unless (($CS_q_positive/$CS_total_qualified_q) < $CS_reporting_threshold){
					print $CSoutFH "$readID\t$MS_q_positive\t$MS_analyzed_read_c\t$MS_total_qualified_q\t$MS_specific_percent_positive\t$MS_negative_match_percent\t";
					print $CSoutFH "$CS_total_qualified_q\t$CS_q_positive\t$CS_q_negative\t";
					print $CSoutFH "$MMEF\t";
					my ($MSseq, $CSseq);
					if (defined $MS_readID_seq{$readID}){
						$MSseq = join ("", @{$MS_readID_seq{$readID}});
					} else {
						print "MS for $readID not defined\n";
						next;
					}
					if (defined $CS_readID_seq{$readID}){
						$CSseq = join ("", @{$CS_readID_seq{$readID}});
					} else {
						print "CS for $readID not defined\n";
						next;
					}
					my @split_rname = split (/\_/, $readID);
					print "$split_rname[0]\n";
					my $MS_CS_ori_seq = join ("", @{$MS_CS_ori_readID_seq{$split_rname[0]}});
					print $CSoutFH "$MS_optimal_outcome_v\t$CS_optimal_outcome_v\t$MSseq\t$CSseq\t$MS_CS_ori_seq\t";
					#print "@{$MMEF_MS{$sorted_MMEF_MS_k[0]}}\t";
					my $MS_str = join " ",@{$MMEF_MS{$sorted_MMEF_MS_k[0]}};
					#print "$MS_str\n";
					my $MS_out = reduce_str_complexity($MS_str);
					print $CSoutFH $MS_out."\t";
					my $CS_str = join " ",@{$MMEF_CS{$sorted_MMEF_CS_k[0]}};
					#print "@{$MMEF_CS{$sorted_MMEF_CS_k[0]}}\n";
					my $CS_out = reduce_str_complexity($CS_str);
					print $CSoutFH $CS_out."\n";
					#print "\n\n";
					##Stich the sequence back to read <-> Find Viral-Human or Human-Viral
				}
			}
		}
	}
	close ($CSoutFH);
}

sub Assemble_CS {
	my $curr_process =(caller(0))[3];
	print STDERR "Entering $curr_process.\n";
	my $usage = qq{Usage: $curr_process <aligner_runID> <cap3 bin> <clipped seq result file>
	};
	die($usage) if (@_!=4);
	my $aligner_runID = $_[0];
	my $cap3 = $_[1];
	my $CS_RL_result_file = $_[2]; #The result of read-level analysis, where highly confident clipped-seq are reported. We use the CS & MS parts there to reduce noises in de novo assembly
	my $logdir = $_[3];

	##Assemble the MS and CS part of the reads reported by the Clipped-Seq analysis module
	my %MS_CS_contigs_to_anonymous_array_of_simIDlist_simID_hash_to_readID = ();
	my %MS_CS_simID_to_readID = ();
	my $CS_RL_MS_fasta_file = "$CS_RL_result_file.MS.fa";
	my $CS_RL_CS_fasta_file = "$CS_RL_result_file.CS.fa";
	open my $CS_RL_MS_fasta_FH, '>', $CS_RL_MS_fasta_file or die "Can't open $CS_RL_MS_fasta_file for writing : $!";
	open my $CS_RL_CS_fasta_FH, '>', $CS_RL_CS_fasta_file or die "Can't open $CS_RL_CS_fasta_file for writing : $!";
	my @CS_RL_fasta_a = @{file_2_array($CS_RL_result_file)};
	for my $l(@CS_RL_fasta_a){
		my @s_CS_RL = split (/\t/, $l);
		print $CS_RL_MS_fasta_FH '>'.$s_CS_RL[0]."\n".$s_CS_RL[12]."\n";
		print $CS_RL_CS_fasta_FH '>'.$s_CS_RL[0]."\n".$s_CS_RL[13]."\n";
	}
	close ($CS_RL_MS_fasta_FH);
	close ($CS_RL_CS_fasta_FH);
	##MS
	my %readsID_simplified_tracking_MS = %{create_simplified_readsID_fa("$CS_RL_MS_fasta_file", "$CS_RL_MS_fasta_file.sim", "$CS_RL_MS_fasta_file.tracking")};
	my %MSreads_in_contigs = %{cap3_run("$CS_RL_MS_fasta_file.sim", $cap3, $logdir)};
	for my $ra ( keys %MSreads_in_contigs ) { #$ra is the Contig #x, $MSreads_in_contigs hashs the simplied, numerical values
		my @a = @{$MSreads_in_contigs{$ra}};
		for (@a){
			#print "$ra\t$_\t$readsID_simplified_tracking_MS{$_}\n";
			push(@{$MS_CS_contigs_to_anonymous_array_of_simIDlist_simID_hash_to_readID{'MS'}{$ra}}, $_);
			$MS_CS_simID_to_readID{'MS'}{$_}=$readsID_simplified_tracking_MS{$_};
		}
	}
	for my $key ( keys %readsID_simplified_tracking_MS ) { #key is the simplified, numerical values, readsID_simplified_tracking_MS hashes the read name
		#print "$key\t$readsID_simplified_tracking_MS{$key}\n";
	}
	##CS
	my %readsID_simplified_tracking_CS = %{create_simplified_readsID_fa("$CS_RL_CS_fasta_file", "$CS_RL_CS_fasta_file.sim", "$CS_RL_CS_fasta_file.tracking")};
	my %CSreads_in_contigs = %{cap3_run("$CS_RL_CS_fasta_file.sim", $cap3, $logdir)};
	for my $ra ( keys %CSreads_in_contigs ) {
		my @a = @{$CSreads_in_contigs{$ra}};
		for (@a){
			#print "$ra\t$_\n";
			push(@{$MS_CS_contigs_to_anonymous_array_of_simIDlist_simID_hash_to_readID{'CS'}{$ra}}, $_);
			$MS_CS_simID_to_readID{'CS'}{$_}=$readsID_simplified_tracking_CS{$_};
		}
	}
	for my $key ( keys %readsID_simplified_tracking_CS ) {
		#print "$key\t$readsID_simplified_tracking_CS{$key}\n";
	}
	## Return a hash of contigs->readsID
	my @a = (\%MS_CS_contigs_to_anonymous_array_of_simIDlist_simID_hash_to_readID, \%MS_CS_simID_to_readID);
	return \@a;
}

sub tx_targeted_assembly {
	#takes MS, CS result file from CS_read_level_analysis
	my $curr_process =(caller(0))[3];
	print STDERR "Entering $curr_process.\n";
	my $usage = qq{Usage: $curr_process <aligner_runID> <clipped seq result file> <Read-pair method result file> <Forward read> <Reverse Read> <insert size> <sensitive tx denovo 1/0> <logdir>
	};
	die($usage) if (@_ < 13);
	my $aligner_runID = $_[0];
	my $CS_RL_result_file = $_[1];
	my $CS_MS_ori_read_file = $_[2];
	my $RPm_result_file = $_[3];
	my $F_read_file = $_[4];
	my $R_read_file = $_[5];
	my $human_ref = $_[6];
	my $humanDecoy = $_[7];
	my $threads = $_[8];
	my $bwa_bin = $_[9];
	my $SSAKE_bin = $_[10];
	my $samtools = $_[11];
	my $isize_supplied = $_[12];
	my $sensitive_tx_denovo_assembly = $_[13];
	my $logdir = $_[14];

	my $range_boundary_threshold = 500; ## Proximity to RPm hits (bp): Parameter used for extraction of read pair for senstiive denovo
	########################
	#Temporary files
	########################
	my $CS_RL_MS_CS_fullread_fasta_file = "$CS_RL_result_file.fullread.fa";
	##CS_RL to temp fasta file
	my @CS_RL_fasta_a = @{file_2_array($CS_RL_result_file)};
	open my $CS_RL_MS_CS_fullread_fasta_file_FH, '>', $CS_RL_MS_CS_fullread_fasta_file or die "Can't open $CS_RL_MS_CS_fullread_fasta_file for writing : $!";
	for my $l(@CS_RL_fasta_a){
		my @s_CS_RL = split (/\t/, $l);
		print $CS_RL_MS_CS_fullread_fasta_file_FH '>'.$s_CS_RL[0]."\n".$s_CS_RL[14]."\n";
	}
	close ($CS_RL_MS_CS_fullread_fasta_file_FH);
	##RP
	my %RP_vicinity_reads = ();
	my $RPm_F_fastq_out = $RPm_result_file."_1.fa";
	my $RPm_R_fastq_out = $RPm_result_file."_2.fa";
	open my $F_fastq_out_FH, '>', $RPm_F_fastq_out or die "Can't open $RPm_F_fastq_out for writing : $!";
	open my $R_fastq_out_FH, '>', $RPm_R_fastq_out or die "Can't open $RPm_R_fastq_out for writing : $!";
	my @RPm_fasta_a = @{file_2_array($RPm_result_file)}; # have to reduce readID
	for (@RPm_fasta_a){
		my @s_RPm = split (/\t/, $_);
		my $v_readID = $s_RPm[3];
		my $h_chr = $s_RPm[-20];
		my $h_read_st = $s_RPm[-19];
		my $h_read_ed = $s_RPm[-18];
		my $range_minus = ($h_read_st - $range_boundary_threshold);
		my $range_plus = ($h_read_ed + $range_boundary_threshold);
		my @range_in = ($range_minus..$range_plus);
		push (@{$RP_vicinity_reads{$h_chr}}, \@range_in);
		my $h_readID = $s_RPm[-17];
		if ($v_readID =~ /\/1/){
			print $F_fastq_out_FH '>'.$v_readID."\n".$s_RPm[-21]."\n";
		} else {
			print $R_fastq_out_FH '>'.$v_readID."\n".$s_RPm[-21]."\n";
		}
		if ($h_readID =~ /\/1/){
			print $F_fastq_out_FH '>'.$h_readID."\n".$s_RPm[-1]."\n";
		} else {
			print $R_fastq_out_FH '>'.$h_readID."\n".$s_RPm[-1]."\n";
		}
	}
	close ($F_fastq_out_FH);
	close ($R_fastq_out_FH);
	##SSAKE
	if ($isize_supplied == 0) {
		print STDERR "No insert size supplied by user. VFS is going to find insert size\n";
		$isize_supplied = bwa_pipeline($human_ref, "$aligner_runID.human.isize", $threads, $F_read_file, $R_read_file, $bwa_bin, 1); #<reference genome> <SAM prefix> <threads> <Forward fastq> <Reverse fastq> <obtain IS 1/0>
		print STDERR "Estimated Insert size: $isize_supplied\n";
	} else {
		print STDERR "Skipped Insert Size estimation. Supplied insert size: $isize_supplied\n";
	}
	#ssake_reads_prep("$aligner_runID.RPm.out_1.fa","$aligner_runID.RPm.out_2.fa", "fasta", $isize_supplied, "$aligner_runID.RPm.ssake.in");
	ssake_reads_prep($RPm_F_fastq_out, $RPm_R_fastq_out, "fasta", $isize_supplied, "$aligner_runID.RPm.ssake.in");

	unless ($sensitive_tx_denovo_assembly  == 1){
		ssake_run($SSAKE_bin, "$aligner_runID.RPm.ssake.in", $CS_RL_MS_CS_fullread_fasta_file, "$aligner_runID.targeted.assembly", $logdir); # <SSAKE bin> <fa in> <seed in> <fileout prefix>
	}
	else {
	    ###################################
	    ## Sensitive tx denovo
	    ###################################
		#return CS_min_mates_seq ($CS_MS_ori_read_file.mates.fa)
		ssake_reads_prep($CS_MS_ori_read_file, "$CS_MS_ori_read_file.mates.fa", "fasta", $isize_supplied, "$aligner_runID.CS_rp.ssake.in");

		#return Readpair with one end mapped in close proximity to RPm human part, but antoher end unmapped.
		bwa_pipeline($humanDecoy, "$aligner_runID.human_targeted", $threads, $F_read_file, $R_read_file, $bwa_bin, 0);
		##parse $vfs_prefix.$runID.human_targeted.sorted.bam
		my $need_parse = 1;
		if ($need_parse == 1){ ##For development use only
			my $bam_file = "$aligner_runID.human_targeted.sorted.bam";
			open my $BAMrFH, " $samtools view $bam_file |" or die "Can't open pipe : $!";
			print STDERR "1st pass on $bam_file\n";
			my %vicinity_rp = ();
			while (<$BAMrFH>){
				chomp;
				my @s = split (/\t/, $_);
				my ($bitwise_flag, $curr_read_chr, $curr_read_st, $curr_read_len) = ($s[1], $s[2], $s[3], $s[9]);
				my $len = length ($curr_read_len);
				my $curr_read_ed = ($curr_read_st + $len);
				unless (($bitwise_flag & 0x4) == 0x4){ #this read IS mapped
					if (($bitwise_flag & 0x8) == 0x8){ #another end unmapped
						if (defined $RP_vicinity_reads{$curr_read_chr}){
							my @match_range_ref_in_a = @{$RP_vicinity_reads{$curr_read_chr}};
							for my $match_range_ref (@match_range_ref_in_a){
								my @match_range = @{$match_range_ref};
								if (in_range(\@match_range, $curr_read_st, $curr_read_ed)){
									my $str_in = '@'.$s[0];
									$vicinity_rp{$str_in} = "";	#readname core
								}
							}
						}
					}
				}
			}
			close ($BAMrFH);

			#open my $FfastqFH, '<', $F_read_file; #change to support gz
			open my $FfastqFHwrote, '>', "$F_read_file.vicinity";
			my $is_gzip_file = DetermineFileType ($F_read_file);
			my $FfastqFH = OpenFastqFileRead ($F_read_file, $is_gzip_file);
			print STDERR "Extracting RP vicinity ($range_boundary_threshold bp) $F_read_file from $bam_file\n";
			my $line; my $seq; my $aux; my $qual;
			#while (1){

			while (<$FfastqFH>){
				my @split = split (/ /, $_);
				$seq = <$FfastqFH>; chomp $seq;
				$aux = <$FfastqFH>; chomp $aux;
				my $qual = <$FfastqFH>; chomp $qual;
				if (defined $vicinity_rp{$split[0]}){
					print $FfastqFHwrote $split[0]."\n".$seq."\n".$aux."\n".$qual."\n";
				}
			}

			#}
			#close ($FfastqFH);
			CloseFastqFile ($FfastqFH);
			close ($FfastqFHwrote);


			#open my $RfastqFH, '<', $F_read_file; #change to support gz
			open my $RfastqFHwrote, '>', "$R_read_file.vicinity";
			$is_gzip_file = DetermineFileType ($R_read_file);
			my $RfastqFH = OpenFastqFileRead ($R_read_file, $is_gzip_file);
			print STDERR "Extracting RP vicinity ($range_boundary_threshold bp) $R_read_file from $bam_file\n";
			#while (1){

			while (<$RfastqFH>){
				my @split = split (/ /, $_);
				$seq = <$RfastqFH>; chomp $seq;
				$aux = <$RfastqFH>; chomp $aux;
				$qual = <$RfastqFH>; chomp $qual;
				if (defined $vicinity_rp{$split[0]}){
					print $RfastqFHwrote $split[0]."\n".$seq."\n".$aux."\n".$qual."\n";
				}
			}
				#close ($RfastqFH);
			CloseFastqFile ($RfastqFH);
			close ($RfastqFHwrote);
			#}
		}
		ssake_reads_prep("$F_read_file.vicinity", "$R_read_file.vicinity", "fastq", $isize_supplied, "$aligner_runID.RPm.vicinity.ssake.in");
		##merge to super-sensitive.ssake_in
		system ("cat $aligner_runID.CS_rp.ssake.in $aligner_runID.RPm.ssake.in $aligner_runID.RPm.vicinity.ssake.in >$aligner_runID.RPm.plus.vicinity.CS_FnR.ssake.in");
	}
	ssake_run($SSAKE_bin, "$aligner_runID.RPm.plus.vicinity.CS_FnR.ssake.in", $CS_RL_MS_CS_fullread_fasta_file, "$aligner_runID.targeted.assembly.sensitive", $logdir); # <SSAKE bin> <fa in> <seed in> <fileout prefix>
}

sub hash_fasta_ID_2_seq {
	my $file_in = $_[0];
	open my $rFH, '<', $file_in or die "Can't open $file_in for reading : $!";
	my $readID; my $seq;
	my %h = ();
	while (<$rFH>){
		chomp;
		if (/^>/){
			$readID = substr $_, 1;
		}
		else {
			push (@{$h{$readID}}, $_);
		}
	}
	return \%h;
}

sub reduce_str_complexity {
	my @words=  split(" ", $_[0]);
	my %done = ();
	my @newwords = ();
	for (@words) {
		push(@newwords,$_) unless $done{$_}++;
	}
	my $str = join ' ', @newwords;
	return $str;
}

sub max {
	my ($max, $next, @vars) = @_;
	return $max if not $next;
	return max( $max > $next ? $max : $next, @vars );
}

sub min {
    my ($min, $next, @vars) = @_;
    return $min if not $next;
    return min( $min < $next ? $min : $next, @vars );
}

sub blast_to_data_str {
	my @a = @{$_[0]};
	my %h = ();
	for (@a){
		my @s = split (/\t/, $_);
		my $readID = $s[0];
		shift (@s);
		my $str = join ("\t", @s);
		push (@{$h{$readID}}, $str);
	}
	return \%h;
}



sub ssake_reads_prep {
	my $curr_process =(caller(0))[3];
	print $curr_process."\n";
	my $usage = qq{Usage: $curr_process <F> <R> <fasta / fastq> <insert size> <out file>
	};
	die($usage) if (@_!=5);
	my $F_file = $_[0];
	my $R_file = $_[1];
	my $file_type = $_[2];
	my $isize = $_[3];
	my $out_file = $_[4];
	my $fastq = 0; #Default assume fasta
	if ($file_type =~ /q/g){
		$fastq = 1;
	}
	#printf "RP_1: -$F_file- [%u]\n", (-s $F_file);
	#printf "RP_2: -$R_file- [%u]\n", (-s $F_file);
	#print "File type as specified: $file_type\nOutput file: $out_file\nIf FastQ: $fastq\n";
	#It's feasible to hash entire input file because the input files are RPm output (i.e. not original fastq)
	if ($fastq == 1){ #illumina 1.8 fastq
		print STDERR "Assuming fastq file\n";
		#hash R read
		my %Rreadh = ();
		my $i = 2;
		print STDERR "Hashing $R_file\n";
		open my $Rfh, '<', $R_file or die "Can't open $R_file for reading : $!";
		while (<$Rfh>){
			chomp;
			my @s = split (/ /, $_);
			my $readname = $s[0];
			$readname =~ s/://g;
			my $seq = <$Rfh>; chomp $seq;
			my $aux = <$Rfh>; chomp $aux;
			my $qual = <$Rfh>; chomp $aux;
			$Rreadh{$readname} = $seq;
		}
		close($Rfh);
		#process the reads
		open my $Lfh, '<', $F_file || die "Can't open $F_file for reading : $!";
		open my $wFH, '>', $out_file or die "Can't open $out_file for writing : $!";
		print STDERR "Processing $F_file and $R_file\n";
		my $q_readname; my $R_seq;
		while (<$Lfh>){
			chomp;
			if (/^@/){
				my @s = split (/ /, $_);
				$q_readname = $s[0];
				$q_readname =~ s/://g;
				$R_seq = $Rreadh{$q_readname};
			}
			else {
				unless (($_ =~ /n/ig)||($R_seq =~ /n/ig)){
					print $wFH ">$q_readname:$isize\n";
					print $wFH "$_:$R_seq\n";
				}
				my $aux = <$Lfh>; chomp $aux;
				my $qual = <$Lfh>; chomp $aux;
			}
		}
	} else { #fasta
		print STDERR "Assuming fasta file\n";
		#hash R read
		my %Rreadh = ();
		printf STDERR "Hashing $R_file\n";
		open my $Rfh, '<', $R_file || die "Can't open $R_file : $!";;
		while (<$Rfh>){
			chomp;
			if (/^>/){
				my @s = split (/\//, $_);
				my $readname = $s[0];
				if ($readname =~/:/){
					$readname =~ s/://g;
				}
				my $seq = <$Rfh>; chomp $seq;
				$Rreadh{$readname} = $seq;
			}
		}
		close($Rfh);
		open my $Lfh, '<', $F_file or die "Can't open $F_file for reading : $!";
		open my $wFH, '>', $out_file or die "Can't open $out_file for writing : $!";
		print STDERR "Processing $F_file and $R_file\n"; sleep 4;
		my $q_readname; my $R_seq;
		while (<$Lfh>){
			chomp;
			if (/^>/){
				my @s = split (/\//, $_);

				if ($s[0] =~ /_/){
					my @ss = split (/_/, $s[0]);
					$s[0] = $ss[0];
				}

				$q_readname = $s[0];
				if ($q_readname =~/:/){
					$q_readname =~ s/://g;
				}
				$R_seq = $Rreadh{$q_readname};
			}
			else {
				unless (($_ =~ /n/ig)||($R_seq =~ /n/ig)){ #SSAKE only takes [ATGC]
					print $wFH "$q_readname:$isize\n";
					print $wFH "$_:$R_seq\n";
				}
			}
		}
	}
}

sub ssake_run {
	my $curr_process =(caller(0))[3];
	my $usage = qq{Usage: $curr_process <fa in> <seed in> <fileout prefix>
	};
	die($usage) if (@_!=5);
	#my $dir = getcwd;
	my $ssake = $_[0];
	my $fa_in = $_[1];
	my $seed_in = $_[2];
	my $out_prefix = $_[3];
	my $logdir = $_[4];

	my @cmd = split (/ /, "$ssake -f $fa_in -w 1 -s $seed_in -p 1 -o 3 -r 0.9 -b $out_prefix 2>>$logdir/ssake.log"); #-b denotes basename of output file
	if ($DEBUG_MODE) {
		print join (" ", @cmd), "\n";
	}
	ExecCmd (\@cmd, $DEBUG_MODE);
}

sub clean_tempfiles {
	my $curr_process =(caller(0))[3];
	my $dir_prefix = $_[0];
	my $cmd;
	my $curr_dir = getcwd;
	#print "$dir\n";
	my $temp_dir = $dir_prefix."_temp";
	if (-d $temp_dir){ #exists
		if (-w $temp_dir) { #writable
			system("mv $dir_prefix.* $temp_dir/");
		}
	} else {
		mkdir ($temp_dir);
		system("mv $dir_prefix.* $temp_dir/");
	}
	##Move result files to current root directory
	system ("mv $temp_dir/*.RPm.out $curr_dir/");
	system ("mv $temp_dir/*.CSm.out $curr_dir/");
	system ("mv $temp_dir/*.targeted.*.contigs $curr_dir/");
}

 sub file_2_array_skip_header {
 	my $curr_process =(caller(0))[3];
 	my $usage = qq{Usage: $curr_process <file in>
	};
	die($usage) if (@_!=1);
	my $file_in = $_[0];
	print STDERR "Reading $file_in\n";
	my @out = ();
	my $i = 1;
	open my $fh, '<', $file_in or die "Can't open $file_in for reading : $!";
	while (<$fh>){
		chomp;
		unless ($i == 1){
			push (@out, $_);
		}
		$i = 0;
	}
	close ($fh);
	return @out;
}

sub file_2_array {
	my $file_in = $_[0];
	my @out = ();
	open my $file2_a_fh, '<', $file_in or die "Can't open $file_in for reading : $!";
	while (<$file2_a_fh>){
		chomp;
		push (@out, $_);
	}
	close ($file2_a_fh);
	return \@out;
}

sub rtn_fastq_phred {
	my $curr_process =(caller(0))[3];
	#my $dir_prefix = $_[0];
	my @in_a = @{$_[0]};
	my $sampling_size = 100000;			# sample size
	my $cnt;
	my $default_phred = 64; #default phred encoding: 64
	my $id;
	my $seq;
	my $comm;
	my $qual;
	my @rtn = ();
	my $is_gzip_input = 0;

	#if (@_ > 0) {
	if (@in_a > 0) {
		#for my $f (@_) {
		for my $fn (@in_a) {
			my $FH;
			$is_gzip_input = DetermineFileType ($fn);

			if ($is_gzip_input) {
				open ($FH, "<:gzip", $fn) or die "Can't open $fn (gzipped) for reading : $!";
			}
			else {
				open $FH, '<', $fn or die "Can't open $fn for reading : $!";
			}
			$id = <$FH>;

			if (($id =~ /^\@/) && ($id !~ /^\@SQ\t/)) { #ID
				scalar <$FH>; #seq
				$comm = scalar <$FH>; #aux
				$qual = <$FH>; #qual
				chomp $qual;
			}
			while($qual) {
				for (my $i =length($qual)/2; $i < length($qual); ++$i) {
					if (ord(substr($qual,$i,1)) < 64) {
						$default_phred = 33;
						$cnt=$sampling_size;	# last
						last;
					}
				}
				$qual = '';
				last if ++$cnt >= $sampling_size;
				# fastq
				scalar <$FH>;		# id
				scalar <$FH>;		# read
				scalar <$FH>;		# comment
				$qual = <$FH>;
				chomp $qual;
			}
			push (@rtn, "$fn\t$default_phred");

			close ($FH);
		}
	}

	return \@rtn;
}

1;


