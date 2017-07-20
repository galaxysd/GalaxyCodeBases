#!/usr/bin/perl
###########################################################################
##  ViralFusionSeq
##  Software for discovering and annotating viral integration event and
##    fusion transcript
##  
##  Version 2.0
##  
##  Copyright (C) 2016 by Jing-Woei Li & Raymond Wan, All rights reserved.
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

##  $LastChangedDate: 16th Aug 2016 $
##  $LastChangedRevision: 1289 $

package RPmethod;

use diagnostics;
use strict;
use warnings;

$| = 1; #disable buffer

use Exporter;
use AppConfig;
use PerlIO::gzip;

use FindBin;
use lib "$FindBin::Bin/include";
use BWAauto; #BWA sampe/samse pipeline
use General; # Module for general subroutines

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = qw(exec_RPmethod);
our @EXPORT_OK   = qw(exec_RPmethod);
our %EXPORT_TAGS = ( DEFAULT => [qw(&exec_RPmethod)]);

sub exec_RPmethod {
	my $curr_process =(caller(0))[3];
	print STDERR "executing $curr_process.\n";
	my $vfs_prefix = $_[0];
	my $runID = $_[1];
	my $RPm_result_file = $_[2];
	my $threads = $_[3];
	my $forward_fastq = $_[4];
	my $reverse_fastq = $_[5];
	my $human_ref_fa = $_[6];
	my $viral_ref_fa_prefix = $_[7];
	my $bwa_bin = $_[8];
	my $samtools = $_[9];
	my $BEDtoolPATH = $_[10];
	my $orfRef = $_[11];
	my $logdir = $runID;

	##  Set up the ORF arrays
	my %orfs;
	my %orf_status;

	for (my $i = 0; $i < CountORFs ($orfRef); $i++) {
		my $orf_name = GetNameORFs ($i, $orfRef);
		my $orf_start = GetStartORFs ($i, $orfRef);
		my $orf_end = GetEndORFs ($i, $orfRef);

		##  Initialize the ORF range
		print "ORF: $orf_name\t$orf_start\t$orf_end\n";
		my @a_in = ($orf_start..$orf_end);
		#$orfs{$orf_name} = ($orf_start..$orf_end);
		$orfs{$orf_name} = \@a_in;

		##  Initialize the status to false
		$orf_status{$orf_name} = 0;
	}

	### pre-defined ###
	my $geneBED = "$vfs_prefix.$runID.orphan.human.sorted.mapped.bed.gene.anno";
	my $repeatBED = "$vfs_prefix.$runID.orphan.human.sorted.mapped.bed.repeat.anno";
	my $NM_BED = "$vfs_prefix.$runID.orphan.human.sorted.mapped.NM.bed";
	my $original_BED = "$vfs_prefix.$runID.orphan.human.sorted.mapped.bed";
	my $viral_BED = "$vfs_prefix.$runID.viral.mapped.bed";
	my $viral_NM_BED = "$vfs_prefix.$runID.viral.mapped.NM.bed";

	my @suffix = ("amb", "ann", "bwt", "pac", "sa");
	my $p = 1;
	for (@suffix){
		my $ck = $human_ref_fa.".".$_;
		unless ((-e $ck) && (-s $ck)){
			warn "$ck is either non-exist or its file size is 0 byte";
			$p = 0;
		}
	}
	for (@suffix){
		my $ck = $viral_ref_fa_prefix.".".$_;
		unless ((-e $ck) && (-s $ck)){
			warn "$ck is either non-exist or its file size is 0 byte";
			$p = 0;
		}
	}
	die ("Please make sure all references files are present in the directory") if ($p != 1);

	# (1) Paired-end mapping to viral decoy
	bwa_pipeline($viral_ref_fa_prefix, "$vfs_prefix.$runID.viral", $threads, $forward_fastq, $reverse_fastq, $bwa_bin, 0, $samtools);

	# (2) Extract one end mapped reads, extract another end
	## readname and sequence
	my $BWA_viral_bam_file = "$vfs_prefix.$runID.viral.sorted.bam";
	#hash sequence reads
	my %fastqH = %{reads_factory($samtools, $BWA_viral_bam_file, "$vfs_prefix.$runID", $forward_fastq, $reverse_fastq, $viral_ref_fa_prefix)};

	# (3) map onto human decoy sequences
	bwa_pipeline($human_ref_fa, "$vfs_prefix.$runID.orphan.human", $threads, "$vfs_prefix.$runID.orphan_F.fq", "$vfs_prefix.$runID.orphan_R.fq", $bwa_bin, 0, $samtools);
	my @cmd = split (/ /, "$samtools view -F 0x4 $vfs_prefix.$runID.orphan.human.sorted.bam | $samtools view -Sb -T $human_ref_fa - >$vfs_prefix.$runID.orphan.human.sorted.mapped.bam");
	ExecCmd (\@cmd, $DEBUG_MODE);
	#annotation step
	#BAMtoBED
	#Viral
	@cmd = split (/ /, "$BEDtoolPATH/bamToBed -cigar -i $vfs_prefix.$runID.viral.sorted.bam >$viral_BED");
	ExecCmd (\@cmd, $DEBUG_MODE);
	@cmd = split (/ /, "$BEDtoolPATH/bamToBed -ed -i $vfs_prefix.$runID.viral.sorted.bam >$viral_NM_BED");
	ExecCmd (\@cmd, $DEBUG_MODE);
	#Human
	@cmd = split (/ /, "$BEDtoolPATH/bamToBed -cigar -i $vfs_prefix.$runID.orphan.human.sorted.mapped.bam >$original_BED");
	ExecCmd (\@cmd, $DEBUG_MODE);
	@cmd = split (/ /, "$BEDtoolPATH/bamToBed -ed -i $vfs_prefix.$runID.orphan.human.sorted.mapped.bam >$NM_BED");
	ExecCmd (\@cmd, $DEBUG_MODE);
	#add_aux_human("$vfs_prefix.$runID.orphan.human.sorted.mapped.bed", $repeatBED, $geneBED);
	#my $bed_file = $_[0];
	#my $repeatBED = $_[1];
	#my $geneBED = $_[2];
	@cmd = split (/ /, "$BEDtoolPATH/closestBed -a $vfs_prefix.$runID.orphan.human.sorted.mapped.bed -b annotation/hg19.repeatmasker.anno -d >$repeatBED");  #add in repeat family and CLASS
	ExecCmd (\@cmd, $DEBUG_MODE);
	@cmd = split (/ /, "$BEDtoolPATH/closestBed -a $vfs_prefix.$runID.orphan.human.sorted.mapped.bed -b annotation/hg19.gene.info.in.f -d >$geneBED");#add NM
	ExecCmd (\@cmd, $DEBUG_MODE);
	my @gene = @{file_2_array($geneBED)};
	my @repeat = @{file_2_array($repeatBED)};
	my @NM = @{file_2_array($NM_BED)};
	#readID to NM
	my %readNM = ();
	for (@NM){
		my @s = split (/\t/, $_);
		$readNM{$s[3]}=$s[4];
	}
	#readID to repeat info
	my %readrepeat = ();
	for (@repeat){
		my @s = split (/\t/, $_);
		$readrepeat{$s[3]}="$s[8]\t$s[9]\t$s[10]\t$s[11]";
	}
	#infer exon/intron from RefSeq coordinate
	@cmd = split / /, "$BEDtoolPATH/closestBed -a $original_BED -b annotation/hg19.refseq.exon.bed -d >$original_BED.refseq.exon.distance";
	ExecCmd (\@cmd, $DEBUG_MODE);
	my @a = @{file_2_array("$original_BED.refseq.exon.distance")};
	#readID to exon distance
	my %readExonDistance = ();
	for (@a) {
		my @s = split (/\t/, $_);
		unless (defined $readExonDistance{$s[3]}){
			$readExonDistance{$s[3]}="$s[-1]";
		}
		else {
			if ($s[-1] < $readExonDistance{$s[3]}){
				$readExonDistance{$s[3]}="$s[-1]";
			}
		}
	}

	#merge human info
	open my $HinfowFH, ">", "$vfs_prefix.$runID.human.info" or die "Can't open $vfs_prefix.$runID.human.info : $!";
	for (@gene){
		my @s = split (/\t/, $_);
		my $exon_info;
		if ($readExonDistance{$s[3]} == 0){
			$exon_info = "T\t$readExonDistance{$s[3]}";
		}
		else {
			$exon_info = "F\t$readExonDistance{$s[3]}";
		}
		print $HinfowFH $_."\t".$readNM{$s[3]}."\t".$exon_info."\t".$readrepeat{$s[3]}."\n";
	}
	close ($HinfowFH);

	#parse HBV NM
	my %viral_readNM =();
	my @viral_NM_BED = @{file_2_array($viral_NM_BED)};
	for (@viral_NM_BED){
		my @s = split (/\t/, $_);
		$viral_readNM{$s[3]}=$s[4];
	}

	#parse HBV ORF info
	open my $VinfowFH, ">", "$vfs_prefix.$runID.viral.info" or die "Can't open $vfs_prefix.$runID.viral.info : $!";
	my @viral_BED = @{file_2_array($viral_BED)};;
	for (@viral_BED){
		my @s = split (/\t/, $_);
		my ($Vst, $Ved) = ($s[1], $s[2]);
		my $vNM = $viral_readNM{$s[3]};
		print $VinfowFH "$_\t$vNM";
		foreach my $hash_key (sort (keys %orfs)) {
			$orf_status{$hash_key} = in_range ($orfs{$hash_key}, $Vst, $Ved);
			printf $VinfowFH "\t%u", $orf_status{$hash_key};
		}
		print $VinfowFH "\n";
	}
	close ($VinfowFH);

	my @viral = @{file_2_array("$vfs_prefix.$runID.viral.info")};
	my @human = @{file_2_array("$vfs_prefix.$runID.human.info")};

	my %h = ();
	for (@viral){
		my @s = split (/\t/, $_);
		if ($s[3] =~ /(\V*)\/(\d*)/){
		$h{$1}= $_; #readID to viral_info
		}
	}
	my %human_h = ();
	my ($viral_read, $human_read);
	open my $result_outFH, '>', $RPm_result_file or die "Can't open $RPm_result_file for writing : $!";
	#1st pass to get the readname required for annotaiton
	my %readIDcore_required = ();
	for (@human){
		my $core;
		my @s = split (/\t/, $_);
		if ($s[3] =~ /(\V*)\/(\d*)/){
			$core = $1;
			my $h_read_strand = $2;
			my $q = '@'.$core;
			$readIDcore_required{$q} = "";
		}
	}
	my %F_read = %{hash_readID_seq_only_reads_req($forward_fastq, \%readIDcore_required)};
	my %R_read = %{hash_readID_seq_only_reads_req($reverse_fastq, \%readIDcore_required)};

	#merge viral and human information with respective read's sequence
	for (@human){
		my @s = split (/\t/, $_);
		my $core;
		if ($s[3] =~ /(\V*)\/(\d*)/){
			$core = $1;
			my $h_read_strand = $2;
			my $q = '@'.$core;
			if (defined $h{$core}){ #viral part has information for this read in another end
				if ($h_read_strand =~ /1/){
					($viral_read, $human_read) = ($R_read{$q}, $F_read{$q});
				}
				if ($h_read_strand =~ /2/){
					($viral_read, $human_read) = ($F_read{$q}, $R_read{$q});
				}
				print $result_outFH "$h{$core}\t$viral_read\t$_\t$human_read\n";
				$human_h{$q}="";
			}
		}
	}
	close ($result_outFH);

	#Fix the duplicate entries in output file
	my @raw_out = @{file_2_array($RPm_result_file)};
	my %v_h_readID_ck = ();
	my @result_order = ();
	my %seen = ();
	my @dedup_a = ();
	my %dedup = ();
	for my $raw_line (@raw_out){
		unless (defined $dedup{$raw_line}){
			push (@dedup_a, $raw_line);
			$dedup{$raw_line} = "";
		}
	}
	for my $dedup_line (@dedup_a){
		my @s = split (/\t/, $dedup_line);
		my $vReadID = $s[3];
		my $hDesc = $s[-10];
		my $hReadID = $s[-17];
		my $str = "$vReadID\t$hReadID";
		$str =~ s/://g;
		push (@result_order, $str); #record the order of readID as in the raw result file
		unless (defined $v_h_readID_ck{$str}){
			$v_h_readID_ck{$str} = $dedup_line;
		} else { #previous record already there
			my @s_already_result = split(/\t/, $v_h_readID_ck{$str});
			my $current_human_annotation = $s[-10];
			my $rep = join (", ", $s_already_result[-10], $current_human_annotation);
			$s_already_result[-10] = $rep;
			$v_h_readID_ck{$str} = join ("\t", @s_already_result)
		}
	}
	open my $fix_result_FH, '>', $RPm_result_file or die "Can't open $RPm_result_file for writing : $!";
	for my $q (@result_order){
		unless (defined $seen{$q}){
			print $fix_result_FH $v_h_readID_ck{$q}."\n";
			$seen{$q} = "";
		}
	}
	close ($fix_result_FH);

	#output remaining fastq for checking of viral mapped end, another side viral unmapped but still unmapped in human
	my %fastq = %{hash_fastq("$vfs_prefix.$runID.4h_map.fq")};
	open my $remainFH, ">", "$vfs_prefix.$runID.4h_map.fq.unmapped.fq" or die "Can't open $vfs_prefix.$runID.4h_map.fq.unmapped.fq for writing : $!";
	while ( my ($q, $data) = each(%fastq) ) {
		unless (defined $human_h{$q}) {
			print $remainFH $fastq{$q}."\n";
		}
	}
}



#sub in_range {
#	my @boundary = @{$_[0]};
#	my $start = $_[1];
#	my $end = $_[2];
#	my $status = 0;
#	@boundary = sort {$a <=> $b} @boundary; #min: $array[0]\nmax: $array[-1]
#	if (($start >= $boundary[0]) && ($end <= $boundary[-1])){
#		$status = 1;
#	}
#	#print "start: $start\tend: $end\tMax: $boundary[0]\tMin: $boundary[-1]\t$status\n";
#	return $status;
#}

sub hash_readID_seq_only_reads_req {
	my $file = $_[0];
	my %readID_is_needed = %{$_[1]};
	print STDERR "Hashing $file\n";
	my %h = ();

	my $is_gzip_file = DetermineFileType ($file);
	my $rFH = OpenFastqFileRead ($file, $is_gzip_file);

	while (1) {
		my $record = GetNextFastqRecord ($rFH);
		if (length ($record) == 0) {
			last;
		}
		my ($readID, $seq, $aux, $qual) = split /\n/, $record;
		chomp $seq;
		my @s_readID = split(/ /, $readID);
		my $readID_in = $s_readID[0];
		if (defined $readID_is_needed{$readID_in}){
			$h{$readID_in} = $seq;
		}
	}

	#CloseFastqFile ($Hfastq_rFH);
	CloseFastqFile ($rFH);

	return \%h;
}

sub hash_readID_seq {
	my $file = $_[0];
	print STDERR "Hashing $file to obtain sequences.\n";
	my %h = ();

	my $is_gzip_file = DetermineFileType ($file);
	my $rFH = OpenFastqFileRead ($file, $is_gzip_file);

	while (1) {
		my $record = GetNextFastqRecord ($rFH);
		if (length ($record) == 0) {
			last;
		}
		my ($readID, $seq, $aux, $qual) = split /\n/, $record;
		chomp $seq;
		$h{$readID} = $seq;
	}

	CloseFastqFile ($rFH);

	return \%h;
}


sub reads_factory {
	my $samtools = $_[0];
	my $bam_file = $_[1]; #viral sampe's BAM
	my $prefix = $_[2];
	my $forward_fastq = $_[3];
	my $reverse_fastq = $_[4];
	#my %fastqH = %{$_[5]};
	my $reference_file = $_[5];
	my $hds_read_file = "$prefix.4h_map.fq";
	my $hds_SAM_file = "$prefix.viral.mapped.sam";
	my $orphan_F_file = "$prefix.orphan_F.fq";
	my $orphan_R_file = "$prefix.orphan_R.fq";

	my $BAMrFH;
	if ($bam_file=~ /.bam/i){ #use samtools
		open $BAMrFH, " $samtools view $bam_file |" or die "Can't open pipe : $!";
	}
	if ($bam_file=~ /.sam$"/i){
		open $BAMrFH, "<", $bam_file or die "Can't open $bam_file for reading : $!";
	}
	open my $wFH, ">", $hds_read_file or die "Can't open $hds_read_file for writing : $!";
	open my $vBAM_FH, ">", $hds_SAM_file or die "Can't open $hds_SAM_file for writing : $!";
	open my $orfFH, ">", $orphan_F_file or die "Can't open $orphan_F_file for writing : $!";
	open my $orrFH, ">", $orphan_R_file or die "Can't open $orphan_R_file for writing : $!";
	my %required_fastq_entries = ();
	#1st pass to get what readID's fastq entries are needed
	while (<$BAMrFH>){
		chomp;
		my @s = split (/\t/, $_);
		my $flag = $s[1];
		#this read mapped, while another end unmapped
		unless ($flag & 0x4){ # i.e. viral-end is mapped
			if ($flag & 0x8){ # mate (another end of viral) unmapped in viral
				my $readID_q = '@'.$s[0];
				#print "need $readID_q\n";
				$required_fastq_entries{$readID_q}= "";
			}
		}
	}
	close ($BAMrFH);

	my %fastqH = %{hash_fastqs_only_reads_req($forward_fastq, $reverse_fastq, \%required_fastq_entries)};

	if ($bam_file=~ /.bam/i){ #use samtools
		open $BAMrFH, " $samtools view $bam_file |" or die "Can't open pipe : $!";
	}
	if ($bam_file=~ /.sam$"/i){
		open $BAMrFH, "<", $bam_file or die "Can't open $bam_file for reading : $!";
	}
	while (<$BAMrFH>){
		chomp;
		my @s = split (/\t/, $_);
		my $flag = $s[1];
		#this read mapped, while another end unmapped
		unless ($flag & 0x4){ # itself mapped
			if ($flag & 0x8){ # mate unmapped
				print $vBAM_FH "$_\n";
				my $readID_q = '@'.$s[0];
				my $readstrand_q = "";
				$readstrand_q = "R" if ($flag & 0x40); #itself is first fragment, so find reverse fragment
				$readstrand_q = "F" if ($flag & 0x80);
				if(defined $fastqH{$readID_q}{$readstrand_q}){
					my $record_entry = $fastqH{$readID_q}{$readstrand_q};
					my @s = split (/\n/, $record_entry);
					my $f = $s[0]."_".$readstrand_q;
					$s[0] = $f;
					my $fix = join ("\n", @s);
					print $wFH $fix."\n"; #record the other end to be mapped
					print $orfFH $fastqH{$readID_q}{"F"}."\n";
					print $orrFH $fastqH{$readID_q}{"R"}."\n";
				}
			}
		}
	}
	#SAM->BAM
	my @cmd = split (/ /, "$samtools view -Sb -T $reference_file $hds_SAM_file >$prefix.viral.mapped.bam");
	ExecCmd (\@cmd, $DEBUG_MODE);
	close ($wFH);
	close ($vBAM_FH);
	close ($orfFH);
	close ($orrFH);
	return \%fastqH;
}


sub hash_fastq { #NOT for hashing standard fastq, ONLY for RPmethod
	my $file = $_[0];
	print STDERR "Hashing $file\n";
	my %h = ();

	my $is_gzip_file = DetermineFileType ($file);
	my $Hfastq_rFH = OpenFastqFileRead ($file, $is_gzip_file);

	while (1) {
		my $record = GetNextFastqRecord ($Hfastq_rFH);
		if (length ($record) == 0) {
			last;
		}
		my ($readID, $seq, $aux, $qual) = split /\n/, $record;
		$h{$readID} = $record;
	}

	CloseFastqFile ($Hfastq_rFH);

	return \%h;
}

##  Read in a FASTQ file into a hash and keep only reads that are in the supplied hash
sub hash_fastqs_only_reads_req {
	my $forward_file = $_[0];
	my $reverse_file = $_[1];
	my %required_fastq_entries = %{$_[2]};
	print STDERR "Hashing $forward_file and $reverse_file\n";
	my %h = ();

	##  Process $forward_file
	my $is_gzip_file = DetermineFileType ($forward_file);
	my $fh = OpenFastqFileRead ($forward_file, $is_gzip_file);

	while (1) {
		my $record = GetNextFastqRecord ($fh);
		if (length ($record) == 0) {
			last;
		}
		my ($readID, $seq, $aux, $qual) = split (/\n/, $record);
		my @s_readID = split(/ /, $readID);
		my $readID_in = $s_readID[0];
		my $record_in = join ("\n", $readID_in, $seq, $aux, $qual);
		if (defined ($required_fastq_entries{$readID_in})){
			$h{$readID_in}{"F"} = $record_in;
		}
	}

	CloseFastqFile ($fh);

	##  Process $reverse_file
	$is_gzip_file = DetermineFileType ($reverse_file);
	$fh = OpenFastqFileRead ($reverse_file, $is_gzip_file);

	while (1) {
		my $record = GetNextFastqRecord ($fh);
		if (length ($record) == 0) {
			last;
		}
		my ($readID, $seq, $aux, $qual) = split (/\n/, $record);
		my @s_readID = split(/ /, $readID);
		my $readID_in = $s_readID[0];
		my $record_in = join ("\n", $readID_in, $seq, $aux, $qual);
		if (defined ($required_fastq_entries{$readID_in})){
			$h{$readID_in}{"R"} = $record_in;
		}
	}

	CloseFastqFile ($fh);
	print STDERR "done\n";
	return \%h;
}

##  Read in a FASTQ file into a hash and keep all reads
sub hash_fastqs {
	my $forward_file = $_[0];
	my $reverse_file = $_[1];
	my %required_fastq_entries = %{$_[2]};
	print STDERR "Hashing $forward_file and $reverse_file\n";
	my %h = ();

	##  Process $forward_file
	my $is_gzip_file = DetermineFileType ($forward_file);
	my $fh = OpenFastqFileRead ($forward_file, $is_gzip_file);

	while (1) {
		my $record = GetNextFastqRecord ($fh);
		if (length ($record) == 0) {
			last;
		}
		my ($readID, $seq, $aux, $qual) = split /\n/, $record;

		$h{$readID}{"F"} = $record;
	}

	CloseFastqFile ($fh);

	##  Process $reverse_file
	$is_gzip_file = DetermineFileType ($reverse_file);
	$fh = OpenFastqFileRead ($reverse_file);

	while (1) {
		my $record = GetNextFastqRecord ($fh);
		if (length ($record) == 0) {
			last;
		}
		my ($readID, $seq, $aux, $qual) = split /\n/, $record;

		$h{$readID}{"R"} = $record;
	}

	CloseFastqFile ($fh);
	print STDERR "done\n";
	return \%h;
}


sub file_2_array {
	my $curr_func = (caller(0))[3];
	my $count = @_;
	my $file_in = $_[0];
	my @out = ();
	open my $fh, '<', $file_in or die "Can't open $file_in for reading : $!";
	while (<$fh>){
		chomp;
		unless (/^\s*$/){ #push unless blank line
			unless (/^#/){
				push (@out, $_);
			}
		}
	}
	close ($fh);
	return \@out;
}


1;
