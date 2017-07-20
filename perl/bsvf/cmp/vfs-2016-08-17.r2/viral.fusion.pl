#!/usr/bin/perl
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


##  $LastChangedDate: 2013-12-13 17:56:17 +0800 (Fri, 13 Dec 2013) $
##  $LastChangedRevision: 1288 $
use diagnostics;
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/include";

##  Modules included with this program
use Essential; #VFS's clipseq module
use ReadPreprocess; #Read preprocess module
use RPmethod; #Read pair (F&R) module
use BWAauto; #BWA sampe/samse pipeline
use General; # Module for general subroutines

##  Modules not provided with this program
use Bio::DB::Sam;
use AppConfig;
use AppConfig::Getopt;
use Pod::Usage;

DisplayTime ();

########################
##  Global variables  ##
########################

##  Indicate whether or not threads are available
$| = 1;
my $threads_usable = eval 'use threads; 1';
if ($threads_usable) {
        use threads;
        print "Perl threading enabled\n";
} else {
        print "No threading is possible. Please install perl module: threads\n";
}

## 
my $global_sleepTime = 0;
my $DEBUG_MODE = 1;
##############################################
##   Set up AppConfig; parameters used in   ##
##     command-line and configuration file  ##
##############################################

##  Create AppConfig and AppConfig::Getopt objects
my $config = AppConfig -> new ({
	CASE      => 1,
	GLOBAL => {
	DEFAULT => undef,     ##  Default value for new variables
	}
});

##  Binary (true/false) arguments
$config -> define ("ReadPreprocess!", {
	DEFAULT => 0,
});
$config -> define ("SCmethod!", {
	DEFAULT => 0,
});
$config -> define ("ViralSCmapping!", {
	DEFAULT => 0,
});
$config -> define ("analyzeSCfiles!", {
	DEFAULT => 0,
});
$config->define("SCprepparse!", {
	DEFAULT => 0,
});
$config->define("SCparse!", {
	DEFAULT => 0,
});
$config->define("RPmethod!", {
	DEFAULT => 0,
});
$config->define("readlevelAnalysis!", {
	DEFAULT => 0,
});
$config->define("AssembleSC!", {
	DEFAULT => 0,
});
$config->define("doTargetedAssembly!", {
	DEFAULT => 0,
});
$config -> define ("cleanup!", {
	DEFAULT => 0,
});
$config -> define ("verbose!", {
	DEFAULT => 0,
});  #default no verbose
$config -> define ("help!", {
	DEFAULT => 0,
});                            ##  Help screen

##  Positive integer
$config->define("thread", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => 8,
	ARGS => "=i"
});

##  Positive integer if insert size known; non-integer if unknown and it
##  is up to BWA to determine it. #The detection step is in Essential::tx_targeted_assembly
$config->define("insertSIZE", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "NA",
	ARGS => "=s"
});

##  Software paths (also takes into account the user's ${PATH} variable).
$config -> define ("bwa", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "bwa",
	ARGS => "=s",
});  #default to use system-wide bwa (need BWA 0.6 series)
$config -> define ("samtools", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "samtools",
	ARGS => "=s",
});
$config -> define ("blast", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "blastall",
	ARGS => "=s",
});
$config -> define ("cap3", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "cap3",
	ARGS => "=s",
});
$config -> define ("ssake", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "ssake",
	ARGS => "=s",
});

$config -> define ("minLEN", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => 10,
	VALIDATE => \&CheckPositiveInt,
	ARGS => "=i",
});  #default to require clipped sequence to be >= 10bp
$config -> define ("phredQ", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "NA",
	VALIDATE => \&CheckPhredQ,
	ARGS => "=s",
});
$config -> define ("desiredQ", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => 10,
	VALIDATE => \&CheckPositiveInt,
	ARGS => "=i",
});
$config -> define ("emitThreshold", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => 35,
	VALIDATE => \&CheckPositiveInt,
	ARGS => "=i",
});

##  Paths to files
$config -> define ("viralFa", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "hbv.decoyseq/hbv4.fa",
	ARGS => "=s",
});
$config -> define ("ntDB", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "include/blast.db/nt",
	ARGS => "=s",
});
$config -> define ("humanFa", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "/expt/raid0-A/marco/hbv.fusion/vfs/hg19.clean/bwa060/hg19.clean.gatk.bwa061",
	ARGS => "=s",
});
$config -> define ("humanDecoy", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "/expt/raid0-A/marco/hbv.fusion/vfs/1kg.decoyseq/hs37d5.fa",
	ARGS => "=s",
});
$config -> define ("bedtoolPATH", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "/home/marco/tools/bedtools/BEDTools-Version-2.16.2/bin",
	ARGS => "=s",
});

##  List of values
$config -> define ("clippedSeqKeywords", {
	ARGCOUNT => AppConfig::ARGCOUNT_LIST,
	ARGS => "=s@",
});
$config -> define ("mappedSeqKeywords", {
	ARGCOUNT => AppConfig::ARGCOUNT_LIST,
	ARGS => "=s@",
});
$config -> define ("orf", {
	ARGCOUNT => AppConfig::ARGCOUNT_LIST,
	ARGS => "=s@",
});

##  Path to configuration file
$config -> define ("config", {
	ARGCOUNT => AppConfig::ARGCOUNT_ONE,
	DEFAULT => "vfs.conf", #VFS's default config file is vfs.conf
	ARGS => "=s",
});

############################
##   Parse config file    ##
############################
my @original_ARGV = @ARGV;
#Find if custom config file is given
$config -> getopt ();
if (!-e $config -> get ("config")) {
	printf STDERR "EE\tConfiguration file \'%s\' does not exist!\n", $config -> get ("config");
	exit (1);
} else {
	my $conf_file_used = $config -> get ("config");
	print STDERR "Config file: $conf_file_used\n";
}
$config->file ($config -> get ("config"));

################################################################
##              Process the command-line options              ##
##   Command-line specificed parameters gain higher priority  ##
################################################################
$config -> getopt (\@original_ARGV);
#$do_SCmethod = $config->get('SCmethod');

##   Get update from command-line, if any
my $do_preprocess = $config->get('ReadPreprocess'); # 1: to do reads preprocessing
my $do_SCmethod = $config->get('SCmethod'); #Master switch for CS method
my $do_viral_SC_mapping = $config->get('ViralSCmapping');
my $do_analyze_SC_files = $config->get('analyzeSCfiles'); #contains $do_SC_prep_parse and $do_SC_parse
my $do_SC_prep_parse = $config->get('SCprepparse');
my $do_SC_parse = $config->get('SCparse');
my $do_RPmethod = $config->get('RPmethod');
my $do_read_level_analysis = $config->get('readlevelAnalysis');
my $do_assemble_CS = $config->get('AssembleSC');
my $do_targeted_assembly = $config->get('doTargetedAssembly');
my $do_cleanup = $config->get('cleanup');
my $verbose = $config->get('verbose');
my $thread = $config->get('thread');
my $isize_supplied = 0;  ##  Default is no insert size supplied
if ($config -> get ('insertSIZE') =~ /^\d+$/) {
	$isize_supplied = $config -> get ('insertSIZE');
}
my $bwa = $config->get('bwa');
my $samtools = $config->get('samtools');
my $blast = $config->get('blast');
my $cap3 = $config->get('cap3');
my $SSAKE_bin = $config->get('ssake');
my $minLEN = $config->get('minLEN');
my $phredQ = $config->get('phredQ');
my $desiredQ = $config->get('desiredQ');
my $emitThreshold = $config->get('emitThreshold');
my $viral_fa = $config->get('viralFa');
my $human_fa = $config->get('humanFa');
my $humanDecoy = $config->get('humanDecoy');
my $ntDB = $config->get('ntDB');
my $BEDtoolPATH= $config->get('bedtoolPATH');
my $clippedSeqKeywordsRef = $config->get('clippedSeqKeywords');
my $mappedSeqKeywordsRef = $config->get('mappedSeqKeywords');
my $orfRef = $config->get('orf');

##   Show help screen if --help given
if ($config -> get ("help")) {
	pod2usage (-exitval => 1, -verbose => 1);
}

###################################
##   Validate AppConfig options  ##
###################################

#Have to check the exact path for alias (e.g. samtools but not home/tools/samtools)

##  Retrieve the full path of an executable if no path (neither absolute nor relative)
##    was given
$bwa = ExpandProgramPath ($bwa);
$samtools = ExpandProgramPath ($samtools);
$blast = ExpandProgramPath ($blast);
$cap3 = ExpandProgramPath ($cap3);
$SSAKE_bin = ExpandProgramPath ($SSAKE_bin);

##  Check paths to executable programs
if (!CheckProgram ($bwa, __FILE__, __LINE__, $DEBUG_MODE)) {
  exit (1);
}
if (!CheckProgram ($samtools, __FILE__, __LINE__, $DEBUG_MODE)) {
  exit (1);
}
if (!CheckProgram ($blast, __FILE__, __LINE__, $DEBUG_MODE)) {
  exit (1);
}
if (!CheckProgram ($cap3, __FILE__, __LINE__, $DEBUG_MODE)) {
  exit (1);
}
if (!CheckProgram ($SSAKE_bin, __FILE__, __LINE__, $DEBUG_MODE)) {
  exit (1);
}

##  Check paths to files
if (!CheckFile ($viral_fa, __FILE__, __LINE__, $DEBUG_MODE)) {
  exit (1);
}
if (!CheckFile ($human_fa, __FILE__, __LINE__, $DEBUG_MODE)) {
  exit (1);
}
if (!CheckFile ($humanDecoy, __FILE__, __LINE__, $DEBUG_MODE)) {
  exit (1);
}
if (!CheckPath ($BEDtoolPATH, __FILE__, __LINE__, $DEBUG_MODE)) {
  exit (1);
}

##  Check paths to indexes
my @IndexSuffixes = qw (amb ann bwt pac sa);
foreach my $tmp_index (@IndexSuffixes) {
	my $index_check = 0;
	$index_check = CheckFile ($viral_fa.".".$tmp_index, __FILE__, __LINE__, $DEBUG_MODE);
	if (!$index_check) {
		##  Run bwa
		bwa_index ($bwa, $viral_fa, 1);
	}

	$index_check = CheckFile ($human_fa.".".$tmp_index, __FILE__, __LINE__, $DEBUG_MODE);
	if (!$index_check) {
		bwa_index ($bwa, $human_fa, 0);
	}

	$index_check = CheckFile ($humanDecoy.".".$tmp_index, __FILE__, __LINE__, $DEBUG_MODE);
	if (!$index_check) {
		bwa_index ($bwa, $humanDecoy, 0);
	}
}

##  Check path to NT database
@IndexSuffixes = qw (nhd nhi nhr nin nnd nni nog nsd nsi nsq);
foreach my $tmp_index (@IndexSuffixes) {
	my $index_check = 0;
	$index_check = CheckFile ($ntDB.".00.".$tmp_index, __FILE__, __LINE__, $DEBUG_MODE);
	if (!$index_check ) {
		printf STDERR "Error:  NT database not found in \'%s\'.  Database file with the suffix \'.00.%s\' could not be found.  Please re-download the NT database from NCBI.\n", $ntDB, $tmp_index;
		exit (1);
	}
}


#########################################
##   Process remaining options         ##
##     (Not preceded by "--" or "-".)  ##
#########################################

if (@ARGV < 2) {
	pod2usage (-exitval => 1, -verbose => 0);
}

my $runID = $ARGV[0];
my $F_readfile = $ARGV[1];
my $R_readfile = $ARGV[2];

CheckFile ($F_readfile, __FILE__, __LINE__, $DEBUG_MODE);
if (defined $R_readfile) {
	CheckFile ($R_readfile, __FILE__, __LINE__, $DEBUG_MODE);
}

my @reads = ();
push @reads, $F_readfile;
push @reads, $R_readfile if (defined $R_readfile);
#my @fastq_phred_a = @{rtn_fastq_phred($F_readfile, $R_readfile)};
my @fastq_phred_a = @{rtn_fastq_phred(\@reads)};

##########################################
#Intermediate files suffix: Do NOT change#
##########################################
my $vfs_prefix = "vfs_dev"; #"vfs_dev";
my $CS_MS_ori_read_file = "$vfs_prefix.$runID.CSmethod.fa";
my $CS_file = "$vfs_prefix.$runID.CS";
my $MS_file = "$vfs_prefix.$runID.MS";

##############################################################
#Parameters that should ONLY be modified by experienced users#
##############################################################

my $sensitive_tx_denovo_assembly = 1; # This mode is "ON" by default

#####################################
##  Display essential information  ##
#####################################
print STDERR "######################################\n";
print STDERR "##Sequence reads related infromation #\n";
print STDERR "######################################\n";
print "Process ID of $0: $$\n";
print "Run ID: $runID\n";
print STDERR "Forward read: $F_readfile\n";
unless (!defined $R_readfile){
	print STDERR "Reverse read: $R_readfile\n";
}
else {
	print STDERR "No Reverse read is used. Read-Pair method will be skipped\n";
	print STDERR "No Reverse read is used. Targeted assembly method will be skipped\n";
	$do_RPmethod = 0;
	$do_targeted_assembly = 0;
}
print STDERR "Verbosity: $verbose\n";
print STDERR "Number of thread to use for third-party tools: $thread\n";
print STDERR "External Fastq insert size: ";
if ($isize_supplied != 0) {
	print STDERR $isize_supplied, "\n";
}
else {
	print STDERR "N/A (use bwa)\n";
}
print STDERR "Read preprocessing - Fastq Phred quality: $phredQ\n";

## Detect the Phred coding of fastq files, correct the user's supplied parameter, if discrepency is found.
my %phred_c = ();
for my $fq_phred (@fastq_phred_a){
	my @s = split (/\t/, $fq_phred);
	print "Phred quality of $s[0] detected by VFS: $s[1]\n";
	$phred_c{$s[1]} = "";
}
my $phred_cnt = keys %phred_c;
if (($phred_cnt > 1)){ #The 2 fastq files have different fastq quality
	print STDERR "2 fastq files has different Phred quality.Please check the input files\n";
	exit (1);
}
else { #The 2 fastq files have same Phred Quality, but the Phred Quality specified by the user is wrong
	for my $k (keys %phred_c){ # 1 key only
		if ($phredQ =~ /NA/g){
			print STDERR "User is not sure about the Phred Quality of the fastq files: $phredQ\nGood news is VFS detected a consensus Phred Quality: $k\n";
			print STDERR "VFS will use the detected Phred quality $k during reads pre-processing\n";
			$phredQ = $k;
		}
		else {
			print STDERR "PhredQ specificed by user: $phredQ\nVFS detected consensus Phred Quality: $k. VFS will use the detected Phred quality $k instead of $phredQ during reads pre-processing\n" if ($k != $phredQ );
			print STDERR "PhredQ specificed by user: $phredQ\nVFS detected consensus Phred Quality: $k. No discrepency detected.\n" if ($k == $phredQ );
		}
	}
}

print STDERR "Read preprocessing - Re-implemented BWA's -q base-call quality: $desiredQ\n";
print STDERR "Read preprocessing - Emit threshold: $emitThreshold\n";
print STDERR "Clipped sequence module - minLEN: $minLEN\n";
print STDERR "###################################################\n";
print STDERR "##Modules to be executed. (1)Execute; (0) Skipped #\n";
print STDERR "###################################################\n";
print STDERR "Read Preprocessing: $do_preprocess\n";
print STDERR "Clipped-Seq module: $do_SCmethod\n";
if ($do_SCmethod == 1){
	print STDERR "SC - mapping: $do_viral_SC_mapping\n";
	print STDERR "SC - analyzeSCfiles: $do_analyze_SC_files\n";
	print STDERR "SC - SCprepparse: $do_SC_prep_parse\n";
	print STDERR "SC - SCparse: $do_SC_parse\n";
	print STDERR "SC - readlevelAnalysis: $do_read_level_analysis\n";
	print STDERR "SC";
} else {
	print STDERR "All methods under Clipped-Seq module are to be skipped because Clipped-Seq module is skipped: $do_SCmethod\n";
}
print STDERR "Read-pair module: $do_RPmethod\n";
print STDERR "Targeted assembly: $do_targeted_assembly\n";
print STDERR "Temporary files cleanup: $do_cleanup\n";
print STDERR "##################################################\n";
print STDERR "##Paths to reference files and third-party tools #\n";
print STDERR "##################################################\n";
print STDERR "BWA bin: $bwa\n";
print STDERR "Samtools bin: $samtools\n";
print STDERR "BLAST bin: $blast\n";
print STDERR "BedTools bin: $BEDtoolPATH\n";
print STDERR "CAP3 bin: $cap3\n";
print STDERR "SSAKE bin: $SSAKE_bin\n";
print STDERR "Viral fasta file: $viral_fa\n";
print STDERR "Human fasta file: $human_fa\n";
print STDERR "Human decoy file: $humanDecoy\n";
print STDERR "nt database: $ntDB\n";
print STDERR "#############################################\n";
print STDERR "## Keywords used for Clipped-Seq searching ##\n";
print STDERR "#############################################\n";
print STDERR "Keywords for clipped sequences:  ";
if (scalar (@$clippedSeqKeywordsRef) == 0) {
	printf STDERR "None\n";
	printf STDERR "Warning: Clipped-Seq module may not execute properly. VFS will terminate now\n";
 	exit (1);
}
else {
	printf STDERR lc (join (", ", @$clippedSeqKeywordsRef)); #Display only
}
printf STDERR "\n";

print STDERR "Keywords for mapped segments:  ";
if (scalar (@$mappedSeqKeywordsRef) == 0) {
	printf STDERR "None\n";
	printf STDERR "Warning: Mapped-Seg module may not execute properly. VFS will terminate now\n";
	exit (1);
}
else {
	printf STDERR lc (join (", ", @$mappedSeqKeywordsRef));#Display only
}
printf STDERR "\n";

print STDERR "Viral ORFs to check:  ";
if (scalar (@$orfRef) == 0) {
	printf STDERR "None\n";
}
else {
	printf STDERR "\n";
	for my $orf_tmp (@$orfRef) {
		##  Returns 1 if there was an error in the ORF format
		my $result = CheckORF ($orf_tmp);

		##  Error found; print no error message since CheckORF would have
		if ($result) {
			exit (1);
		}

		my ($orf_name, $orf_start, $orf_end) = split (/,/, $orf_tmp);

		printf STDERR "\t<%s [%u...%u]>\n", $orf_name, $orf_start, $orf_end;
	}
}

#Let user review the parameters used
sleep $global_sleepTime;

############################
#  Program set-up
############################

##  Create log directory
CreateLogDirectory ("$vfs_prefix.$runID.log", $DEBUG_MODE);

############################
#Pre-processing sequence reads, if necessary
############################
my ($F_readfile_f, $R_readfile_f); #Always use the pre-processed/renamed file in any down-steram analysis
if ($do_preprocess == 1){
	my @out_files =@{preprocess_fastq($F_readfile, $R_readfile, $phredQ, $desiredQ, $emitThreshold)};
	$F_readfile_f = $out_files[0];
	$R_readfile_f = $out_files[1] if (defined $out_files[1]);
}
else {
	$F_readfile_f = $F_readfile;
	$R_readfile_f = $R_readfile if (defined $R_readfile);
	print STDERR "Reads preprocessing skipped (As instructed)\n";
}

if ($do_SCmethod == 0){
	print STDERR "Clipped-Seq method skipped (As instructed)\n";
}
else {
	############################
	#Mapping step of Soft-clipped method
	############################
	# BWA-SW
	my $aligner_runID_bam;
	if ($do_viral_SC_mapping == 1){
		if (defined $R_readfile) {
			$aligner_runID_bam = bwa_sw("$vfs_prefix.$runID", $bwa, $F_readfile_f, $R_readfile_f, $viral_fa, $thread, "$vfs_prefix.$runID.log");
		} else {
			$aligner_runID_bam = bwa_sw("$vfs_prefix.$runID", $bwa, $F_readfile_f, "", $viral_fa, $thread, "$vfs_prefix.$runID.log");
		}
	}
	############################
	#Analyze alignment files for viral fusion
	############################
	if ($do_analyze_SC_files == 1){
		if ($do_SC_prep_parse == 1){
			#rtn_CS($aligner_runID_bam, $viral_fa, $CS_file, $MS_file); #extract clipped sequences from alignemnt file
			rtn_CS("$vfs_prefix.$runID.bam", $viral_fa, $CS_file, $MS_file, $CS_MS_ori_read_file, $F_readfile_f, $R_readfile_f);
			extract_fa_by_len($CS_file, "$CS_file.min$minLEN", $minLEN, $MS_file, "$MS_file.min$minLEN"); #enhance specificity

			##CS core: viral part
			map_2_genome($bwa, "$CS_file.min$minLEN", $human_fa, $thread, "$vfs_prefix.$runID.CSrd.sam", "$vfs_prefix.$runID.log");
			blast_run($blast, "$CS_file.min$minLEN", $ntDB, "$vfs_prefix.$runID.nt.CSrd.out", "5", $thread, "$vfs_prefix.$runID.log");

			##MS core: quality checking of Viral part
			#map_2_genome($bwa, $MS_file, $viral_fa, $thread, "$vfs_prefix.$runID.MS.sam", "$vfs_prefix.$runID.log");
			#blast_run($blast, $MS_file, $ntDB, "$vfs_prefix.$runID.nt.MS.out", "5", $thread, "$vfs_prefix.$runID.log");
			map_2_genome($bwa, "$MS_file.min$minLEN", $viral_fa, $thread, "$vfs_prefix.$runID.MS.sam", "$vfs_prefix.$runID.log");
			blast_run($blast, "$MS_file.min$minLEN", $ntDB, "$vfs_prefix.$runID.nt.MS.out", "5", $thread, "$vfs_prefix.$runID.log");
		}
		##Parse for output
		if ($do_SC_parse == 1){
			print STDERR "Parsing results ";
			if ($threads_usable){
				print STDERR "with threads\n";
				my $thr1 = threads->create(\&blast_parse, "$vfs_prefix.$runID.nt.CSrd.out");
				my $thr2 = threads->create(\&blast_parse, "$vfs_prefix.$runID.nt.MS.out");
				my @ReturnData = $thr1->join();
				my $rtn = join(', ', @ReturnData);
				print STDERR "Parsed $vfs_prefix.$runID.nt.CSrd.out\n" if ($rtn eq "1");
				my @ReturnData2 = $thr2->join();
				my $rtn2 = join(', ', @ReturnData2);
				print STDERR "Parsed $vfs_prefix.$runID.nt.MS.out\n" if ($rtn2 eq "1");
			} else {
				print STDERR "with single thread\n";
				blast_parse("$vfs_prefix.$runID.nt.CSrd.out");
				print STDERR "Parsed $vfs_prefix.$runID.nt.CSrd.out\n";
				blast_parse("$vfs_prefix.$runID.nt.MS.out");
				print STDERR "Parsed $vfs_prefix.$runID.nt.MS.out\n";
			}
		}
		#Clipped Seq read_level_analysis
		if ($do_read_level_analysis == 1){
			CS_read_level_analysis("$vfs_prefix.$runID", "$vfs_prefix.$runID.CSm.out", "$vfs_prefix.$runID.nt.MS.out.parsed", "$vfs_prefix.$runID.nt.CSrd.out.parsed", $MS_file, "$CS_file.min$minLEN", $CS_MS_ori_read_file, $clippedSeqKeywordsRef, $mappedSeqKeywordsRef);
		} else {
			print STDERR "Read-level analysis of CS module skipped (As instructed)\n";
		}
		if ($do_assemble_CS == 1){
			my @a = @{Assemble_CS("$vfs_prefix.$runID", $cap3, "$vfs_prefix.$runID.CSm.out", "$vfs_prefix.$runID.log")}; ###Need to work on this part
			###The following results will be used to merge with the output file by read-level anlaysis (the read-level result with assembly result is more specific) - START
			my %contigs_to_anonymous_array_of_simIDlist_simID_hash_to_readID = %{$a[0]};
			my %MS_CS_simID_to_readID = %{$a[1]};
			for my $CS_or_MS ( keys %contigs_to_anonymous_array_of_simIDlist_simID_hash_to_readID) {
				my %contig_hash = %{$contigs_to_anonymous_array_of_simIDlist_simID_hash_to_readID{$CS_or_MS}};
				for my $contig ( keys %contig_hash) {
					print "$CS_or_MS\t$contig\t";
					for my $sim_ID (@{$contig_hash{$contig}}){
						print "$sim_ID: $MS_CS_simID_to_readID{$CS_or_MS}{$sim_ID} ";
					}
					print "\n";
				}
			}
			###The following results will be used to merge with the output file by read-level anlaysis (the read-level result with assembly result is more specific) - END
		} else {
			print STDERR "No assembly will be done on clipped sequences in CS module (As instructed)\n";
		}
	} #Analysis ends
}

############################
#Read-pair method
############################
if ($do_RPmethod == 1){
	exec_RPmethod($vfs_prefix, $runID, "$vfs_prefix.$runID.RPm.out", $thread, $F_readfile_f, $R_readfile_f, $humanDecoy, $viral_fa, $bwa, $samtools, $BEDtoolPATH, $orfRef);
} else {
	print STDERR "Read-pair method skipped (As instructed)\n";
}

############################
#Extension phase - targeted de novo assembly(viral seed)
############################
if ($do_targeted_assembly == 1){
	tx_targeted_assembly("$vfs_prefix.$runID", "$vfs_prefix.$runID.CSm.out", $CS_MS_ori_read_file, "$vfs_prefix.$runID.RPm.out", $F_readfile_f, $R_readfile_f, $human_fa, $humanDecoy, $thread, $bwa, $SSAKE_bin, $samtools, $isize_supplied, $sensitive_tx_denovo_assembly, "$vfs_prefix.$runID.log"); #(1)Clipped-Seq result fle, (2) Read-pair result file, (3) Full-set Forward read, (4) Full-set Reverse Read
} else {
	print STDERR "VFS will not attempt to re-constructed fusion transcript (As instructed)\n";
}

############################
#Clean-up files
############################
if ($do_cleanup == 1){
	clean_tempfiles("$vfs_prefix.$runID");
	#CleanLogDirectory ($runID, $DEBUG_MODE);
} else {
	print STDERR "No temporary file clean up will be done (As instructed)\n";
}
DisplayTime ();


=pod

=head1 NAME

viral.fusion.pl -- Discover viral integration events and fusion transcripts using high-throughput sequencing data.

=head1 SYNOPSIS


 __      ___           _ ______         _              _____              ___  
 \ \    / (_)         | |  ____|       (_)            / ____|            |__ \ 
  \ \  / / _ _ __ __ _| | |__ _   _ ___ _  ___  _ __ | (___   ___  __ _     ) |
   \ \/ / | | '__/ _` | |  __| | | / __| |/ _ \| '_ \ \___ \ / _ \/ _` |   / / 
    \  /  | | | | (_| | | |  | |_| \__ \ | (_) | | | |____) |  __/ (_| |  / /_ 
     \/   |_|_|  \__,_|_|_|   \__,_|___/_|\___/|_| |_|_____/ \___|\__, | |____|
                                                                     | |       
                                                                     |_|       

B<viral.fusion.pl> [OPTIONS] <runID> <Forward read> [Reverse read]

Readname ID's format has to adhere to illumina v.1.8 standard. The script "prep.fastq.2.Illumina.1.8.gz.pl" under misc/ will help you on this

Reads files must be in SAME directory of viral.fusion.pl and could be gzipped with file suffix .gz 

For detailed help. "perldoc viral.fusion.pl"

=head1 DESCRIPTION

This script implements a pipeline for using high-throughput sequencing data to detect viral fusions.  Various programs are used by this pipeline, which should exist in the user's path or have their full paths given to this script.

The script has two mandatory options that must appear I<last>.  The first of these is the run_ID.  Afterwards is the path to the forward reads.  A third option specifies the path to the reverse reads.  This option can be omitted if single-ended reads are being processed.

Various options can appear before these mandatory options.  They are all preceded by "--".  These options can also appear in a configuration file.  They are B<case sensitive> and are as follows:

=head1 OPTIONS

=over 5

=item --config F<file>

Path to the configuration file.  Default is F<vfs.conf> in the current directory.

=item --ReadPreprocess

Indicate that sequence reads need to be quality-trimmed before subjected to fusion discovery

=item --SCmethod

Indicate that the whole SC method should be executed. Sub-processes should be specificed. Those includes --ViralSCmapping; --analyzeSCfiles; --SCprepparse; --SCparse and --readlevelAnalysis

=item --ViralSCmapping

Indicate that viral SC mapping should be performed. --SCmethod has to be enabled

=item --analyzeSCfiles

Indicate that SC files should be analyzed. --SCmethod has to be enabled

=item --SCprepparse

Indicate that the SC files should be prepared for parse. --SCmethod has to be enabled

=item --SCparse

Indicate that the SC files should be parsed. --SCmethod has to be enabled

=item --readlevelAnalysis

Perform read level analysis. --SCmethod has to be enabled

=item --AssembleSC

Assemble the clipped sequences found by --SCmethod.

=item --RPmethod

Indicate that the Read-Pair method should be run.

=item --doTargetedAssembly

Perform targeted assembly.

=item --cleanup

Clean up afterwards.

=item --verbose

Verbose output.

=item --thread I<integer>

Indicate number of threads.  Must be a value larger than 0.

=item --insertSIZE I<size>

Provide the insert size, if known, as an integer.  If unknown, then provide a non-integral value and bwa will be used to determine it.

=item --bwa F<path>

Full system path to the bwa binary (needs to be the version 0.6 series).

=item --samtools F<path>

Full system path to the samtools binary

=item --blast F<path>

Path to the blast binary

=item --cap3 F<path>

Full system path to the CAP3 binary

=item --ssake F<path>

Full system path to the SSAKE binary

=item --minLEN I<integer>

Minimum sequence length of clipped sequences. Should be >= 10 bp

=item --phredQ I<integer>

Parameter for read-preprocessing. Phred encoding scheme for fastq files.  Should be either 33/64. Use "NA" if you are not sure

=item --desiredQ I<integer>

desired is a parameter for read-preprocessing. This parameter is the same as bwa trimming algorithm –q: argmax_x{\sum_{i=x+1}^l(INT-q_i)}

=item --emitThreshold I<integer>

emitThreshold is a parameter for read-preprocessing. The minimal length (bp) of either end of trimmed sequence reads required to return both ends

=item --viralFA F<file>

Full system path to the viral genome reference file

=item --ntDB F<file>

Full system path to the nt database. Make sure the database has been built / extracted successfully. You can download the nt database at ftp://ftp.ncbi.nlm.nih.gov/blast/db/

=item --humanFA F<file>

Full system path to the human genome reference file. It is a single file comprising all chromosomes

=item --humanDecoy F<file>

Full system path to the human decoy reference file. It is a single file comprising all chromosomes.

=item --bedtoolPATH F<path>

Full system path to the BED tools /bin directory

=item --clippedSeqKeywords I<string 1> 

One keyword for clipped sequences. If more than 1 keyword is to be specified in command line, repeat the argument and the keyword multiple times
e.g. --clippedSeqKeywords Keyword1 F<-->clippedSeqKeywords Keyword 2

=item --mappedSeqKeywords I<string 1> 

One keyword for mapped segments (more than one is possible).
e.g. --mappedSeqKeywords Keyword1 F<-->mappedSeqKeywords Keyword 2


=back

=head1 EXAMPLE

=over 5

./viral.fusion.pl run1 forward.fq reverse.fq

=back

=head1 AUTHOR

Jing-Woei Li (Marco) <marcoli@cuhk.edu.hk>

=head1 COPYRIGHT

Copyright (C) 2012, Marco Li, All rights reserved.


