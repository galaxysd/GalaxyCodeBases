#!/usr/bin/perl -w
###################################################################################
#
#  VirusFinder, a fully automatic pipeline for efficient and accurate detection of viruses,
#               viral mutations, and viral integration sites in host genomes through next
#               generation sequencing data.
#
#  VirusFinder is free software
#
#  Version 2
#  Last update: 02/18/2014
#
#  Contact:       Zhongming Zhao: zhongming.zhao@vanderbilt.edu
#                 Qingguo Wang:   qingguo.wang@vanderbilt.edu
#
#  Organization:  Bioinformatics and Systems Medicine Laboratory
#                 Vanderbilt University Medical Center
#                 Nashville, Tennessee, USA
#
#  detect_mutation.pl is part of VirusFinder 2. It detects snps and indels in virus genome.
#
###################################################################################
# Changelog:
# 02/18/2014   Software release
###################################################################################

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use IO::Handle();
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use Mosaic;

my @usage;
push @usage, "\nUsage:  detect_mutation.pl  <-c configuration file>  [options] \n\n";

push @usage, "Options:\n";
push @usage, "  -h, --help     Displays this information.\n";
push @usage, "  -c, --config   Configuration file <required>.\n";
push @usage, "  -o, --output   Full path of a directory to store results, default is current working directory.\n";
push @usage, "  -v, --virus    Sequence file (in fasta format) of the virus. User can also provide virus's ID (or prefix\n";
push @usage, "                 of ID) in the virus database, e.g. gi_310698439_ref_NC_001526.2__Human_papillomavirus_type_16.\n";
push @usage, "                 This argument is mandatory unless user provides a file results-virus-top1.fa under\n";
push @usage, "                 the directory specified by the user.\n";
push @usage, "  --fq1          Sequencing reads are FASTQ file. This file is mandatory unless there is a file\n";
push @usage, "                 unmapped.1.fq in the directory specified by the user. The file unmapped.1.fq\n";
push @usage, "                 can be created from FASTQ files using preprocess.pl, a script for subtracting\n";
push @usage, "                 reads mapped to the host genome.\n";
push @usage, "  --fq2          Sequencing reads are FASTQ file. If this file is not provided and there isn't a\n";
push @usage, "                 file unmapped.2.fq in the directory specified by the user. The program will treat\n";
push @usage, "                 the input data as single end. Unmapped.2.fq can be created from FASTQ files\n";
push @usage, "                 using preprocess.pl, a script for subtracting reads mapped to the host genome.\n";
push @usage, "  -m, --markdup  Mark duplicate reads. Accepted inputs: y[es], n[o]. Default value is y[es].\n";
push @usage, "                 Duplicate reads will not be used for variant calling. For ultra-deep amplicon\n";
push @usage, "                 sequencing data, 'no' should be used for this argument.\n\n";
push @usage, "Example:\n";
push @usage, "  perl detect_mutation.pl -c config.txt -o /scratch/kingw/VirusFinder/simulation\n";
push @usage, "  perl detect_mutation.pl -c config.txt -v gi_310698439_ref_NC_001526.2\n\n";
push @usage, "  Before running the first example, user has to place a virus sequence file results-virus-top1.fa\n";
push @usage, "  under the folder /scratch/kingw/VirusFinder/simulation. Otherwise, the program will complain.\n";
push @usage, "  The file results-virus-top1.fa can be created using script detect_virus.pl.\n";
push @usage, "  In the second case, the program will look up virus database to match string 'gi_310698439_ref_NC_001526.2'.\n";
push @usage, "  The sequence of the matched virus will be used as reference for read mapping and variant calling.\n\n";
push @usage, "Output:\n";
push @usage, "  The program will create a file mutations.vcf that contains identified mutations, i.e. SNPs\n";
push @usage, "  and indels, in the virus genome.\n\n";


my $help;
my $config_file;
my $output_dir;
my $read_file1;
my $read_file2;
my $virus_sequence;
my $consensus_virus_seq;
my $paired = 1;
my $markdup = 'y';

GetOptions
(
 'h|help|?'    => \$help,
 'config=s'    => \$config_file,
 'output=s'    => \$output_dir,
 'virus=s'     => \$virus_sequence,
 'fq1=s'       => \$read_file1,
 'fq2=s'       => \$read_file2,
 'markdup=s'   => \$markdup,
);

if ($help) {
   print @usage;
   exit(0);
}
if (defined $config_file) {
   if (!-e $config_file){
   	print "\nThe configuration file $config_file does not exist!\n\n";
        print @usage;
	exit;
   }
}else{
    print "Please provide a configuration file!\n\n";
    print @usage;
    exit;
}

if (defined $output_dir) {
    if (!-e $output_dir){
    	print "\nThe output directory $output_dir does not exist!\n\n";
        print @usage;
	exit;
    }
}else{
    $output_dir = getcwd;
}
if (!-e "$output_dir/unmapped.1.fq" || !-s "$output_dir/unmapped.1.fq")
{
    if (defined $read_file1){
        if (!-e $read_file1 || !-s $read_file1){
            print "\nRead file1 $read_file1 does not exist!\n\n";
            print @usage;
            exit;
        }else{
            `rm "$output_dir/unmapped.1.fq"` if (-e "$output_dir/unmapped.1.fq");
            `ln -s $read_file1 $output_dir/unmapped.1.fq`;
        }
    }else{
    	print "\nPlease specify read file1 using '--fq1'!\n\n";
        print @usage;
	exit;
    }
}
if (!-e "$output_dir/unmapped.2.fq" || !-s "$output_dir/unmapped.2.fq")
{
    if (defined $read_file2){
        if (!-e $read_file2 || !-s $read_file2){
            print "\nRead file2 $read_file2 does not exist!\n\n";
            print @usage;
            exit;
        }else{
            `rm "$output_dir/unmapped.2.fq"` if (-e "$output_dir/unmapped.2.fq");
            `ln -s $read_file2 $output_dir/unmapped.2.fq`;
        }
    }else{
        $paired = 0;
    }
}

if(lc($markdup) eq 'no' || lc($markdup) eq 'n'){
    $markdup = 'n' ;
}else{
    $markdup = 'y' ;
}


my $config = new();
$config->read($config_file);

my $bwa_bin        = $config->get_value("bwa_bin");
my $script_dir     = "$FindBin::Bin";
my $icorn_dir      = "$script_dir/icorn";
my $thread_no      = $config->get_value("thread_no");
my $virus_database = $config->get_value("virus_database");


if (!-e "$output_dir/results-virus-top1.fa" || !-s "$output_dir/results-virus-top1.fa")
{
    if (defined $virus_sequence){
        if (!-e $virus_sequence || !-s $virus_sequence){

            print "$virus_sequence is not in the file system. Look up virus database...\n";
      	    PulloutTop1Virus($virus_sequence, $virus_database, "$output_dir/results-virus-top1.fa");

            if (!-e "$output_dir/results-virus-top1.fa" || !-s "$output_dir/results-virus-top1.fa"){
                print "\nError: $virus_sequence is not a file and it does not exist in virus database!\n\n";
                print @usage;
                exit;
            }else{
                print "\nFound $virus_sequence in the virus database. It will be used as reference to call viral mutations.\n\n";
            }
        }else{
            `rm "$output_dir/results-virus-top1.fa"` if (-e "$output_dir/results-virus-top1.fa");
            `ln -s $virus_sequence $output_dir/results-virus-top1.fa`;
        }
    }else{
    	print "\nPlease provide virus reference sequence using argument -v!\n\n";
        print @usage;
	exit;
    }
}


my $start_dir = getcwd;
chdir($output_dir);


################################ Consensus virus reference genome #######################################


$virus_sequence = "results-virus-top1.fa";

if ($paired){
    my $sensitivity_level = 1;
    if ($config->has_value("sensitivity_level")) {
        $sensitivity_level = $config->get_value("sensitivity_level");
    }

   my $corrected_ref_file = sprintf "refseq.fa.%d", $sensitivity_level+1;

   if(!-e 'virus-consensus-seq.fa'){

        `mkdir vfix` if (!-e 'vfix');
        chdir("$output_dir/vfix");

        if (!-e $corrected_ref_file || !-s $corrected_ref_file){
            `ln -s ../unmapped.1.fq .` if (!-e "unmapped.1.fq");
            `ln -s ../unmapped.2.fq .` if (!-e "unmapped.2.fq");
            `ln -s ../results-virus-top1.fa refseq.fa`;

            my $align_file = "$output_dir/../step1/alignment.bam";
            die "Can't find BAM file $align_file!\n" if (!-e $align_file);
            my ($g_lower,$g_upper,$g_mean,$g_std) = GetInsertSize($align_file);

            `$icorn_dir/icorn.start.sh refseq.fa 1 $sensitivity_level unmapped.1.fq unmapped.2.fq $g_lower,$g_upper $g_mean $icorn_dir &> /dev/null`;
        }

        if (-e $corrected_ref_file && -s $corrected_ref_file && !-e "../virus-consensus-seq.fa"){
            my $in  = Bio::SeqIO->new(-file => "$corrected_ref_file",        -format => 'Fasta');
            my $out = Bio::SeqIO->new(-file => ">../virus-consensus-seq.fa", -format => 'Fasta');
            while ( my $seq = $in->next_seq() ) {$out->write_seq($seq); }
        }

    	chdir($output_dir);
    }

    $virus_sequence = 'virus-consensus-seq.fa';
   `cp vfix/$corrected_ref_file.PerBase.stats virus-consensus-variant.txt` if (-e "vfix/$corrected_ref_file.PerBase.stats" && !-e "virus-consensus-variant.txt");
}


################################ Align unmapped reads to virus.fa #######################################


if (!-e '$virus_sequence.bwt'){
    print "Generate BWA index of the virus reference...\n\n";
    system("$bwa_bin index $virus_sequence 2> /dev/null");
}

if (!-e 'alignment.sorted.bam'){
   print "Align reads to virus genome...\n\n";

   my $RG = '@RG\tID:N\tSM:N\tLB:illumina\tPL:illumina';
   my $unmappedfq1 = 'unmapped.1.fq';
   my $unmappedfq2 = 'unmapped.2.fq';

   system("$bwa_bin aln -t $thread_no $virus_sequence $unmappedfq1 > unmapped.1.sai 2> /dev/null");

   if ($paired){
       system("$bwa_bin aln -t $thread_no $virus_sequence $unmappedfq2 > unmapped.2.sai 2> /dev/null");
       system("$bwa_bin sampe -r '$RG' $virus_sequence unmapped.1.sai unmapped.2.sai $unmappedfq1 $unmappedfq2 > alignment.sam 2> /dev/null");
   }else{
       system("$bwa_bin samse -r '$RG' $virus_sequence unmapped.1.sai $unmappedfq1 > alignment.sam 2> /dev/null");
   }

   print "Convert SAM file to BAM format...\n\n";
   `samtools view -bS alignment.sam -o alignment.bam 2> /dev/null`;

   print "Sort the BAM file...\n\n";
   `samtools sort   alignment.bam   alignment.sorted`;
   `samtools index  alignment.sorted.bam`;
   `rm alignment.sam alignment.bam`;
}

################################ Mark duplicates (optional) #######################################


my $analysisBAM = 'alignment.sorted.bam';
if ($markdup eq 'y'){
   print "Mark duplicate reads...\n\n";
   if (!-e 'alignment.sorted.md.bam'){
       if ($paired){
           `samtools rmdup    alignment.sorted.bam alignment.sorted.md.bam`;
       }else{
           `samtools rmdup -s alignment.sorted.bam alignment.sorted.md.bam`;
       }
       `samtools index alignment.sorted.md.bam`;
   }
   $analysisBAM = 'alignment.sorted.md.bam';
}

################################ Realignment #######################################


my $GATK_bin = "$script_dir/bin/GenomeAnalysisTK.jar";
my $picard_bin = "$script_dir/bin/CreateSequenceDictionary.jar";

if (!-e 'alignment.sorted.realigned.bam'){
    print "Perform realignment...\n\n";

    my $dict_file = $virus_sequence;
    $dict_file =~ s/fa/dict/;

    `rm $dict_file` if(-e $dict_file && !-s $dict_file);
    if (!-e $dict_file){
    	`java -jar $picard_bin REFERENCE=$virus_sequence OUTPUT=$dict_file TRUNCATE_NAMES_AT_WHITESPACE=true NUM_SEQUENCES=2147483647 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false`;
    }

    `java -jar $GATK_bin -T RealignerTargetCreator -I $analysisBAM -R $virus_sequence -o $analysisBAM.intervals`;
    `java -jar $GATK_bin -T IndelRealigner --targetIntervals $analysisBAM.intervals -I $analysisBAM -R $virus_sequence --out alignment.sorted.realigned.bam`;
}

die "\nError: GATK failed to realign the alignment file!\n" if (!-e 'alignment.sorted.realigned.bam');


################################ Variant calling #######################################


print "Calling variants...\n\n";

if (!-e 'mutation-raw.vcf' || !-s 'mutation-raw.vcf'){
   #`java -jar $GATK_bin -T UnifiedGenotyper  -R $virus_sequence -glm BOTH -I alignment.sorted.realigned.bam -o raw-mutation.vcf -metrics snps.metrics -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1000 -A DepthOfCoverage -A AlleleBalance`;

   my $dcov = 500;
   $dcov = 10000 if($markdup eq 'n');

   `java -jar $GATK_bin -T UnifiedGenotyper  -R $virus_sequence -glm BOTH -I alignment.sorted.realigned.bam -o mutation-raw.vcf -metrics mutation-raw.metrics -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov $dcov -A AlleleBalance`;
}

if (!-e 'mutation.vcf' || !-s 'mutation.vcf'){
   `java -jar $GATK_bin -T VariantFiltration -R $virus_sequence --variant mutation-raw.vcf -o mutation.vcf --clusterWindowSize 10 --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\" --filterExpression \"DP < 5 \" --filterName \"LowCoverage\" --filterExpression \"QUAL < 30.0 \" --filterName \"VeryLowQual\" --filterExpression \"QUAL > 30.0 && QUAL < 50.0 \" --filterName \"LowQual\" --filterExpression \"QD < 1.5 \" --filterName \"LowQD\" --filterExpression \"SB > -10.0 \" --filterName \"StrandBias\"`;
}

print "Finished variant calling!\n\n";

chdir($start_dir);