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
#  Last update: 06/19/2014
#
#  Contact:       Zhongming Zhao: zhongming.zhao@vanderbilt.edu
#                 Qingguo Wang:   qingguo.wang@vanderbilt.edu
#
#  Organization:  Bioinformatics and Systems Medicine Laboratory
#                 Vanderbilt University Medical Center
#                 Nashville, Tennessee, USA
#
#  VirusFinder.pl is part of VirusFinder. It is the interface of the software.
#
###################################################################################
# Changelog:
# 03/12/2013   Software release
# 04/16/2013   Code added to handle Perl module installation errors
# 07/19/2013   Code added to read Perl module paths specified in user's command line
# 07/25/2013   Version 1.1 released. V1.1 provided two detection modes:
#               (1) normal: V1.1 in this mode does exactly what V1.0 does.
#               (2) sensitive: this mode is unique for V1.1.
#              The variable 'detection_mode' in configuration file specifies the mode to run
# 07/31/2013   The default value of the vairable $flank_region_size changed to 4000
# 10/18/2013   blat V.34 was included in VirusFinder 1.2 release based on user's feedback;
#              New 1.2 version also reports virus insertion sites detected in mitochondria.
# 11/11/2013   A bug fixed. The bug caused VirusFinder not to report sites with >500 coverage
# 12/18/2013   Code modified to report the confidence level, either high or low, of virus
#              integration detection.
# 12/26/2013   The command line interface of other scripts updated.
# 01/17/2014   This program now supports single-end sequencing reads.
# 02/07/2014   Code added to call snps and indels in the virus genome.
# 02/19/2014   Release of VirusFinder 2.
# 03/17/2014   Code added to allow turning off virus integration/mutation detection
######################################################################################

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use Cwd;
use Bio::SeqIO;
use FindBin;
use lib "$FindBin::Bin";
use Mosaic;

my @usage;
push @usage, "Program: VirusFinder, a tool for characterizing intra-host viruses through NGS data.\n";
push @usage, "Version: 2 (06/19/2014)\n\n";
push @usage, "Usage: VirusFinder.pl -c <configuration file> [options]\n\n";
push @usage, "Options:\n";
push @usage, "  -h, --help     Displays this information\n";
push @usage, "  -c, --config   Configuration file <required>\n";
push @usage, "  -v, --virus    The sequence file (in fasta format) of the virus. Not required\n";
push @usage, "  -o, --output   The directory to store software output, default is current working directory\n";
push @usage, "  -m, --markdup  Mark duplicate reads. Accepted inputs: y[es], n[o]. Default value is y[es].\n";
push @usage, "                 Duplicate reads will not be used for variant calling. For ultra-deep amplicon\n";
push @usage, "                 sequencing data, 'no' should be used for this argument.\n\n";


my $help;
my $config_file;
my $virus_sequence;
my $output_dir;
my $mode = 'normal';
my $markdup = 'y';


GetOptions
(
 'h|help|?'    => \$help,
 'config=s'    => \$config_file,
 'virus=s'     => \$virus_sequence,
 'output=s'    => \$output_dir,
 'markdup=s'   => \$markdup,
);


if ($help) {
   print @usage;
   exit(0);
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
if (defined $virus_sequence) {
    if (!-e $virus_sequence){
	print "The virus sequence $virus_sequence does not exist!\n\n";
        print @usage;
	exit;
    }
}

if(lc($markdup) eq 'no' || lc($markdup) eq 'n'){
    $markdup = 'n' ;
}else{
    $markdup = 'y' ;
}


my $config = new();
$config->read($config_file);


# read configuration file

my $script_dir                = "$FindBin::Bin";
my $icorn_dir                 = "$script_dir/icorn";
my $preprocess_script         = "$script_dir/preprocess.pl";
my $detect_virus_script       = "$script_dir/detect_virus.pl";
my $detect_integration_script = "$script_dir/detect_integration.pl";
my $detect_mutation_script    = "$script_dir/detect_mutation.pl";
my $blastn_index_human        = $config->get_value("blastn_index_human");


my $mailto;
if ($config->has_value("mailto")) {
    $mailto = $config->get_value("mailto");
}
if ($config->has_value("detection_mode")) {
    $mode = lc($config->get_value("detection_mode"));
    if($mode ne 'sensitive' && $mode ne 'normal'){
	print "Please correct the configuration file: an unrecognized value was set to the parameter 'detection_mode'!\n\n";
	exit;
    }
}

######################################### Step1 Preprocessing ##########################################

system "date";

print "step 1 preprocessing...\n";

if (!-e "$output_dir/step1"){
    `mkdir $output_dir/step1`;
}
my $ILIBs = GetIdir();
`perl $ILIBs $preprocess_script -c $config_file -o $output_dir/step1`;

if (!-e "$output_dir/step1/unmapped.1.fa" || !-s "$output_dir/step1/unmapped.1.fa"){
    die "Step 1 terminated abnormally! Please delete intermediate files and try again.\n";
}

######################################### Step2 detect viruses #########################################

if (! defined $virus_sequence) {

    print "step 2 detect virus...\n";

    if (!-e "$output_dir/step2"){
        `mkdir $output_dir/step2`;
    }

    `ln -s $output_dir/step1/unmapped.1.fa  $output_dir/step2/` if (!-e "$output_dir/step2/unmapped.1.fa");
    `ln -s $output_dir/step1/unmapped.2.fa  $output_dir/step2/` if (-e  "$output_dir/step1/unmapped.2.fa" && !-e "$output_dir/step2/unmapped.2.fa");

    system("perl $ILIBs $detect_virus_script -c $config_file -o $output_dir/step2");

    if (!-e "$output_dir/step2/results-virus.txt" ||
    	!-e "$output_dir/step2/results-virus-top1.fa" ||
        !-s "$output_dir/step2/results-virus-top1.fa"){
        print "\nFailed to detect viruses. \n\n";
        &mailme(0);
        system "date";
        exit;
    }

    `cp $output_dir/step2/results-virus.txt           $output_dir/virus.txt`;
    `cp $output_dir/step2/results-virus-list.txt      $output_dir/virus-list.txt`;
    `cp $output_dir/step2/results-contig.txt          $output_dir/contig.txt`;
    `cp $output_dir/step2/results-novel-contig-2.fa   $output_dir/novel-contig.fa` if (-e '$output_dir/step2/results-novel-contig-2.fa');

    if (-s "$output_dir/step2/results-virus-top1.fa"){
        my $ret = `awk '{if(NR==2)print\$1}' $output_dir/virus.txt`;
    	print "\nVirus detected: $ret\n\n";
    }
}

################################# Step3 detect virus integration sites #################################

my $detect_integration = 'y';
if($config->has_value("detect_integration")){
   $detect_integration = substr(lc($config->get_value("detect_integration")),0,1);
}

if ($detect_integration ne 'n'){

    print "step 3 detect virus integration sites...\n";

    DetectIntegration();

    if (-s "$output_dir/integration-sites.txt"){
        print "\nVirus integration detected! Please check $output_dir/integration-sites.txt for detailed insertion sites.\n\n";
    }
}


################################### Step4 detect viral mutations ###################################

my $detect_mutation = 'y';
if($config->has_value("detect_mutation")){
   $detect_mutation = substr(lc($config->get_value("detect_mutation")),0,1);
}

if ($detect_mutation ne 'n'){

    print "step 4 detect viral mutations...\n";

    DetectMutation();
}


&mailme(1);
system "date";
print "Done!\n";


############ Detect virus integration sites ###############


sub DetectIntegration {

    if (!-e "$output_dir/step3"){
        `mkdir $output_dir/step3`;
    }

    if (defined $virus_sequence) {

        if (!-e "$output_dir/step3/results-virus-top1.fa"){
            #`ln -s $virus_sequence  $output_dir/step3/results-virus-top1.fa`;
            open (IN,  $virus_sequence);
            open (OUT, ">$output_dir/step3/results-virus-top1.fa");
            while(<IN>){
                if (/>/){
                    print OUT ">chrVirus\n";
                }else{
                    print OUT $_;
                }
            }
            close(IN);
            close(OUT);
        }

    }else{

        if (!-e "$output_dir/step3/results-virus-top1.fa"){
            `ln -s $output_dir/step2/results-virus-top1.fa $output_dir/step3/`;
        }

    }


    `ln -s $output_dir/step1/unmapped.1.fq  $output_dir/step3/` if (!-e "$output_dir/step3/unmapped.1.fq");
    `ln -s $output_dir/step1/unmapped.2.fq  $output_dir/step3/` if (-e  "$output_dir/step1/unmapped.2.fq" && !-e "$output_dir/step3/unmapped.2.fq");
    `ln -s $blastn_index_human.fa $output_dir/step3/hg19.fa`    if (!-e "$output_dir/step3/hg19.fa");

    if ($mode eq 'sensitive' && -e "$output_dir/step3/unmapped.2.fq"){
        RunSensitiveMode();
        `cp $output_dir/step3/virus-corrected-seq.fa $output_dir/virus-consensus-seq.fa`;
    }else{
        system("perl $ILIBs $detect_integration_script -c $config_file -o $output_dir/step3");
    }

    `cp $output_dir/step3/results-virus-loci.txt $output_dir/integration-sites.txt`;

}

############ Detect viral mutations ###############


sub DetectMutation {

    if (!-e "$output_dir/step4"){
        `mkdir $output_dir/step4`;
    }
    if (defined $virus_sequence) {

        if (!-e "$output_dir/step4/results-virus-top1.fa"){
            #`ln -s $virus_sequence  $output_dir/step3/results-virus-top1.fa`;
            open (IN,  $virus_sequence);
            open (OUT, ">$output_dir/step4/results-virus-top1.fa");
            while(<IN>){
                if (/>/){
                    print OUT ">chrVirus\n";
                }else{
                    print OUT $_;
                }
            }
            close(IN);
            close(OUT);
        }

    }else{

        if (!-e "$output_dir/step4/results-virus-top1.fa"){
            `ln -s $output_dir/step2/results-virus-top1.fa $output_dir/step4/`;
        }

    }

    `ln -s $output_dir/step1/unmapped.1.fq  $output_dir/step4/` if (!-e "$output_dir/step4/unmapped.1.fq");
    `ln -s $output_dir/step1/unmapped.2.fq  $output_dir/step4/` if (-e  "$output_dir/step1/unmapped.2.fq" && !-e "$output_dir/step4/unmapped.2.fq");

    if ($mode eq 'sensitive' && -e "$output_dir/step4/unmapped.2.fq" && -e "$output_dir/step3/a-vfix"){
        `ln -s $output_dir/step3/a-vfix $output_dir/step4/vfix`;
        `ln -s $output_dir/step3/virus-corrected-seq.fa $output_dir/step4/virus-consensus-seq.fa` if (!-e "$output_dir/step4/virus-consensus-seq.fa");
    }

    system("perl $ILIBs $detect_mutation_script -c $config_file -o $output_dir/step4 -m $markdup");

    if (-e "$output_dir/step4/mutation.vcf"){
       `cp  $output_dir/step4/mutation.vcf $output_dir/viral-mutation.vcf`;
        print "Please check viral-mutation.vcf for called SNPs and indels.\n\n";
    }
}

############ Sensitive detection mode ###############


sub RunSensitiveMode {

    my $sensitivity_level = 1;
    if ($config->has_value("sensitivity_level")) {
        $sensitivity_level = $config->get_value("sensitivity_level");
    }

    my $flank_region_size = 4000;
    if ($config->has_value("flank_region_size")) {
        $flank_region_size = $config->get_value("flank_region_size");
    }

    my $start_dir = getcwd;

    # Get sequencing insert size

    my $align_file = "$output_dir/step1/alignment.bam";
    die "Can't find BAM file $align_file!\n" if (!-e $align_file);


    my ($g_lower,$g_upper,$g_mean,$g_std) = GetInsertSize($align_file);
    #print "\nInsert size distribution: low:$g_lower; high:$g_upper; mean:$g_mean\n\n";


    # Correct reference virus sequence
    print "\nStep3-a-vfix Correcting virus reference sequence...\n\n";

    my $work_dir = "$output_dir/step3/a-vfix";
    if (!-e $work_dir){
       `mkdir $work_dir`;
    }
    chdir($work_dir);

    my $corrected_ref_file = sprintf "refseq.fa.%d", $sensitivity_level+1;
    if (!-e $corrected_ref_file || !-s $corrected_ref_file){
        `ln -s ../unmapped.1.fq .` if (!-e "unmapped.1.fq");
        `ln -s ../unmapped.2.fq .` if (!-e "unmapped.2.fq");
        `ln -s ../results-virus-top1.fa refseq.fa`;
        `$icorn_dir/icorn.start.sh refseq.fa 1 $sensitivity_level unmapped.1.fq unmapped.2.fq $g_lower,$g_upper $g_mean $icorn_dir &> /dev/null`;
    }

    if (-e $corrected_ref_file && -s $corrected_ref_file && !-e "../virus-corrected-seq.fa"){
        my $in  = Bio::SeqIO->new(-file => "$corrected_ref_file",        -format => 'Fasta');
        my $out = Bio::SeqIO->new(-file => ">../virus-corrected-seq.fa", -format => 'Fasta');
        while ( my $seq = $in->next_seq() ) {$out->write_seq($seq); }
    }

    `cp $corrected_ref_file.PerBase.stats ../virus-snp-indel.txt` if (-e "$corrected_ref_file.PerBase.stats" && !-e "../virus-snp-indel.txt");

    # Run SVDetect to obtain virus-harboring regions in human reference
    print "\nStep3-b-align Detect virus integration using corrected virus sequence...\n\n";

    $work_dir = "$output_dir/step3/b-align";
    if (!-e $work_dir){
       `mkdir $work_dir`;
    }

    my $SVDetect_outfile = "$work_dir/SVDetect/results-virus-loci.txt";
    #if (!-e $SVDetect_outfile || !-s $SVDetect_outfile){
        chdir($work_dir);
        `ln -s ../unmapped.1.fq .` if (!-e "unmapped.1.fq");
        `ln -s ../unmapped.2.fq .` if (!-e "unmapped.2.fq");

        if (!-e "results-virus-top1.fa"){
            if (!-e "../a-vfix/$corrected_ref_file" || !-s "../a-vfix/$corrected_ref_file"){
               `ln -s ../results-virus-top1.fa .`;
            }else{
#                my $in  = Bio::SeqIO->new(-file => "../a-vfix/$corrected_ref_file", -format => 'Fasta');
#                my $out = Bio::SeqIO->new(-file => ">results-virus-top1.fa",        -format => 'Fasta');
#                while ( my $seq = $in->next_seq() ) {$out->write_seq($seq); }
               `ln -s ../virus-corrected-seq.fa ./results-virus-top1.fa`;
            }
        }

        `ln -s $blastn_index_human.fa hg19.fa` if (!-e "hg19.fa");
        chdir($start_dir);
        `perl $ILIBs $detect_integration_script -c $config_file -o $work_dir`;

        if (!-e $SVDetect_outfile || !-s $SVDetect_outfile){
	    `cp $work_dir/results-virus-loci.txt $output_dir/step3`;
            return;
        }
    #}

    if (-e "$output_dir/step3/b-align/crest/results.predSV.txt" && -s "$output_dir/step3/b-align/crest/results.predSV.txt"){
    	`awk \'{if(\$1==\"chrVirus\" && \$5==\"chrVirus\")print \$0}\' $output_dir/step3/b-align/crest/results.predSV.txt > $output_dir/step3/virus-SV.txt`;
    }


    my ($chrs, $lowbounds, $upbounds) = GetRefSubSeq($SVDetect_outfile, $flank_region_size);

    ($chrs, $lowbounds, $upbounds) = FilterRefSubSeq(\@$chrs,
    				                     \@$lowbounds,
                                                     \@$upbounds,
                                                     "$work_dir/crest/results-virus-loci.txt",
                                                     "$work_dir/crest/alignment.sorted.bam.sclip.txt");

    if (scalar(@$chrs) == 0){
    	`cp $work_dir/results-virus-loci.txt $output_dir/step3`;
        chdir($start_dir);
        return;
    }


    # Correct virus-harboring sub sequence of human reference
    print "\nStep3-c-rfix Correcting human reference sequence...\n\n";

    $work_dir = "$output_dir/step3/c-rfix";
    if (!-e $work_dir){
       `mkdir $work_dir`;
    }
    chdir($work_dir);

    if (!-e "$work_dir/$corrected_ref_file" || !-s "$work_dir/$corrected_ref_file"){

        if (!-e 'refseq.fa'){
            WriteRefSubSeq(\@$chrs,\@$lowbounds,\@$upbounds,"$blastn_index_human.fa",'refseq.fa');
        }

        # Recuit aligned reads from the BAM file
        `ln -s ../unmapped.1.fq .` if (!-e "unmapped.1.fq");
        `ln -s ../unmapped.2.fq .` if (!-e "unmapped.2.fq");

        my $unmappedfq1 = 'unmapped.1+.fq';
        my $unmappedfq2 = 'unmapped.2+.fq';
        my $count = 0;
        if  (!-e 'unmapped.1+.fq' || !-e 'unmapped.2+.fq' || !-s 'unmapped.1+.fq' || !-s 'unmapped.2+.fq'){

            $count = ExtractMappedReads($align_file,'unmapped.1+.fq','unmapped.2+.fq',\@$chrs,\@$lowbounds,\@$upbounds);
            if ($count > 0){
               `cat ../unmapped.1.fq >> unmapped.1+.fq`;
               `cat ../unmapped.2.fq >> unmapped.2+.fq`;
            }
        }

        if (!-e 'unmapped.1+.fq' || !-e 'unmapped.2+.fq' || !-s 'unmapped.1+.fq' || !-s 'unmapped.2+.fq'){
           `$icorn_dir/icorn.start.sh refseq.fa 1 $sensitivity_level unmapped.1.fq  unmapped.2.fq  $g_lower,$g_upper $g_mean $icorn_dir &> /dev/null`;
        }else{
           `$icorn_dir/icorn.start.sh refseq.fa 1 $sensitivity_level unmapped.1+.fq unmapped.2+.fq $g_lower,$g_upper $g_mean $icorn_dir &> /dev/null`;
        }
    }


    ## Detect virus integration sites
    print "\nstep3-d-align Detect virus integration using corrected reference sequence...\n\n";

    $work_dir = "$output_dir/step3/d-align";
    if (!-e $work_dir){
       `mkdir $work_dir`;
    }

    my $crest_outfile = "$work_dir/crest/results-virus-loci.txt";
    if (!-e $crest_outfile || !-s $crest_outfile){
        chdir($work_dir);
        `ln -s ../unmapped.1.fq .` if (!-e "unmapped.1.fq");
        `ln -s ../unmapped.2.fq .` if (!-e "unmapped.2.fq");

        if (!-e "results-virus-top1.fa"){
            `ln -s ../b-align/results-virus-top1.fa .`;
        }

        if (!-e "hg19.fa"){
            my $in  = Bio::SeqIO->new(-file => "../c-rfix/$corrected_ref_file", -format => 'Fasta');
            my $out = Bio::SeqIO->new(-file => ">hg19.fa",                      -format => 'Fasta');
            while (my $seq = $in->next_seq() ) {$out->write_seq($seq);}
        }

        chdir($start_dir);
        `perl $ILIBs $detect_integration_script -c $config_file -o $work_dir -m sensitive`;
    }

    if (!-e $crest_outfile || !-s $crest_outfile){
    	`cp $output_dir/step3/b-align/results-virus-loci.txt $output_dir/step3`;
        chdir($start_dir);
        return;
    }

    ####### Pullout virus insertion sites

    ## Read integration loci detected in step3/d-align
    open (IN, $crest_outfile);
    my @lines = <IN>;
    close(IN);

    ## Read the coverage of each loci
    my %cov_hash;
    for (my $i=1; $i<scalar(@lines); $i++){
        my @data = split(/\t/, $lines[$i]);

        my $cov = eval $data[6];
        my $key = $data[0];
        $key = $data[3] if ($data[0] eq 'chrVirus');

        if(exists $cov_hash{$key}){
            $cov_hash{$key} = $cov if ($cov_hash{$key} < $cov);
        }else{
            $cov_hash{$key} = $cov;
        }
    }

    ## Combine the two detections
    my $consensus_output = "$output_dir/step3/results-virus-loci.txt";
    open (OUT, ">$consensus_output");
    print OUT $lines[0];

    ## Write integration loci detected in step3-b-align
    my $step3_b_file = "$output_dir/step3/b-align/results-virus-loci.txt";
    if (-e $step3_b_file && -s $step3_b_file){
    	my $flag = 0;
        open (IN, $step3_b_file);
        while(<IN>){
            if ($flag==0){
                $flag++;
            }else{
                print OUT $_;
            }
        }
        close(IN);
    }

    ## Write integration loci detected in step3-d-align
    for (my $i=1; $i<scalar(@lines); $i++){
        my @data = split(/\t/, $lines[$i]);
        my $cov = eval $data[6];
        if (scalar(@data)<8){
           chomp $data[6];
           push @data, "\n";
        }

        if ($data[0] ne 'chrVirus'){
            next if ($cov_hash{$data[0]}>$cov);
            my @temp = split(/_/, $data[0]);
            my $pos = $temp[1]+$data[1]-1;
            print OUT "$temp[0]\t$pos\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]";
        }else{
            next if ($cov_hash{$data[3]}>$cov);
            my @temp = split(/_/, $data[3]);
            my $pos = $temp[1]+$data[4]-1;
            print OUT "$data[0]\t$data[1]\t$data[2]\t$temp[0]\t$pos\t$data[5]\t$data[6]\t$data[7]";
        }

    }
    close(OUT);
    chdir($start_dir);
}


############################################ e-mail format #############################################

sub mailme {
    return if not defined $mailto;
    my $flag = shift;
    my $text;

    if ($flag == 1){
    	$text = "VirusFinder has succefully analyzed the following data,\n\n";
    }else{
    	$text = "Error occurred when analyzing the following data,\n\n";
    }

    if ($config->has_value("alignment_file")) {
        $text .= $config->get_value("alignment_file")."\n";
    }
    if ($config->has_value("fastq1")) {
        $text .= $config->get_value("fastq1")."\n";
    }
    if ($config->has_value("fastq2")) {
        $text .= $config->get_value("fastq2")."\n";
    }
    if ($flag == 1){
    	system "echo '$text' | mail -s 'VirusFinder complete notice' $mailto";
    }else{
    	system "echo '$text' | mail -s 'VirusFinder error notice' $mailto";
    }
}