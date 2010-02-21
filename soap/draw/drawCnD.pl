#!/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);

my $coverage_bin = "$Bin/subBin/soap.coverage";
my $plot_bin = "perl $Bin/subBin/plotDepthDistribution.pl";
my $plot_chr_bin = "perl $Bin/subBin/plotChrDepthDistribution.pl";
my $call_ngap_bin = "perl $Bin/subBin/callNgap.pl";

my ($pro_dir, $ref_fa, $pro_name, $out_dir, $input_soap_list);
my ($is_effective, $ngap_file, $gc_file, $chrorder_file);
my ($step, $mem);
my ($help);
GetOptions (
	"proDir:s" => \$pro_dir,
	"il:s" => \$input_soap_list,
	"proName:s" => \$pro_name,
	"ref:s" => \$ref_fa,
	"outDir:s" => \$out_dir,
	"effective:i" => \$is_effective,
	"ngapFile:s" => \$ngap_file,
	"gcFile:s" => \$gc_file,
	"chrid:s" => \$chrorder_file,
	"step:i" => \$step,
	"mem:s" => \$mem,
	"help" => \$help
);
if (!$pro_dir || !$pro_name || !$ref_fa || $help) {
        print "\n";
	print "Mandatory:\n";
	print "\t-il\t*.soap/*.single file list\n";
	print "\t-proDir\tproject directory\n";
	print "\t-proName\tproject name\n";
	print "\t-ref\treference fa file\n";
	print "\t-outDir\t[proDir]\n";
        print "\t-chrid\tthe file contains the chromosome id[/share/raid010/resequencing/resequencing/tmp/pub/Genome/Rice/Genome_IRGSP/IRGSP.chrorder]\n";
	print "Option:\n";
	print "\t-gcFile\tthe file contains the gc percentage of each chromosome(file format: chr total_genome_size effective_genome_size gc_size gc_percentage)if not provide this file, one gc file will be generated automatically by this program.\n";
	print "\t-effective\t0:not; 1:yes[0]\n";
        print "\t-ngapFile\tngap file\n";
	print "\t-mem\tmemory to run soap.coverage[6g for human]\n";
	print "\t-step\t3 steps: coverage(1)->read depth(2)->draw depth distribution(3)\t[1]\n";
	print "Example:\n\tperl $0 -proDir /share/tmp/resequencing/human/Arab/ -proName Arab -ref /share/tmp/pub/Genome/Human/human.fa -outDir /share/tmp/resequencing/human/Arab/ -effective 1 -step 2\n";
	print "\tperl $0 -proDir /share/raid11/zhanghao/software/postProcess/test/Rice/IRGSP -il soap.l -ref /share/raid1/database/BGI/rice/IRGSP_chromosomes_build04.fa -outDir /share/raid11/zhanghao/software/postProcess/test/Rice/IRGSP -ngapFile /share/raid11/zhanghao/software/postProcess/test/Rice/IRGSP/IRGSP.ngap -mem 4g -effective 1  -step 3 -proName IRGSP\n";
	print "\tperl $0 -proDir /share/raid11/zhanghao/software/postProcess/test/Rice/9311 -il soap.l -ref /share/raid1/database/BGI/rice/9311_main_chromosomes.fa -outDir /share/raid11/zhanghao/software/postProcess/test/Rice/9311 -ngapFile /share/raid11/zhanghao/software/postProcess/test/Rice/9311/9311.ngap -mem 4g -effective 1 -step 3 -proName 9311\n";
	print "\tAuther: Hao Zhang\tTime: 23:22 25/05/2009\n";
	print " modify by HuXuesong @ 2010\n";

        exit 0;
}
$out_dir ||= $pro_dir;
mkdir $out_dir unless (-d $out_dir);
$input_soap_list ||= '';
$ngap_file ||= '';	#"$Bin/subBin/pub_data/human.ngap";
$gc_file ||= '';
$chrorder_file ||= "$Bin/subBin/pub_data/human.chrnum";
$is_effective ||= 0;
$mem ||= "6g";
$step ||= 1;
#my $pe_dir = "$pro_dir/PE";
my $pe_dir = "$pro_dir/soap/PE";
my $depth_dir = "$out_dir/depth";
mkdir $depth_dir unless (-d $depth_dir);
my $coverage_dir = "$out_dir/coverage";
mkdir $coverage_dir unless (-d $coverage_dir);
my $distribution_dir = "$depth_dir/distribution";
mkdir $distribution_dir unless (-d $distribution_dir);

my $wait_second = 30;

my %hDepth = ();        #{depth} => num
my %hChrDepth = ();     #{chr}{depth} => num
my %hChrGeSize = ();    #{chr} => size
my %hChrMapBase = ();   #{chr} => basenum
my %hChrEffLen = ();	#{chr} => effective length

my $total_list = $input_soap_list;
my $total_base_count = 0;


if (1 == $step) {
	soapCoverage($total_list, \$is_effective, \$out_dir, \$coverage_dir, \$depth_dir);
}
elsif (2 == $step) {
        readDepth(\$total_base_count);
	calculateMeanDepthAndPoisson(\$total_base_count);
	calculateGCpercentage(\$ref_fa, \$gc_file);
	calculateGCpercentageForSoap($total_list);	
	drawDepthDistribution();

	my $gc_file2="$out_dir/$pro_name.soap_gc.info";

	drawChrDepthDistribution2($gc_file2);
}
elsif (3 == $step) {
	#calculateGCpercentage(\$ref_fa, \$gc_file);
	#calculateGCpercentageForSoap($total_list);
	my $gc_file2="$out_dir/$pro_name.soap_gc.info";
	print "$gc_file2\n";

	drawDepthDistribution();
	drawChrDepthDistribution2($gc_file2);
}

elsif(4== $step)
{

	readDepth2(\$total_base_count);
       	calculateMeanDepthAndPoisson(\$total_base_count);
	drawDepthDistribution();

}

sub soapCoverage {
	my ($total_list, $is_effective, $out_dir, $coverage_dir, $depth_dir) = @_;
	######generate list######
	my $soap_list = '';
	unless (-f $total_list) {	#if there is not soap.list available, create one list
		$soap_list = "$$depth_dir/soap.list";
		my $single_list = "$$depth_dir/single.list";
		$total_list = "$$depth_dir/file.list";
		system ("cp $soap_list $soap_list.old") if (-f $soap_list);
	        system ("cp $single_list $single_list.old") if (-f $single_list);
        	#system ("cp $$total_list $$total_list.old") if (-f $$total_list);
		#`ls $pe_dir/*.soap >$soap_list`;
		#`ls $pe_dir/*.single >$single_list`;
		#`cat $soap_list $single_list >$total_list`;

		system "ls $pe_dir/*.soap >$soap_list";
		system "ls $pe_dir/*.single >$single_list";
		system "cat $soap_list $single_list >$total_list";
	}


	if($$is_effective) {
		######generate ngap list######
		unless (-f $ngap_file) {
			$ngap_file = "$$depth_dir/$pro_name.ngap";
			unless (-f $ngap_file) {
				my $call_ngap_command = "$call_ngap_bin $ref_fa 1 $ngap_file";
				print "$call_ngap_command\n";
				#`$call_ngap_command`;
				system "$call_ngap_command";
			}
		}

		open SH, ">$$out_dir/run_total_coverage.sh" or die "$!";
        	print SH "#! /bin/sh\n#\$ -S /bin/sh\nexport LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/lib;\n$coverage_bin -cvg -refsingle $ref_fa -il $total_list -o $$coverage_dir/total_coverage.info -depthsingle $$depth_dir/total_depthsingle -addn $ngap_file";
        	close SH;
	}
	else {
                open SH, ">$$out_dir/run_total_coverage.sh" or die "$!";
                print SH "#! /bin/sh\n#\$ -S /bin/sh\nexport LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/lib;\n$coverage_bin -cvg -refsingle $ref_fa -il $total_list -o $$coverage_dir/total_coverage.info -depthsingle $$depth_dir/total_depthsingle";
                close SH;
	}

	my %jobIDs = ();
	my @coverage_shes = `find $$out_dir -name "run_total_coverage.sh"`;
	foreach my $sh (@coverage_shes) {
		my $qsub = `qsub -cwd -l vf=$mem $sh`;
		my $jobID = (split /\s+/, $qsub)[2];
		my $local_time = localtime();
		print "$local_time\nsubmit $jobID\n";
		$jobIDs{$jobID} = 1;
	}

	waitJobsDone2(%jobIDs);
}

sub waitJobsDone2 {
	my (%jobIDs) = @_;
	while(keys %jobIDs){
	foreach my $jobid (keys %jobIDs) {
		my $qstat = `qstat -j $jobid`;

		my %runningJobs = ();
		chomp $qstat;

		my @qstat_lines = split /\n/, $qstat;
		chomp @qstat_lines;

	if(defined($qstat_lines[1]))
	{

#	print "test : $qstat_lines[1]\n";	
		my $job = (split /\s+/,$qstat_lines[1])[1];

#		print "job: $job\n";

		if($job eq $jobid){
				$runningJobs{$job} = 1;
		}
	}	
		foreach my $jobID(keys %jobIDs){
			if(! exists $runningJobs{$jobID}){
				my $local_time = localtime();
				print "\n$local_time\nFinish jobid $jobID\n";
				delete $jobIDs{$jobID};
			}
		}
		sleep $wait_second;
	}

   } #end while
=b  # old 
	readDepth(\$total_base_count);
	calculateMeanDepthAndPoisson(\$total_base_count);
	calculateGCpercentage(\$ref_fa, \$gc_file);
	calculateGCpercentageForSoap($total_list);
	drawDepthDistribution();
	drawChrDepthDistribution();
=cut

	readDepth(\$total_base_count);
        calculateMeanDepthAndPoisson(\$total_base_count);
        calculateGCpercentage(\$ref_fa, \$gc_file);
        calculateGCpercentageForSoap($total_list);
        drawDepthDistribution();

        my $gc_file2="$out_dir/$pro_name.soap_gc.info";

        drawChrDepthDistribution2($gc_file2);

}





sub readDepth2 {


	my ($total_base_count) = @_;
	my $current_chr = '';
	my $effect_base_count = 0;

	open DEPTH, "$depth_dir/total_depthsingle" or die "$!";
	while (my $line = <DEPTH>) {
        	chomp $line;
	        if ($line =~ /^>/) {
        	        last if ($line =~ /scaffold/i);
                	$line =~ s/^>//;
	                $current_chr = $line;
			my $local_time = localtime();
        	        print "$local_time\tprocess $current_chr\n";
        	}
	        else {
        	        my @depths = split /\s+/, $line;
                	foreach my $depth (@depths) {
                        	$hDepth{$depth} += 1;
	                        $hChrDepth{$current_chr}{$depth} += 1;
        	                $hChrGeSize{$current_chr} += 1;
				if ($depth != 65535) {
                	        	$hChrMapBase{$current_chr} += $depth;
					$hChrEffLen{$current_chr} += 1;
					++$effect_base_count;
				}
                        	++$$total_base_count;
	                }
        	}
	} #end while
	close DEPTH;
	

#	calculateMeanDepthAndPoisson(\$total_base_count);

=b

	print "total base: $$total_base_count\neffective base: $effect_base_count\n";

	open TCHR, ">$distribution_dir/all.info" or die "$!";
	print TCHR "chr\tgenome_size\teffective_size\tmappable_base\n";
	foreach my $chr (keys %hChrDepth) {
		print TCHR "$chr\t$hChrGeSize{$chr}\t$hChrEffLen{$chr}\t$hChrMapBase{$chr}\n";

		open CHR, ">$distribution_dir/$chr.depth" or die "$!";
		foreach my $depth (sort {$a<=>$b} keys %{$hChrDepth{$chr}}) {
                	print CHR "$depth\t$hChrDepth{$chr}{$depth}\n";
		}
		close CHR;
	}
	close TCHR;
=cut

}



sub readDepth {

	my ($total_base_count) = @_;
	my $current_chr = '';
	my $effect_base_count = 0;

	open DEPTH, "$depth_dir/total_depthsingle" or die "$!";
	while (my $line = <DEPTH>) {
        	chomp $line;
	        if ($line =~ /^>/) {
        	        last if ($line =~ /scaffold/i);
                	$line =~ s/^>//;
	                $current_chr = $line;
			my $local_time = localtime();
        	        print "$local_time\tprocess $current_chr\n";
        	}
	        else {
        	        my @depths = split /\s+/, $line;
                	foreach my $depth (@depths) {
                        	$hDepth{$depth} += 1;
	                        $hChrDepth{$current_chr}{$depth} += 1;
        	                $hChrGeSize{$current_chr} += 1;
				if ($depth != 65535) {
                	        	$hChrMapBase{$current_chr} += $depth;
					$hChrEffLen{$current_chr} += 1;
					++$effect_base_count;
				}
                        	++$$total_base_count;
	                }
        	}
	} #end while
	close DEPTH;
	
	print "total base: $$total_base_count\neffective base: $effect_base_count\n";

	open TCHR, ">$distribution_dir/all.info" or die "$!";
	print TCHR "chr\tgenome_size\teffective_size\tmappable_base\n";
	foreach my $chr (keys %hChrDepth) {
		print TCHR "$chr\t$hChrGeSize{$chr}\t$hChrEffLen{$chr}\t$hChrMapBase{$chr}\n";

		open CHR, ">$distribution_dir/$chr.depth" or die "$!";
		foreach my $depth (sort {$a<=>$b} keys %{$hChrDepth{$chr}}) {
                	print CHR "$depth\t$hChrDepth{$chr}{$depth}\n";
		}
		close CHR;
	}
	close TCHR;
}

sub calculateMeanDepthAndPoisson {

	my ($count) = @_;

	#read effective length of reference
	open COVERAGE, "$coverage_dir/total_coverage.info" or die "$!";
	my $ref_base = 0;
	while (my $line = <COVERAGE>) {
		chomp $line;
		if ($line =~ /^Total/) {
			#print "$line\n";
			$ref_base = (split /\:/, $line)[-1];
			print "$ref_base\n";
			last;
		}
	} #end while
	close COVERAGE;

	#count mean depth
	my $sum_mappable_base = 0;
	my $stat_file = "$pro_dir/Statistics.txt";
	if (-f $stat_file) {
		open STAT, "$pro_dir/Statistics.txt" or die "$!";
		my $read_count = 0;
		while (my $line = <STAT>) {
			++$read_count;
			next if ($read_count < 3);
			chomp $line;
			my @stat_info = split /\s+/, $line;
			#next if ($stat_info[0] !~ /\.fq$/);
			$sum_mappable_base += $stat_info[8];
		} #end while
		close STAT;
	}
	else {
		if (-f $total_list) {
			open LIST, "$total_list" or die "$!";
		}
		else {
			open LIST, "$depth_dir/file.list" or die "$!";
		}
		my @soap_files = <LIST>;
		chomp @soap_files;
		close LIST;

		#print "list file $total_list\n";
		#map {print "$_\n"} @soap_files;
		
		foreach my $soap_file (@soap_files) {
			my $local_time = localtime();
			print "$local_time\tprocessing $soap_file\n";
			$sum_mappable_base += `awk '{sum+=\$6}END{print sum}' $soap_file`;
		}
	}
	#print "$sum_mappable_base\n";
	my $mean_depth = $sum_mappable_base / $ref_base;
	print "sum mappable bases: $sum_mappable_base\nreference length: $ref_base\nmean depth: $mean_depth\n";

	open DISTRIBUTION, ">$distribution_dir/$pro_name.distribution.data" or die "$!";
	print DISTRIBUTION "Depth\tnum\tpercentage\tpoisson\n";
	#$hDepth{0} -= (243568484-226441261);
	#$$count -= (243568484-226441261);
	if ($is_effective) {
		$$count -= $hDepth{65535};
	}

print "effective: $$count\n";
	foreach my $depth (sort {$a<=>$b} keys %hDepth) {
        	my $per = $hDepth{$depth} / $$count * 100;

		my $factorial = 1;
	        my $power = 1;
        	for my $index (1..$depth) {
                	$factorial *= $index;
                	$power *= $mean_depth;
        	}
        	my $poisson = $power * exp(0-$mean_depth) /$factorial * 100;
        	print DISTRIBUTION "$depth\t$hDepth{$depth}\t$per\t$poisson\n";
	}
	close DISTRIBUTION;
}

sub calculateGCpercentage {
	my ($fa_file, $gc_file) = @_;
	unless (-f $$gc_file) {
		my $chrName = '';

		my %hChrSeq = ();       #{chr} => sequence
		my @chrNames = ();      #(chr1, chr2,...)

		open FA, $$fa_file or die "$!";
		while (my $line = <FA>) {
        		chomp $line;
        		if ($line =~ /^>(.*)/) {
                		$chrName = $1;
                		push @chrNames, $chrName;
        		} #end if
        		else {
                		$hChrSeq{$chrName} .= $line;
        		} #end else
		} #end while
		close FA;

		open GCOUT, ">$out_dir/$pro_name.ref_gc.info" or die "$!";
		foreach my $chr (@chrNames) {
        		my $chr_len = length $hChrSeq{$chr};
        		my $countN = ($hChrSeq{$chr} =~ s/N/q/ig);
        		my $countGC = ($hChrSeq{$chr} =~ s/[GC]/p/ig);
        		my $effective = $chr_len - $countN;
			my $gc_percentage = $countGC / $effective * 100;
        		print GCOUT "$chr\t$chr_len\t$effective\t$countGC\t$gc_percentage\n";
		} #end foreach
		close GCOUT;

		$$gc_file = "$out_dir/$pro_name.ref_gc.info";

	} #end else
}
sub calculateGCpercentageForSoap{
	my ($soaplist)=@_;
	my @arry_soap=();
	open FI,"$soaplist" or die "can't open $soaplist, in calculateGCpercentageForSoap()!";

print "\n ---- begin statics gc for soap result! \n";

	while(my $line = <FI>)
	{
		chomp $line;

		if($line =~ /\.*.soap/)
		{
			print " soap file: $line\n";
			push @arry_soap, $line;
		}
	}
	

	close FI;


        my %hchr = ();       #{chr} => sequence
	
	foreach my $file (@arry_soap)
	{
		print " stat soapfile $file\n";
                open FA, $file or die "$!";
                while (my $line = <FA>) {
                        chomp $line;
			my @l = split /\s+/, $line;
			if(exists($hchr{$l[7]}))
			{
				${$hchr{$l[7]}}[0] +=$l[1] =~ s/[GC]/p/ig;
				${$hchr{$l[7]}}[1] +=$l[5];
			}
			else
			{
				my $gc_num = $l[1] =~ s/[GC]/p/ig;
				
				push @{$hchr{$l[7]}},$gc_num;
				 push @{$hchr{$l[7]}},$l[5];
	

			}
		



                } #end while
                close FA;

	}
                open GCOUT, ">$out_dir/$pro_name.soap_gc.info" or die "$!";
                foreach my $chr (keys %hchr) {
                        my $chr_len = ${$hchr{$chr}}[1];

                         my $countGC = ${$hchr{$chr}}[0];

                        my $gc_percentage = $countGC / $chr_len * 100;
                        print GCOUT "$chr\t$chr_len\t$countGC\t$gc_percentage\n";
                } #end foreach
                close GCOUT;

}

sub drawChrDepthDistribution2 {

	my ($gc_file2) = @_;

	print "draw chrdepthion !\n";
	
	my @chr_order = ();	#(1..22, "X", "Y");

	my %hChrCoverage = ();	#{chr} => coverage
	## read coverage of each chromosome
	open COVERAGE, " $coverage_dir/total_coverage.info" or die "$!";
	while (my $line = <COVERAGE>) {
        	chomp $line;
        	last if ($line =~ /Overall/i);
        	if ($line =~ /Percentage/i) {
                	my @info = split "\:", $line;
                	#$info[-1] =~ s/\%//;
			#$hChrCoverage{$info[0]} = $info[-1];
			my $covered_length = (split /\//, $info[1])[0];
			$hChrCoverage{$info[0]} = $covered_length;
        	}
	}
	close COVERAGE;


	## count mean depth of each chromosome
	open ALLINFO, "$distribution_dir/all.info" or die "$!";
	my %hMeanDepthChr = ();	#{chr} => mean depth
	my %hMappableChr = ();	#{chr} => mappable bases
	my %hEffChr = ();	#{chr} => effective length
	while (my $line = <ALLINFO>) {
		chomp $line;
		my @info = split /\s+/, $line;
		next if ($info[2] !~ /\d{1,}/);
		$hMeanDepthChr{$info[0]} = $info[3] / $info[2];
		$hMappableChr{$info[0]} = $info[3];
		$hEffChr{$info[0]} = $info[2];
	}
	close ALLINFO;


	## count mode depth of each chromosome
	my %hModeDepthChr = ();	#{chr}{depth} => num;
	my @depth_files = `find $distribution_dir -name "*.depth"`;
	chomp @depth_files;
	foreach my $depth_file (@depth_files) {
		my $basename = basename $depth_file;
		$basename =~ s/\.depth//;
		#$basename =~ s/chromosome//;
		open DEPTH, $depth_file or die "$!";
		while (my $line = <DEPTH>) {
			chomp $line;
			my @info = split /\s+/, $line;
			$hModeDepthChr{$basename}{$info[0]} = $info[1];
		}
		close DEPTH;
	}
	
	## read GC percentage file
	my %hGCperChr = ();	#{chr} => gc percentage
	open GCLIST, $gc_file2 or die "$!";
	while (my $line = <GCLIST>) {
		chomp $line;
		my @info = split /\s+/, $line;
		$hGCperChr{$info[0]} = $info[3];
	}
	close GCLIST;

	open CHRORDER, $chrorder_file or die "$!";
	while (my $line = <CHRORDER>) {
		chomp $line;
		push @chr_order, $line;
	}
	close CHRORDER;

	open CHRDIS, ">$distribution_dir/$pro_name.chrdistribution.data" or die "$!";
	print CHRDIS "chrnum\tchr\tmean_depth\tmode_depth\tgc_percentage\n";
	open STATISTICS, ">$distribution_dir/$pro_name.statistics_by_chr" or die "$!";
	print STATISTICS "chr\teffective_length\tmappable_base\tmean_depth\tmode_depth\tcoverage(%)\n";
	foreach my $chrnum (@chr_order) {
		#my $chr = "canFam2-chr$chrnum";
		#my $chr = "chr$chrnum";
		my $chr = $chrnum;
	#foreach my $chr (sort keys %hModeDepthChr) {
		delete $hModeDepthChr{$chr}{0};
		delete $hModeDepthChr{$chr}{65535};
		foreach my $depth (sort {$hModeDepthChr{$chr}{$b}<=>$hModeDepthChr{$chr}{$a}} keys %{$hModeDepthChr{$chr}}) {
			print CHRDIS "$chrnum\t$chr\t$hMeanDepthChr{$chr}\t$depth\t$hGCperChr{$chr}\n";
			#print STATISTICS "$chr\t$hEffChr{$chr}\t$hMappableChr{$chr}\t$hMeanDepthChr{$chr}\t$depth\t$hChrCoverage{$chr}\n";
			my $coverage = $hChrCoverage{$chr} / $hEffChr{$chr} * 100;
			print STATISTICS "$chr\t$hEffChr{$chr}\t$hMappableChr{$chr}\t$hMeanDepthChr{$chr}\t$depth\t$coverage\n";
			last;
		}
	}
	close STATISTICS;
	close CHRDIS;

	my $local_time = localtime();
	print "$local_time\tdrawing depth distribution by chromosome\n";
	open SH, ">$distribution_dir/run_drawchrdis.sh" or die "$!";
	print SH "$plot_chr_bin $distribution_dir/$pro_name.chrdistribution.data > $distribution_dir/a2 && gnuplot $distribution_dir/a2";
	close SH;

	`sh $distribution_dir/run_drawchrdis.sh`;
}



sub drawChrDepthDistribution {

	print "draw chrdepthion !\n";

	my @chr_order = ();	#(1..22, "X", "Y");

	my %hChrCoverage = ();	#{chr} => coverage
	## read coverage of each chromosome
	open COVERAGE, " $coverage_dir/total_coverage.info" or die "$!";
	while (my $line = <COVERAGE>) {
        	chomp $line;
        	last if ($line =~ /Overall/i);
        	if ($line =~ /Percentage/i) {
                	my @info = split "\:", $line;
                	#$info[-1] =~ s/\%//;
			#$hChrCoverage{$info[0]} = $info[-1];
			my $covered_length = (split /\//, $info[1])[0];
			$hChrCoverage{$info[0]} = $covered_length;
        	}
	}
	close COVERAGE;


	## count mean depth of each chromosome
	open ALLINFO, "$distribution_dir/all.info" or die "$!";
	my %hMeanDepthChr = ();	#{chr} => mean depth
	my %hMappableChr = ();	#{chr} => mappable bases
	my %hEffChr = ();	#{chr} => effective length
	while (my $line = <ALLINFO>) {
		chomp $line;
		my @info = split /\s+/, $line;
		next if ($info[2] !~ /\d{1,}/);
		$hMeanDepthChr{$info[0]} = $info[3] / $info[2];
		$hMappableChr{$info[0]} = $info[3];
		$hEffChr{$info[0]} = $info[2];
	}
	close ALLINFO;


	## count mode depth of each chromosome
	my %hModeDepthChr = ();	#{chr}{depth} => num;
	my @depth_files = `find $distribution_dir -name "*.depth"`;
	chomp @depth_files;
	foreach my $depth_file (@depth_files) {
		my $basename = basename $depth_file;
		$basename =~ s/\.depth//;
		#$basename =~ s/chromosome//;
		open DEPTH, $depth_file or die "$!";
		while (my $line = <DEPTH>) {
			chomp $line;
			my @info = split /\s+/, $line;
			$hModeDepthChr{$basename}{$info[0]} = $info[1];
		}
		close DEPTH;
	}
	
	## read GC percentage file
	my %hGCperChr = ();	#{chr} => gc percentage
	open GCLIST, $gc_file or die "$!";
	while (my $line = <GCLIST>) {
		chomp $line;
		my @info = split /\s+/, $line;
		$hGCperChr{$info[0]} = $info[4];
	}
	close GCLIST;

	open CHRORDER, $chrorder_file or die "$!";
	while (my $line = <CHRORDER>) {
		chomp $line;
		push @chr_order, $line;
	}
	close CHRORDER;

	open CHRDIS, ">$distribution_dir/$pro_name.chrdistribution.data" or die "$!";
	print CHRDIS "chrnum\tchr\tmean_depth\tmode_depth\tgc_percentage\n";
	open STATISTICS, ">$distribution_dir/$pro_name.statistics_by_chr" or die "$!";
	print STATISTICS "chr\teffective_length\tmappable_base\tmean_depth\tmode_depth\tcoverage(%)\n";
	foreach my $chrnum (@chr_order) {
		#my $chr = "canFam2-chr$chrnum";
		#my $chr = "chr$chrnum";
		my $chr = $chrnum;
	#foreach my $chr (sort keys %hModeDepthChr) {
		delete $hModeDepthChr{$chr}{0};
		delete $hModeDepthChr{$chr}{65535};
		foreach my $depth (sort {$hModeDepthChr{$chr}{$b}<=>$hModeDepthChr{$chr}{$a}} keys %{$hModeDepthChr{$chr}}) {
			print CHRDIS "$chrnum\t$chr\t$hMeanDepthChr{$chr}\t$depth\t$hGCperChr{$chr}\n";
			#print STATISTICS "$chr\t$hEffChr{$chr}\t$hMappableChr{$chr}\t$hMeanDepthChr{$chr}\t$depth\t$hChrCoverage{$chr}\n";
			my $coverage = $hChrCoverage{$chr} / $hEffChr{$chr} * 100;
			print STATISTICS "$chr\t$hEffChr{$chr}\t$hMappableChr{$chr}\t$hMeanDepthChr{$chr}\t$depth\t$coverage\n";
			last;
		}
	}
	close STATISTICS;
	close CHRDIS;

	my $local_time = localtime();
	print "$local_time\tdrawing depth distribution by chromosome\n";
	open SH, ">$distribution_dir/run_drawchrdis.sh" or die "$!";
	print SH "$plot_chr_bin $distribution_dir/$pro_name.chrdistribution.data > $distribution_dir/a2 && gnuplot $distribution_dir/a2";
	close SH;

	`sh $distribution_dir/run_drawchrdis.sh`;
}

sub drawDepthDistribution {
	my $local_time = localtime();
	print "\n$local_time\ndrawing depth distribution\n";
	open SH, ">$distribution_dir/run_drawdis.sh" or die "$!";
	print SH "$plot_bin $distribution_dir/$pro_name.distribution.data $distribution_dir 100 10 10 2 $pro_name  > $distribution_dir/a1 && gnuplot $distribution_dir/a1";
	close SH;

	`sh $distribution_dir/run_drawdis.sh`;
}

__END__
nohup perl /panfs/GAG/junli/raid010/pipeline/soap/draw/coverageAndDepthDistribution.V5.SGC2.pl -il soap_result.list -proDir /panfs/GAG/junli/resequencing/project/plant/Soybean_qdj -proName Soybean -ref /share/raid010/resequencing/resequencing/tmp/pub/Genome/Soybean/soybean.fa -outDir /panfs/GAG/junli/resequencing/project/plant/Soybean_qdj -effective 1 -step 2 -chrSuffix chrnum -chrpattern Gm -mem 3g &
