#!/usr/bin/perl -w
###########################################################################
#
#  VirusFinder, a fully automatic pipeline for efficient and accurate detection of viruses,
#               viral mutations, and viral integration sites in host genomes through next
#               generation sequencing data.
#
#  VirusFinder is free software
#
#  Version 2
#  Last update: 02/07/2014
#
#  Contact:       Zhongming Zhao: zhongming.zhao@vanderbilt.edu
#                 Qingguo Wang:   qingguo.wang@vanderbilt.edu
#
#  Organization:  Bioinformatics and Systems Medicine Laboratory
#                 Vanderbilt University Medical Center
#                 Nashville, Tennessee, USA
#
#  sys_check.pl is part of VirusFinder. It checks Java and CPAN modules required by VirusFinder.
#
###########################################################################
# Changelog:
# 03/12/2013   Software release
# 04/16/2013   Code added to check Java version
# 07/19/2013   Code added to check samtools
##############################################################################


print "\nChecking Java version...\n\n";

my $ret = `java -version 2>&1`;
print "$ret\n";

if (index($ret, '1.6') == -1) {
    printf "Warning: The tool Trinity of the Broad Institute may require Java 1.6.\n\n";
}


print "\nChecking SAMtools...\n\n";

$ret = `which samtools 2>&1`;
if (index($ret, 'no samtools') == -1) {
    printf "%-30s\tOK\n\n", 'SAMtools';
}else{
    printf "%-30s\tnot found\n\n", 'SAMtools';
}


my @required_modules = ("Bio::DB::Sam",
                        "Bio::DB::Sam::Constants",
                        "Bio::SeqIO",
                        "Bio::SearchIO",
                        "Carp",
                        "Config::General",
                        "Cwd",
                        "Data::Dumper",
                        "English",
                        "File::Basename",
                        "File::Copy",
                        "File::Path",
                        "File::Spec",
                        "File::Temp",
                        "FindBin",
                        "Getopt::Std",
                        "Getopt::Long",
                        "IO::Handle",
                        "List::MoreUtils",
                        "Pod::Usage",
                        "threads");

print "\nChecking CPAN modules required by VirusFinder...\n\n";
my $count = 0;
for my $module (@required_modules){

	eval("use $module");
	if ($@) {
		printf "%-30s\tFailed\n", $module;
                $count++;
	}
	else {
		printf "%-30s\tOK\n",     $module;
	}
}

if ($count==1){
	print "\n\nOne module may not be installed properly.\n\n";
}elsif ($count > 1){
	print "\n\n$count modules may not be installed properly.\n\n";
}else{
	print "\n\nAll CPAN modules checked!\n\n";
}