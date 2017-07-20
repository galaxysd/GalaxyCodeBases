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

unless (prompt_yn("Do you want to proceed with downloading and indexing VFS references?")){
	exit();
}

my $foundexecuableBWA = 0;
retry:
my $bwa_path_to_bin = prompt("Please provide path to bin of BWA ");
print STDERR "VFS will try to index the downloaded reference files using $bwa_path_to_bin\n";
if (-X "$bwa_path_to_bin"){
	print STDERR "BWA can be executed by $bwa_path_to_bin\n";
	$foundexecuableBWA = 1;
} else {
	print STDERR "Failed to execute $bwa_path_to_bin\n";
}

if ($foundexecuableBWA !=1){
	goto retry;
}

#Check CPAN modules
my $pass_CPAN_modules_ck = 1;
my @failed_CPAN = ();
my @required_modules = ("List::Util qw( min max )", "Bio::DB::Sam", "Bio::SeqIO", "Bio::SearchIO", "Bio::DB::Sam", "AppConfig", "AppConfig::Getopt", "Cwd", "Exporter", "File::Which","FileHandle", "FindBin", "Pod::Usage", "Statistics::Descriptive", "File::Copy", "PerlIO::gzip"); ##include all CPAN modules here
for my $curr_module (@required_modules){
	if (try_load($curr_module)) {
		print STDERR "$curr_module loaded: OK\n";
	}
	else {
		print STDERR "$curr_module fails loading.\n";
		$pass_CPAN_modules_ck = 0;
		push (@failed_CPAN, $curr_module);
	}
}
if ($pass_CPAN_modules_ck == 1){
	print STDERR "All CPAN modules checked\n";
} else{
	for (@failed_CPAN){
		print STDERR "Please check $_ is installed properly\n";
	}
	exit;
}

#Check /annotation for pre-defined files
print STDERR "Going to check VFS's bundled annotation files\n";
my @VFS_default_ann = ("hg19.gene.info.in.f", "hg19.refseq.exon.bed", "hg19.repeatmasker.anno");
unless (-d "annotation"){
	print STDERR "annotation folder missing\n";
	exit;
}
my $bundled_annotation_pass = 1;
for my $file (@VFS_default_ann){
	unless (-e "annotation/$file"){
		print STDERR "$file missing in annotation folder\n";
		my $bundled_annotation_pass = 0;
	} else {
		print STDERR "$file checked\n";
	}
}

if ($bundled_annotation_pass == 0){
	print STDERR "Problem(s) of annotation. Try re-download VFS\n";
	exit;
}

#Download nt database. Move to annotation folder
my $NCBI_nt_base = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz";
my $cmd = "wget $NCBI_nt_base";
print STDERR "Fetching nt database\n";
system ($cmd); #fetch database

use Cwd;
my $pwd = cwd();
print "$pwd\n";
opendir(my $dh, $pwd) || die;
	while(readdir $dh) {
	if (/.tar.gz$/i){
		$cmd = "tar -xvzf $_";
		system ($cmd);
	}
}
system ("mv -f nt.* annotation/nt/");
closedir $dh;

#system ("rm -f nt.0*.gz"); #delete the original files

print STDERR "Place the nt database to annotation/nt/\n";

#Download hg19 sequences. Move to references folder
my $hg19 = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz";
my $cmd = "wget $hg19";
print STDERR "Fetching human hg19\n";
system ($cmd);
$cmd = "tar -xvzf chromFa.tar.gz";
system ($cmd);
system ("rm -f chromFa.tar.gz");
#cat all canonical chromosome as one file
system ("cat chrM.fa chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrX.fa chrY.fa >hg19.fa");
system ("mv -f hg19.fa references/human/");
system ("rm chr*.fa");
print STDERR "Place the human hg19 database to references/human/\n";

#Download human decoy sequences. Move to references folder
my $human_decoy = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz";
my $cmd = "wget $human_decoy";
print STDERR "Fetching human decoy\n";
system ($cmd);
$cmd = "gzip -d hs37d5.fa.gz";
system ($cmd);
system ("rm -f hs37d5.fa.gz");
system ("mv -f hs37d5.fa references/human/");

print STDERR "Place the human decoy sequence hs37d5 to references/human/\n";
print STDERR "Indexing databases\n";
system ("$bwa_path_to_bin index -a is references/viral/hbv4.fa");
system ("$bwa_path_to_bin index -a bwtsw references/human/hs37d5.fa");
system ("$bwa_path_to_bin index -a bwtsw references/human/hg19.fa");
#Promopt user to config the config file.
print STDERR "Please now modify the configuration file using the Full system path\n";
print STDERR "Path of nt database RELATIVE to VFS folder: annotation/nt/\n";
print STDERR "Path of viral RELATIVE to VFS folder: references/viral/hbv4.fa\n";
print STDERR "Path of human hg19 RELATIVE to VFS folder: references/human/hs37d5.fa\n";
print STDERR "Path of human decoy RELATIVE to VFS folder: references/human/hg19.fa\n";

sub try_load {
	my $mod = shift;
	eval("use $mod");
	if ($@) {
		return(0);
	}
	else {
		return(1);
 	}
}

sub prompt_yn {
  my ($query) = @_;
  my $answer = prompt("$query (Y/N): ");
  return lc($answer) eq 'y';
}

sub prompt {
  my ($query) = @_;
  local $| = 1;
  print $query;
  chomp(my $answer = <STDIN>);
  return $answer;
}

