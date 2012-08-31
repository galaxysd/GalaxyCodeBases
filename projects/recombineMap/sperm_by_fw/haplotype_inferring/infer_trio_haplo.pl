#!/usr/bin/perl

=head1 Name



=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  
  perl infer_trio_haplo.pl <child_snp_file> <father_snp_file> <mother_snp_file>
  --score <int>  the score cutoff for parent SNPs
  --depth <int>  the depth cutoff for parent SNPs
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple



=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my ($Parent_score_cutoff, $Min_depth_cutoff, $Verbose,$Help);
GetOptions(
	"score:i"=>\$Parent_score_cutoff,
	"depth:i"=>\$Min_depth_cutoff,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Parent_score_cutoff ||= 25;
$Min_depth_cutoff ||= 5;
die `pod2text $0` if ($Help);

my $child_snp_file = shift;
my $father_cns_file = shift;
my $mother_cns_file = shift;

my %Child;
my %Father;
my %Mother;
my $Chr_name;

my @Bases = ("A","C","G","T");

#chr1    10583   G       R       99      G       34      10      16      A       30      11      11      28      0.676126        1.28571 0
#in child, we only keep the heterozygous SNPs
open IN,$child_snp_file || die "fail $child_snp_file";
while (<IN>) {
	my @t = split /\t/;
	$Chr_name = $t[0];
	my $pos = $t[1];
	my $genotype = "$t[5]$t[9]";
	#print "$pos\t$genotype\n";
	$Child{$pos} = $genotype;
}
close IN;

print STDERR "finished reading child snp\n";

open IN, "gzip -dc $father_cns_file|" || die "fail $father_cns_file";
while (<IN>) {
	my @t = split /\t/;
	my $pos = $t[1];
	next if(!exists $Child{$pos} || $t[4] < $Parent_score_cutoff || $t[13] < $Min_depth_cutoff);
	
	my $cns = $t[3];
	my $genotype;
	if ($cns =~ /[ACGT]/) {
		$genotype = "$t[5]$t[5]";
	}else{
		$genotype = "$t[5]$t[9]";
	}
	$Father{$pos} = $genotype;
}
close IN;

print STDERR "finished reading father snp\n";


open IN, "gzip -dc $mother_cns_file|" || die "fail $mother_cns_file";
while (<IN>) {
	my @t = split /\t/;
	my $pos = $t[1];
	next if(!exists $Child{$pos} || $t[4] < $Parent_score_cutoff || $t[13] < $Min_depth_cutoff);
	
	my $cns = $t[3];
	my $genotype;
	if ($cns =~ /[ACGT]/) {
		$genotype = "$t[5]$t[5]";
	}else{
		$genotype = "$t[5]$t[9]";
	}
	$Mother{$pos} = $genotype;
}
close IN;

print STDERR "finished reading mother snp\n";


##analyze belong to father or mother
my $child_count = 0;
my $equal_count = 0;
my $error_count = 0;
my $correct_count = 0;
my $covered_count = 0;
my $same_count = 0;
my $mis_output;
my $head = "#Chr\tPos\tFather's\tMother's\n";
print $head;
foreach my $pos (sort {$a<=>$b} keys %Child) {
	$child_count ++;
	my $child_gt = $Child{$pos};
	my ($child_base1,$child_base2) = ($1,$2) if($child_gt =~ /(\w)(\w)/);
	my $line = "$Chr_name\t$pos\t";
	my %Parent1GT; ##parent共有8种可能的genotype
	my %Parent2GT;
	foreach my $base (@Bases) {
		$Parent1GT{"$child_base1$base"} = 1;
		$Parent1GT{"$base$child_base1"} = 1;
		$Parent2GT{"$child_base2$base"} = 1;
		$Parent2GT{"$base$child_base2"} = 1;
	}
	
	if (exists $Father{$pos} && exists $Mother{$pos}) {
		$covered_count ++;
	}

	my $father_gt = (exists $Father{$pos}) ? $Father{$pos} : "NN";
	my $mother_gt = (exists $Mother{$pos}) ? $Mother{$pos} : "NN";
	

	##满足推断条件的为以下两种可能
	my $is_satisfied = 0;
	if ( exists $Parent1GT{$father_gt} && exists $Parent2GT{$mother_gt} ) {
		if(compare_genotype($father_gt,$mother_gt)){
			$same_count ++;
			next;
		}
		$line .= "$child_base1\t$child_base2\n";
		$correct_count ++;
		$is_satisfied = 1;
		#print "Correct:\t$child_gt\t$father_gt\t$mother_gt\n";
	}
	elsif( exists $Parent2GT{$father_gt} && exists $Parent1GT{$mother_gt} ) {
		if(compare_genotype($father_gt,$mother_gt)){
			$same_count ++;
			next;
		}
		$line .= "$child_base2\t$child_base1\n";
		$correct_count ++;
		$is_satisfied = 1;
		#print "Correct:\t$child_gt\t$father_gt\t$mother_gt\n";
	}
	else{
		$error_count ++;
		#print "Error:\t$child_gt\t$father_gt\t$mother_gt\n";
	}
	print $line if($is_satisfied);
}
print STDERR "child_count: $child_count\n";
print STDERR "phased_count: $correct_count\n";
print STDERR "sameTrio_count:   $same_count\n";
print STDERR "Unphased_count:   $error_count\n";
print STDERR "covered_count:  $covered_count\n";


################################################################
sub compare_genotype{
	my ($father_gt, $mother_gt) = @_;
	my ($fa_gt1,$fa_gt2) = ($1,$2)  if($father_gt =~ /(\w)(\w)/);
	my $is_equal = 0;
	if($father_gt eq $mother_gt || "$fa_gt2$fa_gt1" eq $mother_gt){ #当父母的genotype相等时无法推测子代allele归属
		$is_equal = 1;
	}
	return $is_equal;
}