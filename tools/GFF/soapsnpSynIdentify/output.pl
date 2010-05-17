#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.2.4;

our $opts='i:o:s:bv';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_d);

our $help=<<EOH;
\t-i Input SQLite data file (_result.sqlite)
\t-s Specie name (human)
\t  [Specie name MUST be ths SAME throughout !]
\t-o Output txt file (_output.txt)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='_output.txt' if ! defined $opt_o;
$opt_i='_result.sqlite' if ! $opt_i;
$opt_s='human' if ! $opt_s;

print STDERR "From [$opt_i] to [$opt_o], Specie:[$opt_s]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0
);
our $dbh = DBI->connect('dbi:SQLite:dbname='.$opt_i,'','',\%attr) or die $DBI::errstr;
my $srdhi = $dbh->prepare("SELECT * FROM res$opt_s ORDER BY chrid,position ASC") or warn $dbh->errstr;
$srdhi->execute;
###################
sub write2file($) {
	my $specname=$_[0];
	my (%utr5,%utr3,%scds,%nscds,%unknown,%rest,%intron,%id,%primary,%err,%gene,%exon);
	open FH,'>',$opt_o or die "Error: $!\n";
	print FH "#name\tChrID\tPosition\tPrimaryINF\tRef_base\tSNP_base\tRNA_chg\tAA_chg\tAA_Changed\tStrand\tReg_start\tReg_end\tGroupINF\n";
	while (my $ary_ref = $srdhi->fetchrow_arrayref) {
		my ($chrid,$position,$ref_base,$snp_base,$snp_q,$rna_chg,$aa_chg,$chged,$primary_inf,$start,$end,$strand,$frame,$name,$groups)=@$ary_ref;
		++$id{$chrid};
		++$primary{$primary_inf};
		if ($primary_inf =~ /CDS$/) {
			unless (defined $aa_chg) {
				$rna_chg=$aa_chg=$chged='Error';++$err{$chrid};goto PRINT;
			}
			if ($chged==1) {$chged='Non_Syn';++$nscds{$chrid};goto PRINT;}
			if ($chged==0) {$chged='Syn';++$scds{$chrid};goto PRINT;}
			if ($chged==-1) {$chged='Unknown';++$unknown{$chrid};goto PRINT;}
		} else {$rna_chg=$aa_chg=$chged='N/A' unless defined $aa_chg;}
		if ($primary_inf =~ /(5|f).*UTR$/i) {++$utr5{$chrid};goto PRINT;}
		if ($primary_inf =~ /(3|t).*UTR$/i) {++$utr3{$chrid};goto PRINT;}
		if ($primary_inf =~ /intron/i) {++$intron{$chrid};goto PRINT;}
		if ($primary_inf =~ /exon/i) {++$exon{$chrid};goto PRINT;}
		if ($primary_inf =~ /gene/i) {++$gene{$chrid};goto PRINT;}
		++$rest{$chrid};
PRINT:
		print FH "$name\t$chrid\t$position\t$primary_inf\t$ref_base\t$snp_base\t$rna_chg\t$aa_chg\t$chged\t$strand\t$start\t$end\t$groups\n";
	}
	my $out="\n__END__\n#ChrID\t5'-UTR\t3'-UTR\tSyn_CDS\tNon-syn_CDS\tUnknown_CDS\tIntron\tRest\tError\tGene\tExon\tSum\n";
	print $out;
	print FH $out;
	#my $print_primary=0;
	for (sort keys %id) {
		$utr5{$_}=0 if ! defined $utr5{$_};
		$utr3{$_}=0 if ! defined $utr3{$_};
		$scds{$_}=0 if ! defined $scds{$_};
		$nscds{$_}=0 if ! defined $nscds{$_};
		$unknown{$_}=0 if ! defined $unknown{$_};
		$rest{$_}=0 if ! defined $rest{$_};
		$intron{$_}=0 if ! defined $intron{$_};
		$err{$_}=0 if ! defined $err{$_};
		$gene{$_}=0 if ! defined $gene{$_};
		$exon{$_}=0 if ! defined $exon{$_};
		#$print_primary=1 if $rest{$_}==0;
		my $out="Chr$_\t$utr5{$_}\t$utr3{$_}\t$scds{$_}\t$nscds{$_}\t$unknown{$_}\t$intron{$_}\t$rest{$_}\t$err{$_}\t$gene{$_}\t$exon{$_}\t$id{$_}\n";
		print $out;
		print FH $out;
	}
	#if ($print_primary) {	# Well, always prints
		my $out="\nPrimary_Info Count:\n";
		print $out;
		print FH $out;
		for (sort keys %primary) {
			$out="\t[$_]:\t$primary{$_}\n";
			print $out;
			print FH $out;
		}
	#}
	close FH;
}

write2file($opt_s);

my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";

print STDERR "\033[32;1m Please use [$opt_s] as Specie name in later steps.\033[0;0m\n";
