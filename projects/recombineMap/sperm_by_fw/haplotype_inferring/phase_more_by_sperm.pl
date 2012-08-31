#!/usr/bin/perl

=head1 Name

The output sperm SNPs are quality filtered, which is equal or larger than the cutoff score

=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  
  perl phase_more_by_sperm.pl <haplo_file> <snp_file>
  --score <int>  the score cutoff for sperm SNPs
  --depth <int>  the depth cutoff for sperm SNPs
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
my ($Min_score_cutoff, $Min_depth_cutoff, $Verbose,$Help);
GetOptions(
	"score:i"=>\$Min_score_cutoff,
	"depth:i"=>\$Min_depth_cutoff,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Min_score_cutoff ||= 20;
$Min_depth_cutoff ||= 1;
die `pod2text $0` if ($Help);

my $diploid_haplo_file = shift;
my $sperm_snp_file = shift;

my $Chr_name;
my %HAPLO;
my @DiSNP;
my @SpSNP;
my @POS;
my $Sperm_num = 99;
my $Check_num = 5;


open IN,$diploid_haplo_file || die "fail $diploid_haplo_file";
while (<IN>) {
	chomp;
	next if(/^\#/);
	my @t = split /\t/;
	my $chr = $t[0];
	$Chr_name = $chr;
	my $pos = $t[1];
	$HAPLO{$pos}{"hap1"} = $t[2];
	$HAPLO{$pos}{"hap2"} = $t[3];
}
close IN;

print STDERR "finished reading $diploid_haplo_file\n";

#print Dumper \%HAPLO;

open IN, $sperm_snp_file || die "fail $sperm_snp_file";
while (<IN>) {
	next if(/^\#/); 
	my @t = split /\t/;
	my $chr = shift @t;
	my $pos = shift @t;
	my $disnp = shift @t;
	my ($base1, $base2, $disnp_depth, $disnp_score) = ($1,$2,$3,$4) if($disnp =~ /(\w)(\w)\),(\d+),(\d+)$/);
	push @POS, $pos;
	push @DiSNP, "$base1$base2";
	my $num_spsnp = 0;
	for (my $i=0; $i<@t; $i++) {
		my ($base, $depth, $score) = ($1,$2,$3) if($t[$i] =~ /(\w+),(\d+),(\d+)$/);
		if ($score >= $Min_score_cutoff && $depth >= $Min_depth_cutoff) {
			$num_spsnp ++;
		}else{
			$base = "N";
		}
		push @{$SpSNP[$i]}, $base;
	}
	
}
close IN;

print STDERR "finished reading $sperm_snp_file\n";



my $head = "#Chr\tPos\tHap1\tHap2\tSpNum";
for (my $i=1; $i<=99; $i++) {
	my $spermId = ($i<=9) ? "S0$i" : "S$i";
	$head .= "\t$spermId";
}
print $head."\n";


##$j: snp order id; $i: sperm order id:
for (my $j=0; $j<@POS; $j++) {
	
	my $pos = $POS[$j];
	
	my $snp_str;
	my $snp_num;
	for (my $i=0; $i<$Sperm_num; $i++) {
		$snp_str .= "\t$SpSNP[$i][$j]";
		$snp_num ++ if ($SpSNP[$i][$j] ne "N");
	}
	
	if (exists $HAPLO{$pos}) { ##已经phased
		my $hap_base1 = $HAPLO{$pos}{"hap1"};
		my $hap_base2 = $HAPLO{$pos}{"hap2"};
		print "$Chr_name\t$pos\t$hap_base1\t$hap_base2\t$snp_num$snp_str\n";
		next;
	}
	
	my %BELONG;
	my ($dibase1, $dibase2) = ($1,$2) if($DiSNP[$j] =~ /(\w)(\w)/);
	
	$BELONG{$dibase1}{"hap1"} = 0;
	$BELONG{$dibase2}{"hap2"} = 0;
	$BELONG{$dibase1}{"hap2"} = 0;
	$BELONG{$dibase2}{"hap1"} = 0;

	#print STDERR "dibases: $dibase1, $dibase2\n";
	my $jj;
	my $loop;
	
	#print STDERR "pos: $j ".$pos."\t$dibase1$dibase2\n";

	##查每一个sperm, 根据每一个sperm判断属于哪个parent haplotype
	for (my $i=0; $i<$Sperm_num; $i++) {
		
		next if($SpSNP[$i][$j] eq "N"); ##忽略是N的sperm, 没有鉴定到sperm base
		
		my $center_sperm_base = $SpSNP[$i][$j];
		my $same_to_hap1 = 0;
		my $same_to_hap2 = 0;
		
		#print STDERR "sperm id: $i  $SpSNP[$i][$j]  $center_sperm_base\n";

		##往前找$Check_num个
		$jj = $j;
		$loop = 0;
		while ($jj > 0 && $loop < $Check_num) {
			$jj --;
			next if($SpSNP[$i][$jj] eq "N" || !exists $HAPLO{$POS[$jj]});
			my $this_sperm_base = $SpSNP[$i][$jj];
			#$same_to_hap1 ++ if ($this_sperm_base eq $HAPLO{$POS[$jj]}{"hap1"});
			#$same_to_hap2 ++ if ($this_sperm_base eq $HAPLO{$POS[$jj]}{"hap2"});
			if ($this_sperm_base eq $HAPLO{$POS[$jj]}{"hap1"}){
				$same_to_hap1 ++;
			}elsif($this_sperm_base eq $HAPLO{$POS[$jj]}{"hap2"}){
				$same_to_hap2 ++;
			}
			
			$loop++;
		}

		##往后找$Check_num个
		$jj = $j;
		$loop = 0;
		while ($jj < @POS - 1 && $loop < $Check_num) {
			$jj ++;
			next if($SpSNP[$i][$jj] eq "N" || !exists $HAPLO{$POS[$jj]});
			my $this_sperm_base = $SpSNP[$i][$jj];
			#$same_to_hap1 ++ if ($this_sperm_base eq $HAPLO{$POS[$jj]}{"hap1"});
			#$same_to_hap2 ++ if ($this_sperm_base eq $HAPLO{$POS[$jj]}{"hap2"});
			if ($this_sperm_base eq $HAPLO{$POS[$jj]}{"hap1"}){
				$same_to_hap1 ++;
			}elsif($this_sperm_base eq $HAPLO{$POS[$jj]}{"hap2"}){
				$same_to_hap2 ++;
			}
			
			$loop++;
		}
		#print STDERR "#$i\t$same_to_hap1\t$same_to_hap2\n";
		$BELONG{$center_sperm_base}{"hap1"} ++ if ($same_to_hap1 > $same_to_hap2 * 3 );
		$BELONG{$center_sperm_base}{"hap2"} ++ if ($same_to_hap2 > $same_to_hap1 * 3 );
	}
	
	#print STDERR "# ".$BELONG{$dibase1}{"hap1"}."\t".$BELONG{$dibase2}{"hap2"}."\t".$BELONG{$dibase1}{"hap2"}."\t".$BELONG{$dibase2}{"hap1"}."\n";
	my $mode1_count = $BELONG{$dibase1}{"hap1"} + $BELONG{$dibase2}{"hap2"};
	my $mode2_count = $BELONG{$dibase1}{"hap2"} + $BELONG{$dibase2}{"hap1"};
	if ( $mode1_count >=2 && $mode1_count > $mode2_count*3 ) {
		print "$Chr_name\t$pos\t$dibase1\t$dibase2\t$snp_num$snp_str\n";
	}
	elsif ( $mode2_count >=2 && $mode2_count > $mode1_count*3 ) {
		print "$Chr_name\t$pos\t$dibase2\t$dibase1\t$snp_num$snp_str\n";
	}
}




