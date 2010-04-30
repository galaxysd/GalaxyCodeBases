#! /usr/bin/perl -w

use strict;
use Getopt::Long;

my %opts = ();
GetOptions(\%opts, "snpDir:s", "indelFile:s", "minDepth:i", "maxDepth:i", "minQ:i", "sex:i", "outDir:s");

unless (defined $opts{snpDir} && $opts{indelFile} && $opts{outDir}){
	print "perl $0 is used to filter the indel result.\n";
	print "\t-snpDir\tthe snp directory.\n";
	print "\t-indelFile\tthe indel file.\n";
	print "\t-outDir\tthe output directory.\n";
	print "\t-minDepth\tthe minimal depth of the indel[default is 3].\n";
	print "\t-maxDepth\tthe maximal depth of the indel[default is 20].\n";
	print "\t-minQ\t\tthe minimal quality of the indel[default is 20].\n";
	print "\t-sex\t\tthe sex of sample 0 for female 1 for male[default is 0].\n ";
	print "Author:\tChenJiehu\n";
	print "Date:\t2009/08/13\n";
	exit 0;
}
my $snpDir = $opts{snpDir};
my $indelFile = $opts{indelFile};
my $outDir = $opts{outDir};
my $minDepth = 3;
   $minDepth = $opts{minDepth} if (defined $opts{minDepth});
my $maxDepth = 20;
   $maxDepth = $opts{maxDepth} if (defined $opts{maxDepth});
my $minQ = 20;
   $minQ = $opts{minQ} if (defined $opts{minQ});
my $sex = 0;
   $sex = $opts{sex} if (defined $opts{sex});
my $dis = 4;
my %snp = ();

#print "$snpDir\n$indelFile\n$minDepth\n$maxDepth\n$minQ\n$sex\n";

print "get snp starting .........\n";
get_snp();
print "get snp finished .........\n";
print "filter starting ..........\n";
filter_indel();
print "filter finished ..........\n";

sub get_snp{
	my @snpFile = `find $snpDir -name '*.add_ref'`;
	if (@snpFile == 0){
		@snpFile = `find $snpDir -name '*.Q18'`;
	}
	foreach my $snp_file (@snpFile){
		chomp $snp_file;
		print "$snp_file\n";
		open SNP, $snp_file or die "Can't read SNP file. $!";
		while (my $line = <SNP>){
			chomp $line;
			my ($chr, $pos) = (split /\s+/, $line)[0,1];
			foreach ($pos-$dis .. $pos+$dis){
				$snp{$chr}{$_} = 1;
			}
		}
		close SNP;
	}
}
#0		 1		 2		 3		 4	     5		 6		 7		 8
#chromosome      position        type&size       indel-string    fr(+,-,*)   hete/homo   average-quality indel-depth     total-depth
#chr1            4834            D1              T               -           homo        27              5               6		
sub filter_indel{
	my $snp_indel = 0;
	my $homo = 0;
	my $hete = 0;
	my $insertion = 0;
	my $deletion = 0;
	my $total = 0;
	my $sum = 0;
	my $indel_name = (split /\//, $indelFile)[-1];
	open INDEL, $indelFile or die "$!";
	open OUT, ">$outDir/$indel_name.filter" or die "Can't out put the filter file: $!";
	open OUTC, ">$outDir/indel_summary.txt" or die "Can't out put indel summary. $!";
	while (my $line = <INDEL>){
		chomp $line;
		$sum ++;
		my @indel = split /\s+/, $line;
		next if ($indel[6] < $minQ);
		next if (($indel[7] < $minDepth) || ($indel[7] > $maxDepth));
		if (exists $snp{$indel[0]}{$indel[1]}){
			$snp_indel ++;
			next;
		}
		$total ++;
		if ($indel[5] eq "homo"){
			$homo ++;
		}else{
			$hete ++;
		}
		if ($indel[2] =~ /^I/){
			$insertion ++;
		}else{
			$deletion ++;
		}
		#print "$total\t";
		print OUT "$line\n";
	}
	print OUTC "total\thomo\thomo%\thete\thete%\tinsertion\tins%\tdeletion\tdel%\tINDEL_SNP\tINDEL_SNP%\n\n";
	printf OUTC ("%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\n", $total, $homo, $homo/$total*100, $hete, $hete/$total*100, $insertion, $insertion/$total*100, $deletion, $deletion/$total*100, $snp_indel, $snp_indel/$sum*100);
	close OUTC;
	close OUT;
	close INDEL;
}
