#!/usr/bin/perl -w
use strict;
use warnings;
use lib '/home/huxuesong/perl/lib64/perl5/5.8.5/x86_64-linux-thread-multi/';
use lib '/home/huxuesong/perl/lib64/perl5/site_perl/5.8.5/x86_64-linux-thread-multi/';
use lib '/home/huxuesong/gentoo/home/PERL5LIB/';
use Time::HiRes qw ( gettimeofday tv_interval );
#use Galaxy::ShowHelp;

my ($oldpos,$types,$lastype,$oldchr)=(-1,0,0);
my %typebit=( five_prime_utr => 1, three_prime_utr => 2, CDS => 4,
			  exon => 8, mRNA => 16, nongene => 32 );
my %result=( 32 => 'nonGene', 12 => 'CDS', 16 => 'Intron',
			 9 => '5-UTR', 10 => '3-UTR', 25 => '5-UTR', 26 => '3-UTR',
			 13 => '5-UTR CDS', 14 => '3-UTR CDS', 28 => 'CDS',
			 11 => '5-UTR 3-UTR', 29 => '5-UTR CDS', 30 => '3-UTR CDS', );
my ($chr,$pos,$type,$gene);
while (<>) {
	chomp;
	($chr,$pos,$type,$gene)=split /\t/;
	next if $type eq 'gene';
	if ($pos == $oldpos and $chr eq $oldchr) {	# OR bits
		$types |= $typebit{$lastype};
#print "$types\t$type\t",$typebit{$type},"\n";
	} elsif (defined $oldchr) {	# output last
		$types |= $typebit{$lastype};
		print "$oldchr\t$oldpos\t",$result{$types},"\t$gene\n";
		print STDERR "$types\t$oldchr\t$oldpos\t$gene\n" unless $result{$types};
		$types=0;
	} else {	# first line
		$types=0;	# not needed
	}
	$oldpos=$pos;$oldchr=$chr;$lastype=$type;
}
$types |= $typebit{$lastype};
print STDERR "$types\t$oldchr\t$oldpos\t$gene\n" unless $result{$types};
print "$oldchr\t$oldpos\t",$result{$types},"\t$gene\n";

__END__
3-UTR
5-UTR
CDS
Intron
nonGene

3-UTR
3-UTR CDS
5-UTR
5-UTR 3-UTR
5-UTR CDS
CDS
Intron
nonGene
