#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump;

#die "Usage: $0 <blood vcf.gz> <sperm vcf.gz>\n" if @ARGV < 2;
die "Usage: $0 <spermA vcf> <spermB vcf> <spermC vcf>\n" if @ARGV < 3;
my ($inA,$inB,$inC)=@ARGV;

my (@FH,%SNP);
for ($inA,$inB,$inC) {
	my $t;
	open $t,'<',$_ or die "$!";
	push @FH,$t;
}

sub getGT($) {
	my $str = $_[0];
	my ($chr,$pos,undef,$ref,$alt,undef,undef,$INFO,undef,$data) = split /\t/,$str;
	my @GeneTypes=split /,/,$alt;
	unshift @GeneTypes,$ref;
	my $GT = (split /:/,$data)[0];
	my @GTs=split /[\/|]/,$GT;
	my @ret = ( $chr, $pos, $GeneTypes[$GTs[0]] . $GeneTypes[$GTs[1]] );
	return \@ret
}

sub readfile($) {	# 1 .. 3
	my $id = $_[0];
	my $fh = $FH[ $id - 1 ];
	while (<$fh>) {
		next if /^#/;
		my ($chr,$pos,undef,$ref,$alt,$QUAL,undef,$INFO,undef,$data) = split /\t/;
		next if $chr =~ /[XYM]/i;
		next if $INFO =~ /INDEL;/;
		next if $QUAL < 20;
		my $GQ = (split /:/,$data)[-1];
		next if $GQ < 20;
		next unless $chr =~ /^chr\d+/;
		my $GT;
		($chr,$pos,$GT) = @{getGT($_)};
		$SNP{$chr}{$pos}{$id} = $GT;
	}
}

readfile($_) for (1 .. 3);
close $_ for @FH;

#ddx \%SNP;
#              141050953 => { 1 => "TT", 2 => "TT" },

my %Stat;
for my $chr (keys %SNP) {
	for my $pos (keys %{$SNP{$chr}}) {
		my %Dat = %{$SNP{$chr}{$pos}};
		my $flag = 0;
		for my $k (keys %Dat) {
			$flag |= 2 ** $k;
		}
		++$Stat{$flag};
#ddx $SNP{$chr}{$pos}; print "$flag\n\n";
	}
}

print "$_\t",$Stat{$_},"\n" for sort {$a<=>$b} keys %Stat;

__END__
perl ./cntsnp3.pl sperm23.5cf.filter sperm24.5cf.filter sperm28.5cf.filter
2	243901
4	268409
6	38623
8	253689
10	33956
12	35961
14	7207

perl ./cntsnp3.pl sperms01.5cf.filter sperms02.5cf.filter sperms03.5cf.filter
2	676931
4	635059
6	60638
8	623031
10	65199
12	58892
14	36960
