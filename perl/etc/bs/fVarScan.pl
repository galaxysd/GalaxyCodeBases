#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <VarScan file> >out.txt\n" if @ARGV < 1;
#@ARGV;
my $IN = shift;

my @thePOS = qw(chrom position);
my @SELECTED = qw(ref normal_reads1 normal_reads2 normal_var_freq normal_gt tumor_reads1 tumor_reads2 tumor_var_freq tumor_gt somatic_status);

open I,'<',$IN or die "Error opening $IN: $!\n";
my ($id,$flag,@tmpstr) = (0,0);
my (%hash,@in);
chomp(my $t = <I>);
my @tt = split("\t",$t);
map { s/_read(\d)$/_reads$1/ } @tt;
while(<I>) {
	chomp;
	@hash{@tt} = split /\t/;
	@in[0,1] = @hash{@thePOS};
	$in[2] = [@hash{@SELECTED}];
	next if $in[2]->[9] ne 'Somatic';
	next if $in[2]->[4] ne $in[2]->[0];
	next if ($in[2]->[1]+$in[2]->[2])<20 or ($in[2]->[5]+$in[2]->[6])<15;
	my $t = $in[2]->[7];
	$t =~ s/%$//;
	next if $t < 20;
	#ddx \@in;
	print join("\t",@in[0,1],join(',',@{$in[2]})),"\n";
}
close I;

__END__
1.normal组织没有一个和ref不一样的，并且深度大于20，并且正负链都有支持ref的read。
2. cancer组织深度大于15，突变碱基的频率大于20%。并且正负链都有支持突变的read。
