#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <reference.faidx> <snp files> >out.txt\n" if @ARGV < 3;
#@ARGV;
my $faiN = shift;

my ($id,%ChrOrder)=(0);
open I,'<',$faiN or die "Error opening $faiN: $!\n";
while(<I>) {
	++$id;
	my $chr = (split /\t/)[0];
	$ChrOrder{$chr} = $id;
}
$ChrOrder{'_EOF_'} = 1 + $id;
close I;
warn "Index:[$faiN], $id chrosomes found.\nCompare ",scalar @ARGV," Files:[",join('],[',@ARGV),"].\n";

my @thePOS = qw(chrom position);
my @SELECTED = qw(normal_reads1 normal_reads2 normal_var_freq normal_gt tumor_reads1 tumor_reads2 tumor_var_freq tumor_gt);

sub readnext($) {
	my $in = $_[0];
	if ( ! eof($in->[0]) ) {
		my $record = readline($in->[0]);
		die "readline failed: $!" unless defined $record;
		chomp($record);
		my %hash;
		@hash{@{$in->[1]}} = split /\t/, $record;
		@{$in}[2,3] = @hash{@thePOS};
		$in->[4] = [@hash{@SELECTED}];
		#ddx \%hash;
		return 1;
	} else {
		@{$in}[2,3] = qw(_EOF_ 0);
		$in->[4] = [qw(0 0 NA NA 0 0 NA NA)];
		$in->[5] = 0;
		return 0;
	}
}

sub initfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.zst$/) {
			open( $infile,"-|","zstd -qdc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	chomp(my $t = <$infile>);
	my @tt = split("\t",$t);
	map { s/_read(\d)$/_reads$1/ } @tt;
	my $ret = [$infile,\@tt,undef,-1,[],1,0];
	readnext($ret) or die "[x]File [$filename] is empty. $!\n";
	return $ret;
}

print "# InFiles:",join(',',@ARGV),"\n";
my @FH;
my ($id,$flag) = (0,0);
for (@ARGV) {
	my $t = initfile($_);
	$t->[6] = 1 << $id;
	++$id;
	push @FH,$t;
	$flag += $t->[5];
}

while($flag) {
	my @SortedFH = sort { $ChrOrder{$a->[2]} <=> $ChrOrder{$b->[2]} || $a->[3] <=> $b->[3] } @FH;
	readnext($SortedFH[0]);
	#ddx \@SortedFH;
	$flag = 0;
	$flag += $_->[5] for @FH;
	#ddx $flag;
	my $maxSame = 0;
	for my $i (1 .. $#FH) {
		if ( $SortedFH[0]->[2] eq $SortedFH[$i]->[2] and $SortedFH[0]->[3] == $SortedFH[$i]->[3] ) {
			++$maxSame;
		} else {
			last;
		}
	}
	my @pDat;
	for my $i (0 .. $maxSame) {
		;
	}
}


#ddx \@FH;
close $_->[0] for @FH;
