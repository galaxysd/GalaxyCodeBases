#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <files> >out.txt\n" if @ARGV < 2;
#@ARGV;
warn "Input ",scalar @ARGV," Files: [",join('],[',@ARGV),"].\n";

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
		return 0;
	}
}

sub initfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.zst$/) {
			open( $infile,"-|","zstd -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	chomp(my $t = <$infile>);
	my @tt = split("\t",$t);
	map { s/_read(\d)$/_reads$1/ } @tt;
	my $ret = [$infile,\@tt,undef,-1,[]];
	readnext($ret) or die "[x]File [$filename] is empty. $!\n";
	return $ret;
}

print "# InFiles:",join(',',@ARGV),"\n";
my @FH;
for (@ARGV) {
	my $t = initfile($_);
	push @FH,$t;
}



ddx \@FH;

close $_->[0] for @FH;
