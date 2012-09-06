#!/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120330
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <fq_base_list> <genome_list> <sam_path> <out>\n" if @ARGV<1;
my $fqlst=shift;
my $reflst=shift;
my $sampath=shift;
my $outf=shift;

sub openfile($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.xz$/) {
	    open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
    return $infile;
}

open L,'<',$fqlst or die;
my @fqpres;
while (<L>) {
	chomp;
	push @fqpres,$_;
}
close L;

open L,'<',$reflst or die;
my @refs;
while (<L>) {
	chomp;
	push @refs,$_;
}
close L;

my %FH;
for my $fq (@fqpres) {
	for my $ref (@refs) {
		my $fh = openfile("$sampath/${fq}_$ref.sam.gz");
		push @{$FH{$fq}},$fh;
		#$FH{$fq}{$ref}=$fh;
	}
}


__END__
./countmsam.pl fqbase.lst genomes.lst ./check pigv2.stat
