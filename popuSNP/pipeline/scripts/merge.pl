#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <merge list> <output prefix>\n";
	exit;
}
# $out.se, $out.log
my ($listf,$outf)=@ARGV;
open LST,'<',$listf or die "[x]Error opening $listf: $!\n";
my @files=<LST>;
chomp @files;
close LST;
my $n=$#files;
print '[!]',$n+1," file(s) to be megred.\n";
if ($n == 0) {
	link $files[0],"$outf.sp" or symlink $files[0],"$outf.sp" or exit 1;
} else {
	open O,'>',"$outf.sp" or die "[x]Error opening $outf.sp: $!\n";
	my @FHL;
	for (0..$n) {
		open $FHL[$_][0],'<',$files[$_] or die "[x]Error opening $files[$_]: $!\n";
	}
	#chomp ($FHL[$_][1]=readline $FHL[$_][0]) for (0..$n);
	for (@FHL) {
		$$_[1]=readline $$_[0];	# no need to chomp
		$$_[2]=(split /\t/,$$_[1])[8];
	}
#print '[',$FHL[$_][1],"]\t[",$FHL[$_][2],"]\n" for (0..$n);
	while ($n>-1) {
		@FHL = sort {unless (defined $$a[2]) {return 1} elsif (! defined $$b[2]) {return -1} $$a[2]<=>$$b[2]} @FHL;
		print O $FHL[0][1];
		$FHL[0][1]=readline $FHL[0][0];
		if (defined $FHL[0][1]) {$FHL[0][2]=(split /\t/,$FHL[0][1])[8];}
		 else {--$n; $FHL[0][2]=undef;}
	}
	#print O $FHL[0][1];
#print '[',$FHL[$_][1],"]\t[",$FHL[$_][2],"]\n" for (0..$n);
	close O;
	close $$_[0] for @FHL;
}
$n=@files;
open LOG,'>',"$outf.log";
print LOG "#$n files megred.\ndone !\n";
close LOG;
