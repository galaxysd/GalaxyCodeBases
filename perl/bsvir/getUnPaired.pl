#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <in_bam> [out_prefix] [fq1] [fq2]\n" if @ARGV < 1;
my ($in,$out,$fq1,$fq2)=@ARGV;
unless (defined $out) {
	$out='grep_'.$in;
	$out =~ s/\.[sb]am(\.gz)?//g;
}
unless (defined $fq2) {
	open( IN,"-|","samtools view -H $in") or die "Error opening $in: $!\n";
	while (my $line = <IN>) {
		if ($line =~ /^\@PG/) {
			$line =~ /CL:"([^"]+)"/ or die "[x]SAM/BAM file header error as [$line]\n";
			my @CLIs = split /\s+/,$1;
			while (my $v = shift @CLIs) {
				if ($v eq '-a') {
					$fq1 = shift @CLIs;
				} elsif ($v eq '-b') {
					$fq2 = shift @CLIs;
				}
			}
		}
	}
	close IN;
}
#die "$in,$out,$fq1,$fq2";

sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.xz$/) {
		open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.bz2$/) {
		open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}
my $FQ1 = openfile($fq1);
my $FQ2 = openfile($fq2);

sub getFQitem($) {
	my $FH = $_[0];
	my @dat;
	for (1 .. 4) {
		my $line = <$FH>;
		die '[x]FQ file not match 1 !' unless defined($line);
		push @dat,$line;
	}
	$dat[0] =~ /^\@([^\/]+)\/(\d+)$/ or die;
	my ($id,$r12) = ($1,$2);
	$id = (split /#/,$id)[0];
	return [$id,$r12,\@dat];
}
sub openpipe($$) {
	my ($cmd,$filename)=@_;
	my $infile;	# undef
	if ( length($filename) ) {
		open( $infile,"|-","$cmd >$filename") or die "Error opening [$cmd],[$filename]: $!\n";
	}
	return $infile;
}
my $OUT1 = openpipe('gzip -9c',"$out.1.fq.gz");
my $OUT2 = openpipe('gzip -9c',"$out.2.fq.gz");

open( IN,"-|","samtools view $in") or die "Error opening $in: $!\n";
while (my $line = <IN>) {
	#my ($id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/,$line;
	#print "$id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize\n";
	my @Dat1 = split /\t/,$line;
	my $line2 = <IN>;
	die '[x]SAM/BAM file not paired !' unless defined($line2);
	my @Dat2 = split /\t/,$line2;
	my $FQ1dat = getFQitem($FQ1);
	my $FQ2dat = getFQitem($FQ2);
	unless ( ($FQ1dat->[0] eq $Dat1[0]) and ($FQ2dat->[0] eq $Dat2[0]) ) {
		die '[x]FQ file not match 2 !',$FQ1dat->[0],'|',$Dat1[0];
	}
	#if ( $ref eq 'chrEBV' or ($flag & 12) ) {
	if ( $Dat1[2] eq 'chrEBV' or ($Dat1[1] & 12) or $Dat2[2] eq 'chrEBV' or ($Dat2[1] & 12) ) {
		if ($Dat1[0] ne $Dat2[0]) {
			die "[x]SAM/BAM file not paired as $Dat1[0] != $Dat2[0] !\n";
		}
		print $OUT1,join('',$FQ1dat->[2]);
		print $OUT2,join('',$FQ2dat->[2]);
	} elsif ( $Dat1[2] eq '*' or $Dat2[2] eq '*' ) {
		warn "-->[$line$line2";
	}
}
close $OUT1;
close $OUT2,
close IN;
close $FQ1;
close $FQ2;
#system("samtools view -bS $out.sam.gz >$out.bam");
