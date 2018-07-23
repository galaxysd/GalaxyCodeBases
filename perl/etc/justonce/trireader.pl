#!/usr/bin/perl
use strict;
use warnings;
use List::Util 1.26 qw(sum0);
use Data::Dump qw(ddx);

my $fhead = 'snp/bam.head';
my @files = qw[snp/SNP.out snp/tomor1.snp.out snp/tumour2.snp.out];

sub getDat($) {
	my $fh = $_[0];
	while(<$fh>) {
		chomp;
		my @dat = split /\t/;
		next unless $dat[6] eq 'PASS';
		my @dep1 = split /,/,$dat[9];
		my @dep2 = split /,/,$dat[10];
		my $sum1 = sum0(@dep1);
		my $sum2 = sum0(@dep2);
		next unless ($sum1+$sum2) > 20;
		return [@dat[0,1,7,6,9,10],$sum1+$sum2];
		return [@dat[0,1,7]];
	}
	return ["\t",-1,"NN"];
}

my (@ChrIDs,%ChrLen);
open H,'<',$fhead or die $!;
while (<H>) {
	chomp;
	my @dat = split /\t/;
	next if $dat[0] ne '@SQ';
	my $cid = (split /:/,$dat[1])[1];
	my $clen = (split /:/,$dat[2])[1];
	push @ChrIDs,$cid;
	$ChrLen{$cid} = $clen;
}
close H;
ddx \@ChrIDs,\%ChrLen;

my ($fha,$fhb,$fhc);
open $fha,'<',$files[0] or die $!;
open $fhb,'<',$files[1] or die $!;
open $fhc,'<',$files[2] or die $!;
<$fha>;<$fha>;<$fhb>;<$fhb>;<$fhc>;<$fhc>;

my @IDs = qw[A B C];
my %FHs = (
	'A' => [$fha,getDat($fha)],
	'B' => [$fhb,getDat($fhb)],
	'C' => [$fhc,getDat($fhc)],
);
#print $FHs{'A'}->[1][0],"--\n";
while (($FHs{'A'}->[1][0] ne "\t") and ($FHs{'B'}->[1][0] ne "\t") and ($FHs{'C'}->[1][0] ne "\t")) {
	ddx \%FHs;
	my %ChrIDs = map { $_ => $FHs{$_}->[1][0] } @IDs;
	my @aID = sort { $FHs{$b}->[1][1] <=> $FHs{$a}->[1][1] } @IDs;	# desc
	ddx \%ChrIDs,\@aID;
}



my $flag=1;
while($flag) {
	my $retA = getDat($fha);
	my $retB = getDat($fhb);
	my $retC = getDat($fhc);
	ddx ($retA,$retB,$retC);$flag = 0;
	;
}
close $fha;
close $fhb;
close $fhc;
