#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <snp file> >out.txt\n" if @ARGV < 1;

no warnings 'qw';
my @thePOS = qw(chrom position);
my @SELECTED = qw(ref var eW eC Number_of_watson[A,T,C,G] Number_of_crick[A,T,C,G]);
my %gOrder;
@gOrder{qw(A T C G)} = qw(0 1 2 3);

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
		$in->[5] = $record;
		#ddx \%hash;
		return 1;
	} else {
		@{$in}[2,3] = qw(_EOF_ 0);
		$in->[4] = [qw(0 0 NA NA 0 0 NA NA NA)];
		$in->[5] = '';
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
	print "$t\n";
	my @tt = split("\t",$t);
	map { s/_read(\d)$/_reads$1/ } @tt;
	push @tt,qw(eW eC eA eB);
	my $ret = [$infile,\@tt,undef,-1,[],[]];
	#readnext($ret) or die "[x]File [$filename] is empty. $!\n";
	return $ret;
}

my $FH = initfile($ARGV[0]);
#ddx $FH;

while($FH->[3]) {
	readnext($FH) or last;
	next unless defined $FH->[3];	# skip empty lines
	if (length($FH->[4][1])>1) {
		print $FH->[5],"\n";
		next;
	}
	my $t1 = $gOrder{$FH->[4][0]};
	my $t2 = $gOrder{$FH->[4][1]};
	my @d1 = (split(',',$FH->[4][2]))[$t1,$t2];
	my @d2 = (split(',',$FH->[4][3]))[$t1,$t2];
	my ($d3) = (split(',',$FH->[4][4]))[$t1];
	my ($d4) = (split(',',$FH->[4][5]))[$t1];
	my $s1 = $d1[0] + $d1[1];
	my $s2 = $d2[0] + $d2[1];
	next if ($s1<5 or $s2<5 or ($s1+$s2)<15);
	next if $d1[1]==0 or $d2[1]==0;
	next if $d3 < 5 or $d4 < 5;
	if (($t1+$t2)==3) {
		if ($t1==1 or $t1==2) { # CT
			next if ($d2[1])/($s2) < 0.2;
		} elsif ($t1==0 or $t1==3) { # AG
			next if ($d1[1])/($s1) < 0.2;
		}
	} else {
		next if ($d1[1]+$d2[1])/($s1+$s2) < 0.2;
	}
	print $FH->[5],"\n";
	#ddx $FH;
}
