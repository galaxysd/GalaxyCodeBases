#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <reference.faidx> <snp files> >out.txt\n" if @ARGV < 3;
#@ARGV;
my $faiN = shift;

my ($cid,%ChrOrder)=(0);
open I,'<',$faiN or die "Error opening $faiN: $!\n";
while(<I>) {
	++$cid;
	my $chr = (split /\t/)[0];
	$ChrOrder{$chr} = $cid;
}
$ChrOrder{'_EOF_'} = 1 + $cid;
close I;
warn "Index:[$faiN], $cid chrosomes found.\nCompare ",scalar @ARGV," Files:[",join('],[',@ARGV),"].\n";

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

my @FH;
my ($id,$flag,@tmpstr) = (0,0);
for (@ARGV) {
	my $t = initfile($_);
	$t->[6] = 1 << $id;
	++$id;
	push @FH,$t;
	$flag += $t->[5];
	push @tmpstr,"$_($t->[6])";
}
print "# InFile(VeenID)s: ",join(',',@tmpstr),"\n# INFO: ",join(',',@SELECTED),"\n";
print join("\t",qw(Chr Pos VeenType),@ARGV),"\n";
while($flag) {
	my @SortedFH = sort { $ChrOrder{$a->[2]} <=> $ChrOrder{$b->[2]} || $a->[3] <=> $b->[3] } @FH;
	#ddx \@SortedFH;
	$flag = 0;
	$flag += $_->[5] for @FH;
	#ddx $flag;
	my $maxSame = 0;
	my $VeenType = $SortedFH[0]->[6];
	for my $i (1 .. $#FH) {
		if ( $SortedFH[0]->[2] eq $SortedFH[$i]->[2] and $SortedFH[0]->[3] == $SortedFH[$i]->[3] ) {
			++$maxSame;
			$VeenType += $SortedFH[$i]->[6];
		} else {
			last;
		}
	}
	my @ChrPos = @{$SortedFH[0]}[2,3];
	my @pDat;
	for my $i (0 .. $maxSame) {
		push @pDat,[$SortedFH[$i]->[4],$SortedFH[$i]->[6]];
		readnext($SortedFH[$i]);
	}
	for my $i ((1+$maxSame) .. $#FH) {
		push @pDat,[[qw(0 0 NA NA 0 0 NA NA)],$SortedFH[$i]->[6]];
	}
	@pDat = sort { $a->[1] <=> $b->[1] } @pDat;
	my @res = map { join(',',@{$_->[0]}) } @pDat;
	print join("\t",@ChrPos,$VeenType,@res),"\n" if $flag;
}


#ddx \@FH;
close $_->[0] for @FH;
