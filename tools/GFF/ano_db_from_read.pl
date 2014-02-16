#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Getopt::Long;
#use GTF;
use Galaxy::IO;
use Data::Dump qw(ddx);

my ($strCDS, $strmRNA) = qw(CDS mRNA);
GetOptions ('CDS=s' => \$strCDS, 'mRNA=s' => \$strmRNA);

die "Usage: $0 [--cds=CDS] [--mrna=mRNA] <gene_file> <db_file>\n" if @ARGV<2;
my $infs = shift;
my $outfs = shift;

warn "From [$infs] to [$outfs] with [$strCDS] in [$strmRNA]\n";

open I,'<',$infs or die;
open O,'>',$outfs or die;
print O "# From [$infs] to [$outfs]\n";

my $secname = ']';
my (%GeneDat,%Gene2Chr,$strand,$seqname, $primary, $start, $end, $frame, $LN, $groups);
while (<I>) {
	next if /^(#|((\s)*$))/;
	#next if /^#/;
	#print $_;
	if (/^\[([^]]*)\] ([+-])$/) {
#ddx $Gene2Chr{$secname};
		$secname = $1;
		$strand = $2;
		#print "[$secname $strand]\n";
		die if length $secname == 0;
	} else {
		chomp;
		($seqname, $primary, $start, $end, $frame, $LN, $groups) = split /\t/,$_;
		if ( $primary eq $strCDS or $primary eq $strmRNA ) {
			push @{$GeneDat{"$secname\n$strand"}{$primary}},[$start, $end, $frame];
			++$Gene2Chr{$secname}{$seqname};
		}
	}
}
close I;

my (%Err,%Total);
for (sort keys %GeneDat) {
	my ($secname,$strand) = split /\n/,$_;
	my @ref = keys %{$Gene2Chr{$secname}};
	if (@ref == 1) {
		#print O "\n[$secname] $strand $ref[0]\n";
		my $flag = 'PhaseOK';
		my $CDSarray = $GeneDat{$_}{$strCDS};
		my $mRNAarray = $GeneDat{$_}{$strmRNA};
		my @CDS;
		if ($strand eq '+') {
			@CDS = sort { $a->[0] <=> $b->[0] } @{$CDSarray};
		} elsif ($strand eq '-') {
			@CDS = sort { $b->[0] <=> $a->[0] } @{$CDSarray};
		}
		my $CDSLen = $CDS[0][1] - $CDS[0][0] + 1;
		$CDS[0][3] = $CDS[0][2];
		if ( $CDS[0][2] != 0) {
			$CDS[0][2] = 0;
			++$Err{'CDS0'};
			$flag = 'PhaseERR0';
		}
		++$Total{'CDS0'};
		for (1 .. $#CDS) {
			my $phase = -$CDSLen % 3;
			$CDS[$_][3] = $CDS[$_][2];
			if ( $CDS[$_][2] != $phase) {
				$CDS[$_][2] = $phase;
				++$Err{$strand};
				$flag = 'PhaseERR1';
				++$Err{'Reverse'} if $CDS[$_][2] + $CDS[$_][3] == 3;
			}
			++$Total{$strand};
			$CDSLen += $CDS[$_][1] - $CDS[$_][0] + 1;
		}
		$$mRNAarray[0][3] = $CDSLen;
		#$CDSarray = \@CDS;
#print "$secname,$strand\n";
#ddx \@CDS;
#ddx $GeneDat{$_};
		if ($CDSLen % 3) {
			$flag .= ' InCompCDS';
			++$Err{'InCompCDS'};
		}
		print O "\n[$secname] $strand $ref[0] $flag\n";
		print O join("\t",'mRNA',@{$$mRNAarray[0]}),"\n";
		for (@CDS) {
			print O join("\t",'CDS',@$_),"\n";
		}
	} else {
		die "[x] MultiRef(",scalar @ref,") in $secname($strand): [@ref].\n";
	}
}
close O;

#ddx \%Gene2Chr;
ddx \%Err;
ddx \%Total;

__END__
perl ano_db_from_read.pl Dasnov3.gene1 Dasnov3.db
