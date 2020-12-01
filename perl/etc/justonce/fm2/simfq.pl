#!/usr/bin/env perl

use strict;
use warnings;
use Text::CSV;

use Data::Dump qw(ddx);

my $defaultExtends = 100;
my %pChrID = (
	'Mt' => 'NC_012920.1'
);

my $csv = Text::CSV->new ({ binary => 1, auto_diag => 1});
open my $fh, "<", "tSTR.txt";
binmode $fh; # for Windows
$csv->header ($fh, { sep_set => [ ";", ",", "|", "\t" ], detect_bom => 1, munge_column_names => "lc" });
while (my $row = $csv->getline_hr ($fh)) {
    ddx $row;
	$csv->is_missing (0) and next;
	my $ChrID = $row->{'chromosome'};
	if (exists $pChrID{$ChrID}) {
		$ChrID = $pChrID{$ChrID};
	} else {
		$ChrID = 'chr'.$ChrID;
	}
	my ($PCRst,$PCRed) =  ($row->{'pcr_start'},$row->{'pcr_end'});
	unless ($PCRst > 0) {
		$PCRst = $row->{'repeat_start'} - $defaultExtends;
		$PCRst = 1 if $PCRst < 1;
	}
	unless ($PCRed > 0) {
		$PCRed = $row->{'repeat_end'} + $defaultExtends;
	}
	print join("\t",$ChrID,$PCRst,$PCRed,@$row{qw(name motif motif_length)}),"\n";
}

