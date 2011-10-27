#!/usr/bin/perl
use strict;

#usage: remove the N bases in the sequence, to get a pure sequnence used for simulation


my $ori_len = 0;
my $pure_len = 0;
my $loss_rate = 0;

$/=">";<>;$/="\n";
while (<>) {
	my $title = $_;
	my $seq_name = $1 if($title =~ /^(\S+)/);
	
	$/=">";
	my $seq=<>;
	chomp $seq;
	$/="\n";
	
	$seq =~ s/\s//g;
	$ori_len = length($seq);
	$seq =~ s/[^acgtACGT]//g;
	$seq = uc($seq);
	$pure_len = length($seq);
	
	Display_seq(\$seq);
	print ">$title$seq";
}

$loss_rate = ($ori_len - $pure_len) / $ori_len if($ori_len);

print STDERR "ori_len:$ori_len\tpure_len:$pure_len\tloss_rate:$loss_rate\n";



#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################

