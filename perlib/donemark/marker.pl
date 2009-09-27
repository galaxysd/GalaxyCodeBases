#!/bin/env perl
use strict;
use warnings;
use File::Basename;
use lib '/share/raid010/resequencing/user/huxs/release/PERL5LIB/';
use Galaxy::ShowHelp;
use Galaxy::SGE::DoneMark 0.01;

$main::VERSION=0.0.3;

our $opts='i:lv';
our ($opt_i, $opt_l, $opt_v);

our $desc='DoneMark Marker';
our $help=<<EOH;
\t-i file to mark
\t-l the file marking is a file list
\t-v show verbose info to STDOUT
EOH

ShowHelp();

unless ($opt_i) {
	warn "[x]Must specify -i 'file to wait' !\n";
	exit 101;
}
my $mode='FILE';
if ($opt_l) {
	$mode='LIST';
	unless (-s $opt_i) {
		warn "[x]List file [$opt_i] is nothing !\n";
		exit 102;
	}
}

print STDERR "Waiting for [$opt_i] in\033[32;1m $mode \033[0;0mmode ...\n";
#BEGIN
my (@Files,$FileCount,%Stat,%Waiting,$waitname);
if ($mode eq 'FILE') {
	push @Files,$opt_i;
	$FileCount=1;
	$waitname=transname($opt_i);
	$Waiting{$opt_i}=$waitname;
	$Stat{$opt_i}='';
	print STDERR "$FileCount\t[$opt_i]\n" if $opt_v;
} else {
	open L,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";
	while (<L>) {
		chomp;
		push @Files,$_;
		$waitname=transname($_);
		$Waiting{$_}=$waitname;
		++$FileCount;
		$Stat{$_}='';
		print STDERR "$FileCount\t[$_]\n" if $opt_v;
	}
}
$|=1;
BEGIN:
for (keys %Waiting) {
	unless (-e $_) {
		$Stat{$_}='N';
		warn "[x]File Not exists $_: $!\n";
		unlink $Waiting{$_} and $Stat{$_}='n';
		next;
	}
	if (-e $Waiting{$_}) {
		$Stat{$_}='*';
	} else {
		open F,'>',$Waiting{$_} or (warn "[x]Error cresating $_: $!\n" and next);
		print F `date`;	# Will a empty file save lots of space ?
		close F;
		$Stat{$_}='C';
	}
	delete $Waiting{$_};
	--$FileCount;
}
if ($opt_v) {
	print '-' x 75,"\n";
	for (@Files) {
		print "$Stat{$_} $_ $Waiting{$_}\n";
	}
	print "$FileCount not marked.\n";
} else {print '.'}
if ($FileCount > 0) {
	warn "[x]Not all files marked. $FileCount remains.\n";
}
exit $FileCount;

#END
__END__
