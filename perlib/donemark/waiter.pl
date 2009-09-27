#!/bin/env perl
use strict;
use warnings;
use File::Basename;
use lib '/share/raid010/resequencing/user/huxs/release/PERL5LIB/';
use Galaxy::ShowHelp;
use Galaxy::SGE::DoneMark 0.01;

$main::VERSION=0.0.2;
my $FallBackSleep=300;	# 300 for Normal, 2 for test.
my $CycleSleep=60;	# 60 for Normal, 1 for test.

our $opts='i:lov';
our ($opt_i, $opt_l, $opt_o, $opt_v);

our $desc='DoneMark Waiter';
our $help=<<EOH;
\t-i file to wait
\t-l the file waiting is a file list
\t-o chech Once and return
\t-v show verbose info to STDOUT
EOH

ShowHelp();

unless ($opt_i) {
	warn "[x]Must specify -i 'file to wait' ! Sleep for ${FallBackSleep}s ...\n";
	sleep $FallBackSleep;
	exit 101;
}
my $mode='FILE';
if ($opt_l) {
	$mode='LIST';
	unless (-s $opt_i) {
		warn "[x]List file [$opt_i] is nothing ! Sleep for ${FallBackSleep}s ...\n";
		sleep $FallBackSleep;
		exit 102;
	}
}

print STDERR "Waiting for [$opt_i] in\033[32;1m $mode \033[0;0mmode ...\n";
print STDERR "Check only Once and return unconfirmed file count.\n" if $opt_o;
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
#for (keys %Waiting) {print "$_\t$Waiting{$_}\n"}
$|=1;
BEGIN:
for (keys %Waiting) {
	if (-e $Waiting{$_}) {
		$Stat{$_}='*';
		delete $Waiting{$_};
		--$FileCount;
	}
}
if ($opt_v) {
	print '-' x 75,"\n";
	for (@Files) {
		print "$Stat{$_} $_\n";
	}
	print "$FileCount not confirmed.\n";
} else {print '.'}
if ($opt_o) {
	print "\n";
	exit $FileCount;
}
if ($FileCount > 0) {
	sleep $CycleSleep;
	goto BEGIN;
} else {print "\n";exit 0;}

#END
__END__
