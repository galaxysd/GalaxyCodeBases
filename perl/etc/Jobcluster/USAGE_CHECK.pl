#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw/abs_path/;
use FindBin qw/$Script $RealScript/;

my ($user_id,$script,$maxloop,$loopsleep,$help);

GetOptions(
	"u=s"	=>\$user_id,
	"s=s"	=>\$script,
	"m=i"	=>\$maxloop,
	"t=i"	=>\$loopsleep,
	"h"		=>\$help,
);

die "Usage:
\tperl $Script <-s <shell_script>> [Options]
Options:
\t-u  user_id
\t-s  shell_script at the end of your 'qsub' commend, should be absolutepath
\t-m  the max loop time when fail to check the shell's job [3]
\t-t  the sleep time [s] when loop [10]
\t-h  show this help\n" if($help || !$script);

chomp($user_id = `whoami`) unless($user_id);
$maxloop ||= 3;
$loopsleep ||= 10;
#$script = abs_path($script);

for (my $i=0;$i<$maxloop;$i++){
	my @job_id = `qstat -u $user_id | grep ' r ' | grep '^[0-9][0-9]' | awk '{print \$1}'`;
	sleep $loopsleep;
	foreach (@job_id){
		chomp($_);
		my $status = `qstat -j $_ 2>&1 | grep -E 'job_number|script_file|usage'`;
		if ($status =~ /script_file.+$script/){
			print STDERR "\n".$status."\n";
			exit();
		}
		sleep 5;
	}
	sleep $loopsleep;
}
die "Cannot check the job of $script!\n";
