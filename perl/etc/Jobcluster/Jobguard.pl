#!/usr/bin/perl -w

use strict;
use FindBin qw/$Script $RealBin/;
use lib "$RealBin/..";
use Jobcluster::Jobguard qw/:DEFAULT $MAX_R_JOBS $DEL_BAD_LOG $CHECK_USER/;
use Getopt::Long;

my ($shell,$memory,$sub_line,$submit_queue,$subsh_prefix,$check_job_time,$check_job_usage,$jump_user_determine);
my ($log_directory,$interval,$reqsub,$max_eqw_time,$qhost_timespan,$queue_limit,@df_array,$save_bad_logs);
my ($mem_check_sign,$mem_cycle,$mem_exceed,$mem_record);
my $Help;

GetOptions(
	"s=s"  => \$shell,
	"m=s"  => \$memory,
	"l=i"  => \$sub_line,
	"q=s"  => \$submit_queue,
	"sp=s" => \$subsh_prefix,
	"ct"   => \$check_job_time,
	"cu"   => \$check_job_usage,
	"ju"   => \$jump_user_determine,

	"d=s"  => \$log_directory,
	"t=i"  => \$interval,
	"r"    => \$reqsub,
	"eq=i" => \$max_eqw_time,
	"qs=i" => \$qhost_timespan,
	"ql"   => \$queue_limit,
	"df=s" => \@df_array,
	"sl"   => \$save_bad_logs,

	"ms"   => \$mem_check_sign,
	"mc=i" => \$mem_cycle,
	"me=s" => \$mem_exceed,
	"mr"   => \$mem_record,

    "mj:i" => \$MAX_R_JOBS,

	"h"    => \$Help
);

my $usage_info = "
Name

    Jobguard.pl

Description

    works for Jobguard.pm

Version & author

    Version: 1.0, 2011-08-15
    Author:  JiaWenlong, jiawenlong\@genomics.cn

    Welcome any question, bug-report or suggestion, TIA.

Options

    NOTE: for detail, please perldoc $RealBin/Jobguard.pm

    -s    [s]   shell script. <required>
    -m    [s]   memory the shell script will use. <required>
    -l    [i]   line-number of the sub shells script splited from your -s. [1]
    -q    [s]   queue where you want your jobs to run. (strict format: *.q) [NA]
    -sp   [s]   prefix of each sub shells splited from the original shell ['work']
    -ct         sign of using system command 'time' to check the time each subshell used. [disabled]
    -cu         sign of using system commands 'qstat -j' and 'grep' to get the 'usage' info of each subshell. [disabled]
    -ju         sign of jumping userID determination. [disabled]
                If you have used this script before, and be sure that userID got by it will be definitely your RealOne, set '-ju' for sving your 15s time.
    -d    [s]   directory where the guard_job_log will creat. [./]
                The directory doesnot need to exist as it will be created automatically.
    -t    [i]   cycle time of guarding steps, and it is in second(s). [120,[300]]
    -r          sign of reqsub jobs run to error. [disabled]
                It is suggestted that you should set -r if you prepare to run a multi-step-flow, for the existence of error-jobs.
    -sl         sign of saving logs of bad-jobs for user to check the errors. [disabled]
    -mj   [i]   maximum number of run/qw jobs, while others are hqw. [50]
    -eq   [i]   the maximum time of dealing Eqw status. [[60],100].
    -qs   [i]   time-span of qhost check. [1]
    -ql         the sign of job-submit-queue-limitation in the reqsubing error-jobs. (effective with '-q:*.q') [disabled]
    -df   [s]   info about the disk of which free space you want to check. [NA]
                format: disk,min_space,rls_space
                This parameters can be set in several times, so as to you may have several disks to check, but you should set all disks' info validly.
    -ms         sign of checking memory when jobs are running. [disabled]
    -mc   [i]   time-span of checking memory. (effective with -ms) [5]
    -me   [s]   limited memory that the job can exceed its required memory(-m:vf_mem). (effective with -ms) [0.2]
                it can be decimals(fraction{[0,1]} of vf_mem) or specific memory(in uint 'M/m/G/g').
    -mr         sign of recording memory guard log. (effective with -ms) [abled]
    -h          show this help.

Usage

    perl Jobguard.pl -s <shell_for_qsub> -m <memory_for_run_shell> [Options: -h]

";

die $usage_info if ($Help);
die "perl $Script -s <shell_for_qsub> -m <memory_for_run_shell> [Options: -h]\n" unless($shell && $memory);
my (@qsubmulti,@jobguard);
&std_para();
#------ qsub --------
qsubmulti(@qsubmulti);
#------ check ---------
jobguard(@jobguard);
#-------- sub-routines -----------
sub std_para{
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#--- check paras of qsubmulti ---
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#------ check the shell --------
	die "Cannot find the shellpath: $shell\n" unless(-e $shell);
	push @qsubmulti,"-s:$shell";
	#------ check the memory -------
	if($memory =~ /^\d+$/){
		print STDERR "WARNING: Memory has been reset as 1G\nBecause your input $memory has no unit!\n";
		$memory = '1G';
	}
	push @qsubmulti,"-m:$memory";
	#------- check the subline -----
	if($sub_line){
		if($sub_line !~ /^\d+$/){
			print STDERR "WARNING: Subline has been reset as 1\nBecause your input $sub_line is not whole-digital!\n";
			$sub_line = 1;
		}
		push @qsubmulti,"-l:$sub_line";
	}
	#------- check the submit queue -------
	if($submit_queue){
		push @qsubmulti,"-q:$submit_queue";
	}
	#------- check the subsh prefix -------
	if($subsh_prefix){
		push @qsubmulti,"-sp:$subsh_prefix";
	}
	#------- check the sign of check usage and time ---------
	if($check_job_time){
		push @qsubmulti,"-ct:1";
	}
	if($check_job_usage){
		push @qsubmulti,"-cu:1";
	}
	#------- check the sign of jump user determine ----------
	if($jump_user_determine){
		$CHECK_USER = 0;
	}
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#--- check paras of jobguard ---
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#------ check the log_dir --------
	$log_directory ||= './';
	`mkdir -p $log_directory` unless(-d $log_directory);
	push @jobguard,"-d:$log_directory";
	#------ check the interval -------
	if($interval){
		push @jobguard,"-t:$interval";
	}
	#------- check the reqsub sign -----
	if($reqsub){
		push @jobguard,"-r:1";
	}
	#------- check the max eqw time -----
	if($max_eqw_time){
		push @jobguard,"-e:$max_eqw_time";
	}
	#------- check the qhost time span -----
	if($qhost_timespan){
		push @jobguard,"-h:$qhost_timespan";
	}
	#------- check job queue limit sign ---------
	if($queue_limit){
		push @jobguard,"-ql:1";
	}
	#------- disk free space -------
	if(@df_array){
		push @jobguard,"-df:$_" for @df_array;
	}
	#------- check mem_check_sign --------
	if($mem_check_sign){
		push @jobguard,"-ms:1";
	}
	#------- check mem_cycle -------
	if($mem_check_sign && $mem_cycle){
		push @jobguard,"-mc:$mem_cycle";
	}
	#------- check mem_exceed -------
	if($mem_check_sign && $mem_exceed){
		push @jobguard,"-me:$mem_exceed";
	}
	#------- check memory_record_sign -------
	if($mem_record){
		push @jobguard,"-mr:1";
	}
	#------- check the save of bad-jobs' logs ------
	$DEL_BAD_LOG = !$save_bad_logs;
}
