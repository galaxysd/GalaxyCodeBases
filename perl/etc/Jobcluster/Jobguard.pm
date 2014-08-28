=head1 NAME && VERSION

    Jobcluster::Jobguard - Qsub single or multi line shell to cluster, and guard their jobs.

    Version:  V1.4.2

    This module is based on some functions of the qsub-sge.pl written by Fan Wei and Hu Yujie, Many Thanks to both of them.
    <Declare> The two names are from script "/share/backup/jiawl/useful/perl_script/qsub-sge.pl".

=head1 SYNOPSIS

    use lib '/share/backup/jiawl/useful/perl_script/module';
    use Jobcluster::Jobguard;

    qsubsingle(<"-s:$shell_path">,<"-m:$memory">,["-c:$change_shell_sign"],["-q:$submit_queue"]);
    qsubmulti(<"-s:$shell_path">,<"-m:$memory">,["-l:$split_shell_line_number"],["-q:$submit_queue"],["-sp:$subsh_prefix"],
              ["-ct:[0/1]"],["-cu:[0/1]"]]);
    jobguard(<"-d:$log_directory">,["-t:$interval"],["-r:$reqsub_sign"],["-n:$logname_prefix"],["-e:$max_eqw_time"],["-h:$qhost_timespan"],
             [<"-ms:$mem_check_sign">,["-mc:$mem_cycle"],["-me:$mem_exceed"],["-mr:$mem_record"]],["-ql:$queue_limit"],
             ["-df:$disk,$min_space,$rls_space"]*);

    * means the parameter is multi tyoe, you can set it several times.

=head1 DESCRIPTION

    Routine 'qsubsingle' allows you to qsub single-line shell script to cluster.
    >> "-s:$shell_path" is the shell(single line) you want to submit to run.
    >> "-m:$memory" is the memory you apply for your jobs. Default as 1G when your input lacks of units.
    >> "-c:$change_shell_sign" defaults as 0, which means this routine couldnot change shell script to the right structure(NOTE<1>).
        If you want your shell change to the structure automatically, input para "-c:1".
    >> "-q:$submit_queue" is the job-queue where you want your job to run. (strict format: *.q)

    Routine 'qsubmulti' allows you to qsub multi-line shell script to cluster.
    >> "-s:$shell_path" is same as above. (multi lines)
    >> "-m:$memory" is same as above.
    >> "-l:$split_shell_line_number" is the Line Number of the small shell script splited from your input $shell_path script. Default as 1.
    >> "-q:$submit_queue" is same as above.
    >> "-sp:$subsh_prefix" is prefix of each sub shells split from the original shell "-s:$shell_path". Default as 'work'.
    >> "-ct:[0/1]" is the sign of using system command 'time' to check the time each subshell used. Default as 0, means disable.
    >> "-cu:[0/1]" is the sign of using system commands 'qstat -j' and 'grep' to get the 'usage' info of each subshell. Default as 0, means disable.

    Routine 'jobguard' helps you to guard status of jobs that are submitted by 'qsubsingle' or 'qsubmulti' routine, and do something right.
    >> "-d:$log_directory" is where the guard_job_log will exist, The directory does Not need to exist, for it will be created by this routine.
    >> "-t:$interval" is the cycle time of guarding jobs, and it is in second(s). Default as 300(s).
    >> "-r:$reqsub_sign" is default as '0', which means make the reqsub-error-job function disabled. Set $reqsub_sign as '1' to enable it.
    >> "-n:$logname_prefix" is the prefix of job_guard_log's name.
        when it is 'multi line shell' or one 'single line shell' job to guard, it will be set as same as the name of shell script in despite of your input.
        when it is more than one 'single line shell' job to guard, it will be set as the variable $logname_prefix; or it will be set as 'jobguard' by default.
    >> "-e:$max_eqw_time" is the max time of Eqw status for each job. Default as 60, but is limited below 100.
    >> "-h:$qhost_timespan" is the time of each qhost check, default as 1, which means one qhost-check one cycle.
    >> "-ms:$mem_check_sign" is the memory guard sign, default as 0, which means disable the memory-guard function.
    >> "-mc:$mem_cycle" is similar to "-h:$qhost_timespan", default as 5. (effective with "-ms:1")
    >> "-me:$mem_exceed" is the limited memory that the job can exceed its required memory(vf_mem). (effective with "-ms:1")
        It can be decimals(fraction{[0,1]} of vf_mem) or specific memory(in uint 'M/m/G/g').
    >> "-mr:$mem_record" is the sign of memory-guarding recording, default as 1, which means enable record function. (effective with "-ms:1")
        If this function is enabled, memory-guard log will be found at the end of LOG file.
    >> "-ql:$queue_limit" is the sign of job-submit-queue-limitation in the reqsubing error-jobs.
        Default as 1, which means keep the original job-queue info. Want to resubmit error-jobs to any job-queues, set it as 0.
    >> "-df:$disk,$min_space,$rls_space" is about the disk of which free space you want to check.
        $disk is path of the disk to check via command 'df -h $disk';
        $min_space is the minimum avail space the disk remains, once less, hold all jobs guarded;
        $rls_space is the minimum avail space the disk remains to release all holded jobs.
        Requirement: $rls_space should be at least 50G more than $min_space for safty, else it will die as an exception.
        The available units are 'K,k,M,m,G,g,T,t'.
        This parameters can be set in several times, so as to you may have several disks to check, but you should set all disks' info validly.

=head1 NOTE

    <1> # guard works method #
        1) All shell qsubed should have the '.... && perl -e 'print STDERR "This-work-is-completed\n"' structure to allow the routine 'jobguard' works.
        2) If single-line shell script doesnot have the structure by itself, then set "-c:1" will change it, or routine 'qsubsingle' will die as an exception.
        3) multi-line shell script is not required to have this structure, for routine 'qsubmulti' will add it at the end of each line automatically.
    <2> # job guard method #
        Use `whoami` command to get the userID for the following " qstat -u userID ".
        Q: IF `whoami` could not get your valid name to " qstat -u userID ", how to solve this problem ?
        A: First, Add "use Jobcluster::Jobguard qw/:DEFAULT $USER/;" in your perl script to export the our variable $USER;
           Then, Change $USER to your valid name to let the " qstat -u $USER " works.
    <3> # reqsub sign warning #
        If you set the sign as default ('0'), which means you disable the reqsub-error-job function, once there exists any error job, It will die as an
        exception-Warning when the jobguard routine finishes its all-job guarding. Because the error jobs cannot give you the right results, and your 
        next step which needs those right results will not succeed. The die just gives you chance to reqsub the error jobs by hand. 
        * absolute path of Shell scripts of all error jobs need to reqsub CAN be found at the end of the related guard_job_log.
    <4> # elapsed time of job-guarding #
        Everytime routine 'jobguard' runs, it will record the elapsed time in second(s). It will help you to know how much time each step takes.
        Check it at the end of guard_job_log of each step ("Jobcluster::Jobguard::jobguard Guard time: XXXXXs"), But you know it isn't the cpu time.
    <5> # memory guard method #
        Use `qstat -j JobID | grep usage` to get the usage-info.
        Try to set the "-mc:$mem_cycle" larger to lighten the burden of SGE.
        When memory-guard function enabled, once there exists memory-error jobs, the routine 'jobguard' will die as an exception when finishs its 
        jobs' guarding work. So it is easy to undertand that any flow will die at the step whose scripts' max_memory turn out error.
    <6> # max run/qw job number method #
        Use "qsub ... -hold_jid XXX ..." to achieve the max run/qw job number control.
        A Variable named as $MAX_R_JOBS, which defaults as 50, is the max number of run/qw jobs.
        Add "use Jobcluster::Jobguard qw/:DEFAULT $MAX_R_JOBS/;" in your perl script to export the our variable $MAX_R_JOBS, then change it as you like.
    <7> # delete logs of bad jobs #
        A our type variable named as $DEL_BAD_LOG, defaults as 0, which means enable the deletion.
        Add "use Jobcluster::Jobguard qw/:DEFAULT $DEL_BAD_LOG/;" in your perl script to export this variable, and you can change it.
    <8> # operations when disk are full #
        As your setting "-df:...." in the Routine 'jobguard', this module can do something right when disk reaches the limitation to avoid writing.
        To the running jobs, use command 'qmod -s' to let them sleep, while using 'qmod -us' to wake up them when disk is ok;
        To the qw/hqw jobs, use command 'qhold' to let them hold, while using 'qrls' to unhold them.

=head1 AUTHOR && CONTACT

    Author :  Jiawenlong at 2011/03/06
    Contact:  jiawenlong@genomics.org.cn

    Welcome any question, bug-report or suggestion, TIA.

=head1 UPDATE LOG

    V1.0  
          achieve the basic function.
    V1.1  
          1) add the WARNING-log for nonexistence of jobs' log.
          2) optimize some perldoc infos.
    V1.2  
          1) update the routines' input para model.
          2) add the single-line shell change function.
          3) use 'qmod' to deal eqw-jobs, max deal time defaults as 60 for each job. (Many Thanks to Ye Rui and Chen Shuisheng)
          4) add the qhost time span para (-h) in routine 'jobguard', and default as 1, which means one qhost-check one cycle.
          5) optimize the deadnode-log.
    V1.3
          1) add four memory guard paras in routine 'jobguard'. (-ms,-mc,-me,-mr)
          2) add the -q para to let user specify the queue to run on.
          3) add the userID feedback step.
          4) add the submit-queue-limit para in the routine 'jobguard'.
          5) optimize the memory-record-log for memory-error-shell-script.
          6) add the max_run_jobs function.
          7) add the function that deal the line prefixed by '#' in shell script.
          8) add the function that check the used time and usage info of each sub shell.
          9) add the parameter to set prefix of sub-shell.
    V1.4
          1) add the check function of disk available space.
          2) add the save fuction of bad-jobs' logs for users to check errors manually.
=cut

package Jobcluster::Jobguard;

use strict;
use warnings;
use File::Basename qw/basename dirname/;
use Cwd qw/abs_path/;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, $VERSION, @EXPORT_OK, %EXPORT_TAGS);
@ISA = qw(Exporter);
@EXPORT = qw(jobguard qsubmulti qsubsingle);
@EXPORT_OK = qw(mendnumber $USER $MAX_R_JOBS $DEL_BAD_LOG $CHECK_USER);
%EXPORT_TAGS = ( DEFAULT => [qw(jobguard qsubmulti qsubsingle)],
                 OTHER   => [qw(mendnumber)]);
#----- version --------
$VERSION = "1.4.2";

#----- required variables ------
our ($USER,%JOB_HASH,$CHECK_USER,$MAX_R_JOBS,$DEL_BAD_LOG);
chomp($USER = `whoami`);
$CHECK_USER = 1; ## default need to be check by user
my $SINGLE_CHECK_SIGN = 1;
$MAX_R_JOBS = 50; ## initial the max_run_job_number as 50, others hqw
my $USAGE_CHECK_PERL_SCRIPT;
my (%DEVICE,$DEVICE_HOLD);
$DEVICE_HOLD = 0;
foreach (@INC){
	$USAGE_CHECK_PERL_SCRIPT = abs_path("$_/Jobcluster/USAGE_CHECK.pl") if(-e abs_path("$_/Jobcluster/USAGE_CHECK.pl"));
}
die "Cannot locate the USAGE_CHECK_PERL_SCRIPT\nIt should exist in the same directory with this module\n" unless($USAGE_CHECK_PERL_SCRIPT);
$DEL_BAD_LOG = 1;
#---------- check the status of submitted jobs until all done ----------
#----------------- if error, then rm deadlogs and resub ----------------
sub jobguard{
	my %para_hash;
	&jobguard_para(\%para_hash,@_);
	#-------- necessary para ---------
	my ($log_dir,$interval,$reqsub,$max_eqw_time,$qhost_timespan) = ($para_hash{'d'},$para_hash{'t'},$para_hash{'r'},$para_hash{'e'},$para_hash{'h'});
	#-------- other para ---------
	my ($mem_check_sign,$mem_cycle,$mem_exceed,$mem_record) = ($para_hash{'ms'},$para_hash{'mc'},$para_hash{'me'},$para_hash{'mr'});
	my $queue_limit = $para_hash{'ql'};
	unless(exists($JOB_HASH{START_TIME})){
		chomp(my $time_anchor = `date +\%s`);
		$JOB_HASH{START_TIME} = $time_anchor;
	}
	my $log = ((exists($JOB_HASH{SHELL}) && @{$JOB_HASH{SHELL}} == 1)?(basename @{$JOB_HASH{SHELL}}[0]):($para_hash{'n'} || 'jobguard')).".$JOB_HASH{START_TIME}.log";
	delete $JOB_HASH{SHELL};
	delete $JOB_HASH{START_TIME};
	die "No jobs have been submitted(no jobs to guard), so you cannot not run Jobcluster::Jobguard::jobguard\n" if(scalar(keys %JOB_HASH) == 0);
	open (LOG,">$log_dir/$log") || die "fail $log_dir/$log: $!\n"; ## open log
	print "Please check $log_dir/$log for jobs' stat-info\n"; ## tell user to check log
	my ($finishjob,$reqsub_number,$error_number,$mem_error_number) = (0,0,0,0);
	my $jobnum = scalar(keys %JOB_HASH);
	my @finished_jobs = ();
	my @need_reqsub = ();
	my %Eqw_Time; # this hash doesnot need to clear real-time
	my %MEM_USE;
	chomp(my $start_second = `date +\%s`);
	my $cycle_time = 0;
	while(($finishjob + scalar(@need_reqsub)/2 + $mem_error_number) != $jobnum){
		$cycle_time++;
		my @need_qmod = ();
		sleep $interval;
		####### get the queue status, find the dead nodes #################
		my @queue_status;
		if($cycle_time % $qhost_timespan == 0){
			HOST:{
				@queue_status=split /\n/,`qhost`;
				if(scalar(@queue_status)>0 && ($queue_status[0]!~/HOSTNAME/ || $queue_status[0]!~/LOAD/)){
					sleep 60;
					print LOG "\n<WARNING>:\tThe shell order 'qhost' gets invaid info! Redoing it Now!\n\n";
					redo HOST;
				}
			}
		}
		else{
			print LOG "\nINFO:\tQhost-time-span is $qhost_timespan, this time($cycle_time) skip qhost check...\n\n";
		}
		my %deadnode=map{
			if(/^compute/){
				my @ele=split;
				my $sign='';
				$sign.=$_ for @ele[3..7];
				if($sign=~/-/){
					($ele[0],1);
				}
				else{
					();
				}
			}
			else{
				();
			}
		} @queue_status;
		#------------ get job status ------------
		my @jobstat;
		STAT:{
			@jobstat=split /\n/,`qstat -u $USER`;
			if(scalar(@jobstat)>0 && ($jobstat[0]!~/job-ID/ || $jobstat[0]!~/name/ || $jobstat[0]!~/state/)){
				sleep 60;
				print LOG "\n<WARNING>:\tThe shell order 'qstat -u $USER' gets invaid info! Redoing it Now!\n\n";
				redo STAT;
			}
		}
		my %jobstat=map{
			if(/^\s*\d/){
				s/^\s+//;
				my @ele=split;
				if(/@(compute-\d+-\d+)\D/){
					if(exists($deadnode{$1})){
						($ele[0],"dqueue:$1");
					}
					else{
						@ele[0,4];
					}
				}
				else{
					@ele[0,4];
				}
			}
			else{
				();
			}
		} @jobstat;
		my @jobID=sort {$a<=>$b} keys %JOB_HASH;
		foreach my $jobID(@jobID){
			if(exists $jobstat{$jobID}){ ## jobs are on the cluster now
				my $qdel_sign = 1; ## default as qdel those jobs
				my $memory_ok = 1; ## default as job's used memory is ok, so can reqsub
				if($jobstat{$jobID}=~/^[dE]/i){ ## d or E status
					if($jobstat{$jobID} =~ /^dqueue:(compute.+)$/){ # dead node
						print LOG "\n<WARNING>:\tJob $jobID is running on a dead node($1)! Try to qdel job $jobID...\n\n";
					}
					elsif($jobstat{$jobID} eq 'Eqw'){ # Error qw status
						print LOG "\n<WARNING>:\tJob $jobID is in $jobstat{$jobID} status! ";
						$Eqw_Time{$jobID} = 0 unless(exists($Eqw_Time{$jobID}));
						$Eqw_Time{$jobID}++; # Eqw time + 1
						if($Eqw_Time{$jobID} > $max_eqw_time){ # This JOB Eqw has exceeeded $max_eqw_time times
							print LOG "But Its Eqw time has exceeded $max_eqw_time times!\n\n";
						}
						else{
							print LOG "Try to qmod -cj job $jobID at the end of this guard-cycle...\n\n";
							push @need_qmod,$jobID; # record the jobs need to be qmod
							$qdel_sign = 0; # not qdel this job
						}

					}
					else{ # other error status, like dr
						print LOG "\n<WARNING>:\tJob $jobID is in $jobstat{$jobID} status! Try to qdel job $jobID...\n\n";
					}
				}
				else{ ## jobs not in SGE-error status, next check memory and cpu-time(if ok)
					#-------- running jobs' memory-check -----------
					if($mem_check_sign && ($cycle_time % $mem_cycle == 0) && $jobstat{$jobID} eq 'r'){ ## running jobs' memory-check do
						my $memory_max_usage;
						($memory_max_usage,$memory_ok,$qdel_sign) = &memory_guard($jobID,$mem_exceed);
						if($mem_record && $memory_max_usage ne '-1'){ ## record the memory_max_used
							my $shell = ${$JOB_HASH{$jobID}}[0];
							${$MEM_USE{$shell}}[0] = ${$JOB_HASH{$jobID}}[1]; ## required memory
							${$MEM_USE{$shell}}[1] = $memory_max_usage; ## max used memory
							${$MEM_USE{$shell}}[2] = 'Memory-Error_Need-modify-and-Reqsub' unless($memory_ok); ## memory-error-record-warning
						}
						print LOG "\n<WARNING>:\tJob $jobID exceeds required memory(${$JOB_HASH{$jobID}}[1]) more than ".(($mem_exceed =~ /m/i)?$mem_exceed:"its $mem_exceed")."! Try to qdel job $jobID...\n\n" if($qdel_sign);
					}
					else{ ## everything is fine...
						$qdel_sign = 0; # not qdel this job
					}
				}
					if($qdel_sign){
						&qdel_job($jobID);  ##  `qdel $jobID`;
						if($DEL_BAD_LOG){
							push my @deadlog , (@{$JOB_HASH{$jobID}}[0].'.o'.$jobID) , (@{$JOB_HASH{$jobID}}[0].'.e'.$jobID);
							&check_rm($_) for @deadlog;
							print LOG "INFO:\tlogs of dead job $jobID have been deleted!\n";
						}
						if($reqsub){
							if($memory_ok){
								print LOG "INFO:\tAnd the new job for running the shell: ".${$JOB_HASH{$jobID}}[0]." is ".&qsubsingle("-s:${$JOB_HASH{$jobID}}[0]","-m:${$JOB_HASH{$jobID}}[1]","-r:2","-q:".(($queue_limit)?${$JOB_HASH{$jobID}}[2]:'NA'))."\n\n";
								$reqsub_number++;
							}
							else{
								print LOG "INFO:\tShell: ${$JOB_HASH{$jobID}}[0] exceeds required memory more than ".(($mem_exceed =~ /m/i)?$mem_exceed:"its $mem_exceed").", So it cannot qsub again!\n";
								$mem_error_number++;
							}
							$error_number++;
						}
						else{
							print LOG "\n<WARNING>:\t".@{$JOB_HASH{$jobID}}[0]." is not reqsubed!\n\n";
							if(&test_store(\@need_reqsub,@{$JOB_HASH{$jobID}}[0])){
								push @need_reqsub,@{$JOB_HASH{$jobID}};
								$error_number++;
							}
						}
						delete $JOB_HASH{$jobID};
					}
			}
			else{ ## jobs are not on the cluster now
				my $log=@{$JOB_HASH{$jobID}}[0].".e$jobID";
				my $sign=(-e $log)?`tail -1 $log`:'The log not exists';
				my $key=0;
				if($sign=~/This-work-is-completed/){
					$key=1;
				}
				else{
					sleep 120; ## give the second chance to check log
					$sign=(-e $log)?`tail -1 $log`:'The log not exists';
					$key=1 if($sign=~/This-work-is-completed/);
				}
				if($key==1){
					$finishjob++;
					print LOG "INFO:\tJob $jobID has finished successfully!\n";
					delete $JOB_HASH{$jobID};
					push @finished_jobs,$jobID;
				}
				else{
					if($sign eq 'The log not exists'){
						print LOG "\n<WARNING>:\tCannot find the output log of job $jobID. ($log)\n";
					}
					else{
						print LOG "\n<WARNING>:\toutput log of job $jobID doesnot have the 'This-work-is-completed' ($log)\n";
					}
					print LOG "<WARNING>:\tJob $jobID ruined by some reason! Try to qdel job $jobID...\n\n";
					&qdel_job($jobID);  ##  `qdel $jobID`;
					if($DEL_BAD_LOG){
						push my @deadlog , (@{$JOB_HASH{$jobID}}[0].'.o'.$jobID) , (@{$JOB_HASH{$jobID}}[0].'.e'.$jobID);
						&check_rm($_) for @deadlog;
						print LOG "INFO:\tlogs of dead job $jobID have been deleted!\n";
					}
					if($reqsub){
						print LOG "INFO:\tAnd the new job for running the shell: ".${$JOB_HASH{$jobID}}[0]." is ".&qsubsingle("-s:${$JOB_HASH{$jobID}}[0]","-m:${$JOB_HASH{$jobID}}[1]","-r:2","-q:".(($queue_limit)?${$JOB_HASH{$jobID}}[2]:'NA'))."\n\n";
						$reqsub_number++;
						$error_number++;
					}
					else{
						print LOG "\n<WARNING>:\t".@{$JOB_HASH{$jobID}}[0]." is not reqsubed!\n\n";
						if(&test_store(\@need_reqsub,@{$JOB_HASH{$jobID}}[0])){
							push @need_reqsub,@{$JOB_HASH{$jobID}};
							$error_number++;
						}
					}
					delete $JOB_HASH{$jobID};
				}
			}
		}
		#------- do qmod tehe Eqw jobs -------
		if(@need_qmod != 0){
			print LOG "Now Try to qmod -cj the below jobs:\n";
			print LOG "$_\n" for @need_qmod;
			`qmod -cj @need_qmod`; # clear the Eqw status
			#------- although the log not exists, but try del them for safe -------
			foreach my $eqw_jobID (@need_qmod) {
				push my @deadlog , (@{$JOB_HASH{$eqw_jobID}}[0].'.o'.$eqw_jobID) , (@{$JOB_HASH{$eqw_jobID}}[0].'.e'.$eqw_jobID);
				&check_rm($_) for @deadlog;
			}
			print LOG "qmod finished\n";
		}
		#-------- check the disk free space -------
		if(scalar(keys %DEVICE) != 0){
			if($DEVICE_HOLD == 0){ ## no job qhold for the disk reason
				unless(&disk_free){ ## never do anything else for remaining the user's freedom
					print LOG "\n<WARNING>:\tDisk Avail space is not ok!\n";
					&hold_jobs;
					$DEVICE_HOLD = 1;
					print LOG "\n<WARNING>:\thold all jobs done!\n";
				}
			}
			else{
				if(&disk_free){ ## never do anything else for remaining the user's freedom
					print LOG "Disk Avail space is ok!\n";
					&release_jobs;
					$DEVICE_HOLD = 0;
					print LOG "release all holded jobs done!\n";
				}
			}
		}
		#-------- output the each-cycle-info ---------
		my @running_jobs = sort {$a<=>$b} keys %JOB_HASH;
		print LOG "\nINFO:\tGuard-$cycle_time-time: all:$jobnum finsih:$finishjob [h]run/[h]qw/s:".scalar(@running_jobs)." error:$error_number ".(($reqsub)?"reqsub:$reqsub_number":("need_reqsub:".(scalar(@need_reqsub)/2)))." memory_error:$mem_error_number".' at '.`date`."\n";
		print LOG "Jobs have been holded for the reason of not-enough-free-space-in-disk\n" if($DEVICE_HOLD == 1);
	}
#	undef %JOB_HASH;
	chomp(my $end_second = `date +\%s`);
	my @running_jobs = sort {$a<=>$b} keys %JOB_HASH;
	die "!ERROR!:\tThe Last number of running jobs is not '0'!\n" if(@running_jobs != 0); 
	delete $JOB_HASH{$_} for keys %JOB_HASH;
	if($mem_check_sign && $mem_record){ ## output the memory recorded of each shell_script
		my @shells = sort keys %MEM_USE;
		print LOG "\n---------- Shell scripts' recorded max-memory-usage ------------\n";
		if(@shells == 0){
			print LOG "No job's memory usage is record, because the function didnot run during the whole guarding time!\n\n";
		}
		else{
			print LOG "Shell scripts' recorded memory usage are as followed:\n(<shell> <vf_memory> <max_used_memory> [warning])\n";
			print LOG "$_\t${$MEM_USE{$_}}[0]\t${$MEM_USE{$_}}[1]".(defined(${$MEM_USE{$_}}[2])?"\t${$MEM_USE{$_}}[2]":'')."\n" for @shells;
		}
	}
	print LOG "\nJobcluster::Jobguard::jobguard Guard time: ".($end_second - $start_second)."s\n\n";
	print LOG "All jobs are finished successfully\n" if($jobnum == $finishjob);
	if(!$reqsub && @need_reqsub/2 > 0){
		print LOG "\n";
		print LOG "------------ As 'reqsub function' is closed -------------\n";
		print LOG "----------- And Error job exists IN this time -----------\n";
		print LOG "------------- So Shell need to be reqsubed --------------\n";
		print LOG "vf=$need_reqsub[($_-1)*2+1]\t$need_reqsub[($_-1)*2]\n" for (1..scalar(@need_reqsub)/2);
		print LOG "\n\n";
		print LOG "As 'reqsub function' is closed, so some Jobs Error lead to ALL Flow Shutdown!\nCHECK this log carefully!!!\n";
	}
	my $final_exception = '';
	$final_exception .= "\nAs 'reqsub function' is closed, so some Jobs Error lead to ALL Flow Shutdown!\n" if(!$reqsub && @need_reqsub/2 > 0);
	$final_exception .= "\nAs some shell scripts' max-used-memory have exceeded their required memory more than ".(($mem_exceed =~ /m/i)?$mem_exceed:"their $mem_exceed").", So they have been qdel.\nThis is the end of FLOW. Please optimize the shell script or modify your memory-paras\n" if($mem_error_number != 0);
	if(length($final_exception) != 0){ ## will die?
		print LOG $final_exception;
		die $final_exception; ## yes, this is the exception.
	}
	close LOG;
}
#-------------------- rm the logs of dead jobs --------------------
sub check_rm{
        my $rmfile=shift;
        `rm $rmfile` if(-e $rmfile);
}
#-------------------- qsub the jobs -------------------
sub qsubsingle{
	&check_user if($CHECK_USER);
	my %para_hash;
	&qsubsingle_para(\%para_hash,@_);
	my ($sh,$memory,$record,$change_sign,$qsub_queue) = ($para_hash{'s'},$para_hash{'m'},$para_hash{'r'},$para_hash{'c'},$para_hash{'q'});
	(my $line_sum) = (`wc -l $sh` =~ /^(\d+)\s+\S+$/); ## get the line number of mutil-shell
	die "Shell $sh is not single-line\n" if($line_sum != 1 && $SINGLE_CHECK_SIGN);
	chomp(my $line_info = `cat $sh`);
	if($line_info =~ /^\#/){
		print STDERR "single line shell but prefixed by '#'!\n";
		exit();
	}
#	if($line_info !~ /&&/ || $line_info !~ /This-work-is-completed/ || $line_info !~ /STDERR/){
	if($line_info !~ /\s+&&\s+\(*perl\s+-e\s+'print\s+STDERR\s+".*This-work-is-completed\\n"'\)*\s*$/){
		die "$sh should ends with '.... && (perl -e 'print STDERR \"This-work-is-completed\\n\"')\n" unless($change_sign);
		#------- change the original shell ---------
		open(SH,">$sh") || die "fail change $sh: $!\n";
		print SH "($line_info) && perl -e 'print STDERR \"This-work-is-completed\\n\"'\n";
		close SH;
	}
    my $bin=dirname $sh;
	my $queue_para = ($qsub_queue eq 'NA')?'':"-q $qsub_queue";
    my $qsub = "qsub -o $bin -e $bin $queue_para -l vf=$memory ";
	## check the number of qsub jobs, max_run_jobs number works
	my @sort_r_jobs = grep {/^\d+$/} (keys %JOB_HASH);
	@sort_r_jobs = sort {$a<=>$b} @sort_r_jobs;
	my $r_jobs_numbers = scalar(@sort_r_jobs);
	if($r_jobs_numbers >= $MAX_R_JOBS){
		$qsub .= ("-hold_jid ".$sort_r_jobs[$r_jobs_numbers - $MAX_R_JOBS]);
	}
    my $jobID;
    LIST: {
    $jobID = `$qsub $sh`;
    print "$jobID";
    $jobID = (split /\s+/,$jobID)[2];
    if ($jobID !~/^\d+$/) {
            print "qsub has not been done for $sh ...\n";
            sleep 120;
            redo LIST;
            }
    }
    push @{$JOB_HASH{$jobID}},$sh,$memory,$qsub_queue;
	push @{$JOB_HASH{SHELL}},$sh if($record == 1);
	return $jobID;
}
#-------------------- qsub mutil line shell -----------------
#---- split the shell into single line shell, and qsub ------
sub qsubmulti{
	&check_user if($CHECK_USER);
	my %para_hash;
	&qsubmulti_para(\%para_hash,@_);
	my ($sh,$memory,$splitline,$queue_para,$subsh_prefix) = ($para_hash{'s'},$para_hash{'m'},$para_hash{'l'},$para_hash{'q'},$para_hash{'sp'});
	chomp(my $time_anchor = `date +\%s`);
	die "Jobcluster::Jobguard::qsubmult has run once, You should run Jobcluster::Jobguard::jobguard\n" if(exists($JOB_HASH{START_TIME}));
	$JOB_HASH{START_TIME} = $time_anchor;
	`mkdir -p $sh.$time_anchor.qsub`;
#	(my $line_sum) = (`wc -l $sh` =~ /^(\d+)\s+\S+$/); ## get the line number of mutil-shell
	(my $line_sum) = (`grep -v '^#' $sh | wc -l` =~ /^(\d+)/); ## get the line number of mutil-shell, but not '^#'
	die "Cannot split $sh, because its executable line number is zero.\n" if($line_sum == 0);
	$splitline = $line_sum if($splitline > $line_sum);
	my $work_num = 1;
	open (OLDSH,"$sh") || die "fail $sh: $!\n";
	my $out_line = '';
	while(<OLDSH>){
		next if(/^\s+$/ || /^\#/);
		s/\s+$//; ## chomp
		$out_line .= ($para_hash{'ct'})?"time ($_) && ":"($_) && ";
		if($. >= ($work_num * $splitline) || $. == $line_sum){
			my $split_shell = "$sh.$time_anchor.qsub/$subsh_prefix".&mendnumber(length($work_num),length(int($line_sum/$splitline)+1)).$work_num.'.sh';
			open (SPLITSH,">$split_shell") || die "fail $split_shell: $!\n";
#			$out_line =~ s/;/ && /g; ## convert all ';' to ' && '
			$out_line .= ($para_hash{'cu'})?"(perl $USAGE_CHECK_PERL_SCRIPT -u $USER -s $split_shell) && ":'';
			print SPLITSH "${out_line}(perl -e 'print STDERR \"\\nThis-work-is-completed\\n\"')\n";
			close SPLITSH;
			$SINGLE_CHECK_SIGN = 0; ## close the singleline check
			&qsubsingle("-s:$split_shell","-m:$memory","-r:2","-q:$queue_para"); ## unrecord in qsubsingle
			$SINGLE_CHECK_SIGN = 1; ## open the singleline check
			$work_num++;
			$out_line = '';
		}
	}
	close OLDSH;
	push @{$JOB_HASH{SHELL}},$sh;
}
#-------------- mend the number with prefix-'0' ---------------
sub mendnumber{
	my ($length_mend,$length_max) = @_;
	die "Wrong: $length_max < $length_mend\n" if($length_max < $length_mend);
	return ('0' x ($length_max - $length_mend));
}
#-------------- qdel job on the cluster ---------------
#------------- qdel first and then check --------------
sub qdel_job{
	my $dID = shift;
	for (my $time = 1;;$time++){
		`qdel $dID`; ## qdel first and then check
		sleep 10; ## wait qdel result
		my $info = `qstat -j $dID 2>&1 | head -2`;
		if($info =~ /Following\s+jobs\s+do\s+not\s+exist/){ ## qdel successfully ## Following jobs do not exist
			print LOG "INFO:\tError job $dID has been qdel\n";
			return;
		}
		if($time > 2){ ## three times 'qdel' failed
			print LOG "!ERROR!:\tError job $dID cannot qdel in $time times\n\n";
			return;
		}
		sleep 60;
	}
}
#------- test wheather store the new error job's shell --------
sub test_store{
	my ($error_array,$error_shell) = @_;
	for (my $i=0; $i < @$error_array; $i+=2) {
		return 0 if($$error_array[$i] eq $error_shell);
	}
	return 1;
}
#------- deal the para of routine jobguard ---------
sub jobguard_para{
	my $phash = shift;
	my %test;
	$test{$_} = 1 for ('d','t','r','n','e','h','ms','mc','me','mr','ql','df'); # $log_dir,$interval,$reqsub,$logname_prefix,$max_eqw_time,$qhost_time_span,$mem_check_sign,$mem_cycle,$mem_exceed,$mem_record,$queue_limit,check-disk-free
	foreach my $para (@_) {
		my ($sign,$value) = ($para =~ /^\s*\-+([^\s:]+)\s*[:=]\s*(\S+)\s*$/);
		die"Para_sign:$sign is not required.\n$para is wrong!\n" unless($sign && $value && exists($test{$sign})); # check the para_sign existence
		if(exists($$phash{$sign})){
			die "$sign is not a multi parameter, you can not give it several values!\n" if($sign ne 'df'); ## special for multi parameter 'df'
			$$phash{$sign} .= ';'.$value;
		}
		else{
			$$phash{$sign} = $value; # upload the para_hash
		}
	}
	#------ check the log_dir --------
	die "Your para should contain '-d:log_directory'\n" unless($$phash{'d'});
	`mkdir -p $$phash{'d'}` unless(-d $$phash{'d'});
	$$phash{'d'} = abs_path($$phash{'d'});
	#------ check the interval -------
	$$phash{'t'} ||= 300;
	$$phash{'t'} =~ s/[a-z][A-Z]//g;
	$$phash{'t'} = 300 if($$phash{'t'} !~ /^\d+$/ || $$phash{'t'} < 120);
	#------- check the reqsub sign -----
	$$phash{'r'} ||= 0; ## default as not reqsub error-jobs
	#------- check the max eqw time -----
	$$phash{'e'} = 60 if(!$$phash{'e'} || $$phash{'e'} !~ /^\d+$/ || $$phash{'e'} > 100);
	#------- check the qhost time span -----
	$$phash{'h'} = 1 if(!$$phash{'h'} || $$phash{'h'} !~ /^\d+$/ || $$phash{'h'} > 10);
	#------- check mem_check_sign --------
	$$phash{'ms'} ||= 0; ## default as close memort check
	#------- check mem_cycle -------
	die "Cannot set memory_check_cycle while memory_check_function closed!\n" if(!$$phash{'ms'} && exists($$phash{'mc'}));
	$$phash{'mc'} = 5 if(!$$phash{'mc'} || $$phash{'mc'} !~ /^\d+$/ || !($$phash{'mc'} > 0));
	#------- check mem_exceed -------
	die "Cannot set memory_exceed while memory_check_function closed!\n" if(!$$phash{'ms'} && exists($$phash{'me'}));
	$$phash{'me'} ||= 0.2; ## default
	if($$phash{'me'} =~ /^[\d\.]+$/){ ## pecentage mode
		$$phash{'me'} =~ s/\.+$//;
		$$phash{'me'} = 0.2 unless($$phash{'me'} >= 0 && $$phash{'me'} <= 1);
	}
	elsif($$phash{'me'} =~ /^([\d\.]+)([mg])$/i){ ## just the specific memory
		$$phash{'me'} = &get_m_memory($$phash{'me'}).'m'; ## in unit 'm' for easily memory-compare
	}
	else{ ## bad memory_exceed input
		die "memory_exceed: $$phash{'me'} is not right form!\nOnly decimals or specific memory in 'M/m/G/g' unit!\n";
	}
	#------- check memory_record_sign -------
	die "Cannot set memory_record while memory_check_function closed!\n" if(!$$phash{'ms'} && exists($$phash{'mr'}));
	$$phash{'mr'} = 1 if(!exists($$phash{'mr'})); ## default as open memory record while memory_check_function opens
	#------- job queue limit sign ---------
	$$phash{'ql'} = 1 unless(exists($$phash{'ql'})); ## default queue-limit
	#------- disk free space -------
	if(exists($$phash{'df'})){
		@{$DEVICE{abs_path((split /,/)[0])}} = ((split /,/)[1,2],0) for (split /;/,$$phash{'df'});
		foreach my $disk (keys %DEVICE) {
			my ($min,$rls) = @{$DEVICE{$disk}}[0,1];
			$min = &transform_space($min);
			$rls = &transform_space($rls);
			if($min + 50 > $rls){
				my @jobs = keys %JOB_HASH;
				`qdel @jobs`;
				die "\nThe rls_disk_space should be more than min_disk_space at least 50G for safty!\nBut your input is not!\ndisk: $disk\nmin_disk_space: ${$DEVICE{$disk}}[0]\nrls_disk_space: ${$DEVICE{$disk}}[1]\nQdel all jobs!\nPlease check your jobs by hand!\n\n";
			}
		}
	}
}
#------- deal the para of routine qsubmulti ---------
sub qsubmulti_para{
	my $phash = shift;
	my %test;
	$test{$_} = 1 for ('s','m','l','q','ct','cu','sp'); # $sh,$memory,$splitline,$qsub_queue,check_job_time,check_job_usage,$subsh_prefix
	foreach my $para (@_) {
		my ($sign,$value) = ($para =~ /^\s*\-+([^\s:]+)\s*[:=]\s*(\S+)\s*$/);
		die "Cannot identify the para: $para\n" unless($sign && $value);
		die "Para_sign:$sign is not required.\n$para is wrong!\n" unless(exists($test{$sign})); # check the para_sign existence
		$$phash{$sign} = $value; # upload the para_hash
	}
	#------ check the shellpath --------
	die "Your para should contain '-s:shell_path'\n" unless($$phash{'s'});
	die "Cannot find the shellpath: $$phash{'s'}\n" unless(-e $$phash{'s'});
	$$phash{'s'} = abs_path($$phash{'s'});
	#------ check the interval -------
	die "Your para should contain '-m:memory'\n" unless($$phash{'m'});
	if($$phash{'m'} =~ /^\d+$/){
		print STDERR "WARNING: Your memory for $$phash{'s'} has been reset as 1G\nBecause your input $$phash{'m'} has no unit!\n";
		$$phash{'m'} = '1G';
	}
	#------- check the splitline -----
	if(!$$phash{'l'} || $$phash{'l'} !~ /^\d+$/){
		$$phash{'l'} = 1;
	}
	#------- check the submit queue -------
	$$phash{'q'} ||= 'NA';
	#------- check the subsh prefix -------
	$$phash{'sp'} ||= 'work';
}
#------- deal the para of routine qsubsingle ---------
sub qsubsingle_para{
	my $phash = shift;
	my %test;
	$test{$_} = 1 for ('s','m','r','c','q'); # $sh,$memory,$record,$change_single_sign,$qsub_queue
	foreach my $para (@_) {
		my ($sign,$value) = ($para =~ /^\s*\-+([^\s:]+)\s*[:=]\s*(\S+)\s*$/);
		die"Para_sign:$sign is not required.\n$para is wrong!\n" unless($sign && $value && exists($test{$sign})); # check the para_sign existence
		$$phash{$sign} = $value; # upload the para_hash
	}
	#------ check the shellpath --------
	die "Your para should contain '-s:shell_path'\n" unless($$phash{'s'});
	die "Cannot find the shellpath: $$phash{'s'}\n" unless(-e $$phash{'s'});
	$$phash{'s'} = abs_path($$phash{'s'});
	#------ check the interval -------
	die "Your para should contain '-m:memory'\n" unless($$phash{'m'});
	if($$phash{'m'} =~ /^\d+$/){
		print STDERR "WARNING: Your memory for $$phash{'s'} has been reset as 1G\nBecause your input $$phash{'m'} has no unit!\n";
		$$phash{'m'} = '1G';
	}
	#------- check the splitline -----
	if(!$$phash{'r'} || $$phash{'r'} !~ /^\d+$/){
		$$phash{'r'} = 1;
	}
	#------- check the change shell sign ------
	$$phash{'c'} ||= 0;
	#------- check the submit queue -------
	$$phash{'q'} = (exists($$phash{'q'}) && $$phash{'q'} =~ /(\S+)\.q$/)?"$$phash{'q'}":'NA';
}
#---------- get the memory in unit 'm' -------
sub get_m_memory{
	my $ele = shift;
	my ($number,$unit) = ($ele =~ /^([\d\.]+)([mg])$/i);
	die "$ele is not the right format of Jobcluster::Jobguard::get_m_memory\n" unless($number && $unit);
	$number =~ s/\.+$//; ## do safe
	if($unit eq 'm'){
		return $number;
	}
	elsif($unit eq 'M'){
		return ($number * 1.024);
	}
	elsif($unit eq 'g'){
		return ($number * 1000);
	}
	elsif($unit eq 'G'){
		return ($number * 1024);
	}
	else{
		die "$ele is not the right format of Jobcluster::Jobguard::get_m_memory\n";
	}
}
#-------- jobs' memory guard --------
sub memory_guard{
	my ($jobID,$mem_exceed) = @_[0,1];
	my $usage;
	my $qj_time = 0;
	UASGE:{
		$usage = `qstat -j $jobID | grep usage`; ## stderr to stdout
		$qj_time++;
		# usage    1:                 cpu=02:23:06, mem=4481.25708 GBs, io=2.67408, vmem=602.082M, maxvmem=602.082M
		unless($usage =~ /^usage/ && $usage =~ /cpu=[^,]+,\s+mem=[^,]+,\s+io=[^,]+,\s+vmem=[^,]+,\s+maxvmem=/){
			if($qj_time > 3){
				print LOG "\n<WARNING>:\tThe shell order 'qstat -j $jobID | grep usage' gets invaid info for 3 times! Skipping it...\n\n";
				return ('-1',1,0); ## memory=-1, memory is ok, not qdel
			}
			if((my $new = `qstat -j $jobID 2>&1`) =~ /Following\s+jobs\s+do\s+not\s+exist/){ ## job not exists
				print LOG "\n<WARNING>:\t$jobID does not exist! Skipping it...\n\n";
				return ('-1',1,0); ## memory=-1, memory is ok, not qdel
			}
			sleep 30;
			print LOG "\n" if($qj_time == 1);
			print LOG "<WARNING>:\tThe shell order 'qstat -j $jobID | grep usage' gets invaid info! Redoing it Now!\n";
			redo UASGE;
		}
	}
	chomp($usage);
	my @usage = split /\n+/,$usage;
	my $max_usage;
	foreach my $usage_info (@usage) {
		($max_usage) = ($usage_info =~ /maxvmem=(\S+)/);
		if($max_usage =~ /N\/?A/i){
			return ('NA',1,0); ## memory=NA, memory is ok, not qdel
		}
		else{
			my $used_max_memory = &get_m_memory($max_usage);
			my $required_memory = &get_m_memory(${$JOB_HASH{$jobID}}[1]);
			my $limit_upper_memory;
			if($mem_exceed =~ /m/i){
				(my $mem_excd_number = $mem_exceed) =~ s/m//i;
				$limit_upper_memory = $required_memory + $mem_excd_number;
			}
			else{
				$limit_upper_memory = $required_memory * (1 + $mem_exceed);
			}
			if($used_max_memory > $limit_upper_memory){
				return ("${used_max_memory}m",0,1); ## memory=used_max_memory, memory is not ok, qdel
			}
		}
	}
	return ($max_usage,1,0); ## memory=used_max_memory, memory is ok, not qdel
}
#-------- feedback the user name, ten seconds to kill this backstage task ------
sub check_user{
	my $wait_time = 15;
	print STDERR "\n".('#' x 26)." <NOTE> ".('#' x 26)."\nUse `whoami` and get your userID is '$USER'.\nIs it right?\nIf not, Your will have $wait_time seconds to kill this running backstage task, OR leave it...\nWant to get your valid(real) userID, please perldoc /share/backup/jiawl/useful/perl_script/module/Jobcluster/Jobguard.pm and see NOTE<2>\n".('#' x 60)."\n";
	$CHECK_USER = 0;
#	sleep $wait_time; ## abandon this
	print STDERR "Counting Down: ";
	for (my $second = $wait_time;$second > 0;$second--) { ## wait to be killed....oh, my god~~~ come on~~~
		print STDERR "$second ";
		sleep 1;
	}
	print STDERR "\nOK! The USERID is confirmed, all works start NOW...\n";
	sleep 2;
}
#------ check the disk free space info -------
sub disk_free{
	my ($anchor) = $_[0];
	my $return = 1;
	foreach my $disk (keys %DEVICE) {
		my $disk_avail = &get_disk_info($disk);
		my $anchor = ${$DEVICE{$disk}}[2];
		if(&compare_space($disk_avail,${$DEVICE{$disk}}[$anchor]) eq 'NA'){ ## compare the space 'NA'
			print LOG "\n<WARNING>:\tAvial space ($disk_avail) of disk $disk is less than the minimum amount(${$DEVICE{$disk}}[$anchor])." if($anchor == 0);
			$return = 0;
			${$DEVICE{$disk}}[2] = 1; ## this disk is not ok
		}
		else{ ## avail space OK
			print LOG "Avial space ($disk_avail) of disk(ever bad) $disk has reached the release amount(${$DEVICE{$disk}}[$anchor])." if($anchor == 1);
			${$DEVICE{$disk}}[2] = 0; ## this disk is ok
		}
		sleep 5; ## have a rest
	}
	return $return;
}
#------ get the disk free space info -------
sub get_disk_info{
	my ($disk) = $_[0];
	my $time = 0;
	my $avail;
	DF: 
	{
		my $df_info = `df -h $disk`;
		$time++;
		if($df_info !~ /Filesystem.+Size.+Used.+Avail/i){
			die "Cannot get the disk $disk space info!\n" if($time == 3);
			redo DF;
		}
		else{
			($df_info) = (split /\n+/,$df_info)[2];
			$df_info =~ s/^\s+//g;
			$df_info =~ s/\s+$//g;
			($avail) = (split /\s+/,$df_info)[2];
			if($avail !~ /\d/){
				die "Cannot get the avail space info of disk $disk!\n" if($time == 3);
				redo DF;
			}
		}
	}
	return $avail;
}
#------ compare the spaces -------
sub compare_space{
	my ($want_bigger,$want_smaller) = @_[0,1];
	$want_bigger = &transform_space($want_bigger);
	$want_smaller = &transform_space($want_smaller);
	if($want_bigger >= $want_smaller){
		return 'OK';
	}
	else{
		return 'NA';
	}
}
#------ transform space based on the unit --------
sub transform_space{
	my ($old) = $_[0];
	my ($num,$unit) = ($old =~ /^([\d\.]+)(.+)?$/);
	$unit ||= 'k'; ## default
	my (%multiple,@multiple);
	push @multiple,(1024**($_)) for (-2 .. 1);
	@multiple{qw/T t G g M m K k/} = @multiple[3,3,2,2,1,1,0,0];
	if(exists($multiple{$unit})){
		$num *= $multiple{$unit};
	}
	else{
		die "Cannot distinguish the unit $unit from $old!\n";
	}
	return $num;
}
#------ hold jobs --------
sub hold_jobs{
	my @running_jobs = sort {$a<=>$b} keys %JOB_HASH;
	my @qhold_jobs = `qmod -s @running_jobs | grep 'can not be' | cut -d ' ' -f10 | cut -d '.' -f1`; ## first qmod -s running jobs
	if(@qhold_jobs){
		chomp(@qhold_jobs);
		`qhold @qhold_jobs`; ## then qhold other jobs in others status
	}
}
#------ release jobs --------
sub release_jobs{
	my @running_jobs = sort {$a<=>$b} keys %JOB_HASH;
	my @qrls_jobs = `qmod -us @running_jobs | grep 'can not be' | cut -d ' ' -f10 | cut -d '.' -f1`; ## first qmod -us s jobs
	if(@qrls_jobs){
		chomp(@qrls_jobs);
		`qrls @qrls_jobs`; ## then qhold other jobs in others status
	}
}

1; ## tell the perl script the successful access of this module.
