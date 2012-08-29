#!/bin/sh
#
#  http://web.sevensdoor.com/download/fmsbackup.sh.txt
#
#  FileMaker Server Backup Script ver.1.5.2
#  Backup FileMaker database files latest n hours and latest n days.
#  This script is written for FileMaker Server 7-11 on Mac OS X.
#
#  Copyright (C) 2007 Koji Takeuchi
#
#  This script is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#  
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  Koji Takeuchi <takeuchi@nemoux.com>
#  2007.09.20
#
#
#  To use this script:
#  0. Set permission of this file to executable.
#     ex. $ chgrp fmsadmin fmsbackup.sh
#         $ chmod 750 fmsbackup.sh
#  1. You can edit 'set variables' section of this file directly, then
#     run this script.
#     OR
#  2. You can separate 'set variables' section to other files as a 'conf' file.
#     and run this script with arguments to include variables specified on 
#     separate 'conf' files.
#     ex. $ ./fmsbackup.sh hourly.conf
#         $ /usr/local/bin/fmsbackup.sh ~/Documents/daily.conf
#
#     If you specify the same variables on both this file and the separate file,
#     this file's variables are replaced by variables loaded from the
#     separate file.
#
#  To run this script periodically:
#  1. Put this script on '/Library/FileMaker Server/Data/Scripts', then
#     use 'System-Level Script' schedule on FileMaker Server 9
#     Admin Console.
#     OR
#  2. Install your crontab.
#     ex: $ crontab -e
#         03 * * * * /usr/local/bin/fmsbackup.sh ~/Documents/hourly.conf
#         30 04 * * * /usr/local/bin/fmsbackup.sh ~/Documents/daily.conf
#
#         In this case, hourly schedule will launch every '**:03',
#         and daily schedule will launch every '04:30'.
#         If you lose your mind when you just enter 'crontab -e',
#         see 'man vi'. (or just type ':q!' and forget it.)
#     OR
#  3. If you are Apple Certified System Administrator, use launchd.
#     Figure out a property list file and save as
#     '/Library/LaunchDaemons/com.example.hourlybackup.plist'.
#     And load the property list file.
#     ex: $ sudo launchctl -w load /Library/LaunchDaemons/com.example.hourlybackup.plist
#
#       - contents of property list file should be like this...
#
#     <?xml version="1.0" encoding="UTF-8"?> 
#     <!DOCTYPE plist PUBLIC "-//Apple Computer//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd"> 
#     <plist version="1.0"> 
#     <dict> 
#     	<key>Label</key>	<string>com.example.hourlybackup</string> 
#     	<key>Disabled</key>	<false/> 
#     	<key>ProgramArguments</key> 
#     		<array> 
#     			<string>/usr/local/bin/fmsbackup.sh</string> 
#     			<string>/Users/YOURNAME/Documents/hourly.conf</string> 
#     		</array> 
#     	<key>StartCalendarInterval</key>
#			<dict>
#				<key>Minute</key>  <integer>03</integer>
#			</dict>
#     	<key>RunAtLoad</key>	<true/> 
#     </dict> 
#     </plist>
#
#  To get report:
#     Simply redirect to mail command.
#     You can get (sometimes too many) reports via mail.
#     ex: /usr/local/bin/fmsbackup.sh ~/Documents/daily.conf | mail -s "fms backup report" yourname@example.com
#
#
#
#
################################################
#                set variables                 #
################################################
#############  Required variables  #############
# user name for fms admin console
user=USERNAME
#
# password for fms admin console
pass=PASSWORD
#
# 'data' path ( parent of Database )
data=/Library/'FileMaker Server'/Data
#       
# Database directory path
databases="$data"/Databases
#
# backup directory name
backupdir=Databases
#
# backup directory path
backups=/Library/'FileMaker Server'/Data/Backups
#
# backup Database directory path
backupdatabases="$backups"/"$backupdir"
#
# backup Database directory path for 'fmsadmin backup' command
backupspath="filemac:/'Macintosh HD'$backupdatabases/"
#
# with sub directory? ( y = 1, n = 0 )
subdirectory=1
#
# verify option? ( y = -x, n = "" ) for ver.11 or above
verify=
#
# clone option? ( y = -n, n = "" ) for ver.11 or above
clone=
#
# prefix of archive file name
prefix=servername_fms11
#
# zip or tar.gz? ( zip/tar.gz )
# If you don't need any archiving, just leave it blank
archive=tar.gz
#
#
#
##########  options for archiving   ###########
#  Archive backup files and keep them latest xx hours ( and xx days )
#  If necessary, set 'hours' or 'days' value to n and edit following variables
#
# how many hours do you want to keep shortest terms backups? ( y = n, n = 0 )
#hour_short=5
#
#  If hour=n, set following variables
#
# hourly generation directory name
#shortlydirectory=archives_shortly
#
# hourly generation directory path
#shortlypath="$backups"/"$shortlydirectory"
#
#
# how many hours do you want to keep hourly backups? ( y = n, n = 0 )
#hour=72
#
#  If hour=n, set following variables
#
# hourly generation directory name
#hourlydirectory=archives_hourly
#
# hourly generation directory path
#hourlypath="$backups"/"$hourlydirectory"
#
#
# how many days do you want to keep daily backups? ( y = n, n = 0 )
#day=30
#
#  If day=n, set following variables
#
# daily generation directory name
#dailydirectory=archives_daily
#
# daily generation directory path
#dailypath="$backups"/"$dailydirectory"
#
#
# how many months do you keep monthly backups? ( y = n, n = 0 )
#month=120
#
#  If month=n, set following variables
#
# monthly generation directory name
#monthlydirectory=archives_monthly
#
# monthly generation directory path
#monthlypath="$backups"/"$monthlydirectory"
#
#
#
########  options for external volume  #########
#  Syncronize each backup directory to external backup volume
#  If necessary, set 'ext_*ly' value to 1 and edit following variables
#
# mirror shortly directory to external volume? ( y = 1, n = 0 )
#ext_shortly=0
#
# mirror hourly directory to external volume? ( y = 1, n = 0 )
#ext_hourly=0
#
# mirror daily directory to external volume? ( y = 1, n = 0 )
#ext_daily=0
#
# mirror monthly directory to external volume? ( y = 1, n = 0 )
#ext_monthly=0
#
#  If ext_*ly=1, set following variables
#
# external backup directory path
#external=/Volumes/backup_external/fms/
#
#
#
####  options for alternate backup server 1  #####
#  Syncronize each backup directory to other host via network
#  If necessary, set 'altn_*ly' value to 1 and edit the following variables
#  To use this option, backup server needs the RSA key of FM Server
#
# mirror shortly directory to backup server 1? ( y = 1, n = 0 )
#alt1_shortly=0
#
# mirror hourly directory to backup server 1? ( y = 1, n = 0 )
#alt1_hourly=0
#
# mirror daily directory to backup server 1? ( y = 1, n = 0 )
#alt1_daily=1
#
# mirror monthly directory to backup server 1? ( y = 1, n = 0 )
#alt1_monthly=1
#
#  If altn_*ly=1, set following variables
#
# user name of alternate backup server 1
#host1user=root
#
# host name or IP address of alternate backup server 1
#host1=host.example.com
#
# backup directory path of alternate backup server 1
#host1path=/Volumes/backup_external/foo/fms
#
#
#
####  options for alternate backup server 2  #####
#  Syncronize each backup directory to other host via network
#  If necessary, set 'altn_*ly' value to 1 and edit the following variables
#  To use this option, backup server needs the RSA key of FM Server
#
# mirror shortly directory to backup server 2? ( y = 1, n = 0 )
#alt2_shortly=0
#
# mirror hourly directory to backup server 2? ( y = 1, n = 0 )
#alt2_hourly=0
#
# mirror daily directory to backup server 2? ( y = 1, n = 0 )
#alt2_daily=1
#
# mirror monthly directory to backup server 2? ( y = 1, n = 0 )
#alt2_monthly=1
#
#  If altn_*ly=1, set following variables
#
# user name of alternate backup server 2
#host2user=root
#
# host name or IP address of alternate backup server 2
#host2=host.example.com
#
# backup directory path of alternate backup server 2
#host2path=/Volumes/backup_external/foo/fms
#
#
#
################################################
#                load variables                #
################################################
#######  set variables from other file  ########
#  If run this script with arguments, some of the variables you set above
#  are replaced by contents of file specified by argument
#
# DO NOT EDIT following 5 lines
includes="$1"
#
if [ -r "$includes" ]; then
	. "$includes"
fi
#
################################################
#                fmserver backup               #
################################################
FMSALIVE=$(ps axww | grep fmserverd | grep -v grep)
ISCONSOLEARRIVE=$(/usr/bin/fmsadmin list files -u "$user" -p "$pass" | grep [Ff][Pp]7 | grep -v grep)
if [ -z "$FMSALIVE" ]; then
	echo "FileMaker Server is not working."
	echo "exit script."
	exit 1
fi
if [ -n "$ISCONSOLEARRIVE" ]; then
        echo "fmsadmin is alive."
else
        echo "is fmsadmin dead?\n\n$ISCONSOLEARRIVE" | mail -s "`hostname` fmsadmin failure" takeuchi@splash.jp
        exit 1
fi
echo "-------------------------------------------------------------------"
echo ""
echo "          `hostname` FileMaker Server Backup Report!"
echo ""
echo "-------------------------------------------------------------------"
echo ""
echo ""
echo ""
echo "1.  Backup database files"
#
#   Delete old files in temporary backup directory
rm -rf "$backupdatabases"/*
#
#	Synchronize source directory to temporary backup directory
cd "$data"
/usr/bin/fmsadmin backup "$backupdir"/*.[FfIi][PpNn][7Ss] "$verify" "$clone" -d "$backupspath" -u "$user" -p "$pass"
if [ "$subdirectory" -eq 1 ]; then
	/usr/bin/fmsadmin backup "$backupdir"/*/*.[FfIi][PpNn][7Ss] "$verify" "$clone" -d "$backupspath" -u "$user" -p "$pass"
fi
#
#   Get timestamp for archive file name
shortlyfilename="$prefix"_`date +"%Y_%m_%d_%H%M%S"`
hourlyfilename="$prefix"_`date +"%Y_%m_%d_%H%M%S"`
dailyfilename="$prefix"_`date +"%Y_%m_%d"`
monthlyfilename="$prefix"_`date +"%Y_%m"`
#
echo "    ... Done!"
echo ""
#
#
echo "2.  Deleting old files"
#
#	Delete old files in n hours directory and n days directory
if [ "$hour_short" -gt 0 ]; then
	cd "$shortlypath"
	find . -type f -mmin +`expr "$hour_short" \* 60` | xargs rm -fv
fi
if [ "$hour" -gt 0 ]; then
	cd "$hourlypath"
	find . -type f -mmin +`expr "$hour" \* 60` | xargs rm -fv
fi
if [ "$day" -gt 0 ]; then
	cd "$dailypath"
	find . -type f -mtime +"$day" | xargs rm -fv
fi
echo "    ... Done!"
echo ""
#
#
echo "3.  Archiving files"
#
#	Archive temporary backup directory to zip or tar ball
cd "$backups"
if [ "$archive" = "zip" ]; then
	zip -r "$prefix".zip "$backupdir"
elif [ "$archive" = "tar.gz" ]; then
	tar czvpf "$prefix".tar.gz "$backupdir"
fi
#
#	Copy archive file to 'nhours directory/prefix_yyyy_MM_dd_hhmmss.tar.gz'
#	and 'ndays directory/prefix_yyyy_MM_dd.tar.gz'
#   and 'nmonths directory/prefix_yyyy_MM.tar.gz'
if [ "$hour_short" -gt 0 ]; then
	cp -fp "$prefix"."$archive" "$shortlydirectory"/$shortlyfilename."$archive" 
fi
if [ "$hour" -gt 0 ]; then
	cp -fp "$prefix"."$archive" "$hourlydirectory"/$hourlyfilename."$archive" 
fi
if [ "$day" -gt 0 ]; then
	cp -fp "$prefix"."$archive" "$dailydirectory"/$dailyfilename."$archive"
fi
if [ "$month" -gt 0 ]; then
	cp -fp "$prefix"."$archive"  "$monthlydirectory"/$monthlyfilename."$archive"
fi
echo "    ... Done!"
echo ""
#
#
## options for external volume ##
#
if [ "$ext_shortly" -eq 1 ] || [ "$ext_hourly" -eq 1 ] || [ "$ext_daily" -eq 1 ] || [ "$ext_monthly" -eq 1 ]; then
	echo "4.  Syncronizing archive files to external drive"
#
#
#	Syncronize each backup directory to alternate backup directory
	if [ "$ext_shortly" -eq 1 ]; then
		rsync -vrpo --delete -u --ignore-existing "$shortlydirectory" "$external"
	fi
	if [ "$ext_hourly" -eq 1 ]; then
		rsync -vrpo --delete -u --ignore-existing "$hourlydirectory" "$external"
	fi
	if [ "$ext_daily" -eq 1 ]; then
		rsync -vrpo --delete -u --ignore-existing "$dailydirectory" "$external"
		/bin/cp -vp "$dailydirectory"/$dailyfilename."$archive" "$external"/"$dailydirectory"/
	fi
	if [ "$ext_monthly" -eq 1 ]; then
		rsync -vrpo --delete -u --ignore-existing "$monthlydirectory" "$external"
		/bin/cp -vp "$monthlydirectory"/$monthlyfilename."$archive" "$external"/"$monthlydirectory"/
	fi
	currentnum=4
	echo "    ... Done!"
	echo ""
else
	currentnum=3
fi
#
#
## options for alternate backup server 1 ##
#
if [ "$alt1_shortly" -eq 1 ] || [ "$alt1_hourly" -eq 1 ] || [ "$alt1_daily" -eq 1 ] || [ "$alt1_monthly" -eq 1 ]; then
	echo `expr "$currentnum" \+ 1`".  Transfering archive files to alternate backup server 1";
#
#
#	Syncronize each backup directory to other host via network
	if [ "$alt1_shortly" -eq 1 ]; then
		rsync -vrpo --delete -u --ignore-existing --copy-links -e ssh "$shortlydirectory" "$host1user"@"$host1":"$host1path"/
	fi
	if [ "$alt1_hourly" -eq 1 ]; then
		rsync -vrpo --delete -u --ignore-existing --copy-links -e ssh "$hourlydirectory" "$host1user"@"$host1":"$host1path"/
		/usr/bin/scp -vp "$hourlydirectory"/$hourlyfilename."$archive" "$host1user"@"$host1":"$host1path"/"$hourlydirectory"/
	fi
	if [ "$alt1_daily" -eq 1 ]; then
		rsync -vrpo --delete -u --ignore-existing --copy-links -e ssh "$dailydirectory" "$host1user"@"$host1":"$host1path"/
		/usr/bin/scp -vp "$dailydirectory"/$dailyfilename."$archive" "$host1user"@"$host1":"$host1path"/"$dailydirectory"/
	fi
	if [ "$alt1_monthly" -eq 1 ]; then
		rsync -vrpo --delete -u --ignore-existing --copy-links -e ssh "$monthlydirectory" "$host1user"@"$host1":"$host1path"/
		/usr/bin/scp -vp "$monthlydirectory"/$monthlyfilename."$archive" "$host1user"@"$host1":"$host1path"/"$monthlydirectory"/
	fi
echo "    ... Done!"
echo ""
echo ""
fi
#
#
## options for alternate backup server 2 ##
#
if [ "$alt2_shortly" -eq 1 ] || [ "$alt2_hourly" -eq 1 ] || [ "$alt2_daily" -eq 1 ] || [ "$alt2_monthly" -eq 1 ]; then
        echo `expr "$currentnum" \+ 1`".  Transfering archive files to alternate backup server 2";
#
#
#       Syncronize each backup directory to other host via network
        if [ "$alt2_shortly" -eq 1 ]; then
                rsync -vrpo --delete -u --ignore-existing --copy-links -e ssh "$shortlydirectory" "$host2user"@"$host2":"$host2path"/
        fi
        if [ "$alt2_hourly" -eq 1 ]; then
                rsync -vrpo --delete -u --ignore-existing --copy-links -e ssh "$hourlydirectory" "$host2user"@"$host2":"$host2path"/
				/usr/bin/scp -vp "$hourlydirectory"/$hourlyfilename."$archive" "$host2user"@"$host2":"$host2path"/"$hourlydirectory"/
        fi
        if [ "$alt2_daily" -eq 1 ]; then
                rsync -vrpo --delete -u --ignore-existing --copy-links -e ssh "$dailydirectory" "$host2user"@"$host2":"$host2path"/
		/usr/bin/scp -vp "$dailydirectory"/$dailyfilename."$archive" "$host2user"@"$host2":"$host2path"/"$dailydirectory"/
        fi
        if [ "$alt2_monthly" -eq 1 ]; then
                rsync -vrpo --delete -u --ignore-existing --copy-links -e ssh "$monthlydirectory" "$host2user"@"$host2":"$host2path"/
		/usr/bin/scp -vp "$monthlydirectory"/$monthlyfilename."$archive" "$host2user"@"$host2":"$host2path"/"$monthlydirectory"/
        fi
echo "    ... Done!"
echo ""
echo ""
fi
#
#
## Report what was done ##
#
echo "Now $databases is copied to"
echo "    $backupdatabases"
if [ "$hour_short" -gt 0 ]; then
	echo "and archived to $shortlypath/$shortlyfilename.$archive."
fi
if [ "$hour" -gt 0 ]; then
	echo "and archived to $hourlypath/$hourlyfilename.$archive."
fi
if [ "$day" -gt 0 ]; then
	echo "and archived to $dailypath/$dailyfilename.$archive."
fi
if [ "$month" -gt 0 ]; then
	echo "and archived to $monthlypath/$monthlyfilename.$archive."
fi
if [ "$ext_shortly" -eq 1 ]; then
	echo "and $shortlypath directory is mirrored to"
	echo "    $external"
fi
if [ "$ext_hourly" -eq 1 ]; then
	echo "and $hourlypath directory is mirrored to"
	echo "    $external"
fi
if [ "$ext_daily" -eq 1 ]; then
	echo "and $dailypath directory is mirrored to"
	echo "    $external"
fi
if [ "$ext_monthly" -eq 1 ]; then
	echo "and $monthlypath directory is mirrored to"
	echo "    $external"
fi
if [ "$alt1_shortly" -eq 1 ]; then
	echo "and $shortlypath directory is mirrored to"
	echo "    $host1:$host1path/$shortlydirectory"
fi
if [ "$alt1_hourly" -eq 1 ]; then
	echo "and $hourlypath directory is mirrored to"
	echo "    $host1:$host1path/$hourlydirectory"
fi
if [ "$alt1_daily" -eq 1 ]; then
	echo "and $dailypath directory is mirrored to"
	echo "    $host1:$host1path/$dailydirectory"
fi
if [ "$alt1_monthly" -eq 1 ]; then
	echo "and $monthlypath directory is mirrored to"
	echo "    $host1:$host1path/$monthlydirectory"
fi
if [ "$alt2_shortly" -eq 1 ]; then
        echo "and $shortlypath directory is mirrored to"
        echo "    $host2:$host2path/$shortlydirectory"
fi
if [ "$alt2_hourly" -eq 1 ]; then
        echo "and $hourlypath directory is mirrored to"
        echo "    $host2:$host2path/$hourlydirectory"
fi
if [ "$alt2_daily" -eq 1 ]; then
        echo "and $dailypath directory is mirrored to"
        echo "    $host2:$host2path/$dailydirectory"
fi
if [ "$alt2_monthly" -eq 1 ]; then
        echo "and $monthlypath directory is mirrored to"
        echo "    $host2:$host2path/$monthlydirectory"
fi
echo ""
echo ""
echo ""
echo "FileMaker Server Backup completed."
echo ""
echo ""
echo ""
echo "----------------------  end of Backup Report  ---------------------"
echo ""
echo ""
exit 0
