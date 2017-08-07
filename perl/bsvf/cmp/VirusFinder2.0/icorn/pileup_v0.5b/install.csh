#!/bin/tcsh
# by ygu@sanger.ac.uk
# 24/06/2008
#

set locPath = `pwd`

cd $locPath/ssaha_pileup/ssaha_pileup; make
cd $locPath/ssaha_pileup/other_codes/get_seqreads; make

cd $locPath
cd ssaha2
set arch=`uname -m`
if(! -f ssaha2-2.3_$arch) then
	echo "Error: $arch system is not supported for this release!"
	echo "Please contact ygu@sanger.ac.uk for further information!"
	exit
endif

cd $locPath
awk -v path=$locPath '/^# set ssaha_pileupProgramPath/{printf "set prog=%s/ssaha_pileup ",path}/^# set ssaha2ProgramPath/{printf "set ssaha2Path=%s/ssaha2 ",path}//{print $0}' pileup.csh_src >! pileup.csh
chmod 755 pileup.csh

