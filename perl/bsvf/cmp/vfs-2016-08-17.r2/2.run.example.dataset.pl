#!/usr/bin/perl
###########################################################################
##  ViralFusionSeq
##  Software for discovering and annotating viral integration event and
##    fusion transcript
##  
##  Version 1.0 -- November 15, 2012
##  
##  Copyright (C) 2012 by Jing-Woei Li & Raymond Wan, All rights reserved.
##  Contact:  marcoli@cuhk.edu.hk, rwan@cuhk.edu.hk
##  Organization:  Hong Kong Bioinformatics Centre, School of Life Sciences, The
##                 Chinese University of Hong Kong, Shatin, NT,
##                 Hong Kong SAR
##  
##  This file is part of ViralFusionSeq.
##  
##  ViralFusionSeq is free software; you can redistribute it and/or 
##  modify it under the terms of the GNU General Public License 
##  as published by the Free Software Foundation; either version 
##  3 of the License, or (at your option) any later version.
##  
##  ViralFusionSeq is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public 
##  License along with ViralFusionSeq; if not, see 
##  <http://www.gnu.org/licenses/>.
###########################################################################


##  $LastChangedDate: 2012-11-27 19:50:32 +0800 (Tue, 27 Nov 2012) $
##  $LastChangedRevision: 1094 $

print STDERR "You are about to run a reduced set of the HKCI-5a RNA-Seq data\n";
print STDERR "Please come back after 15 minutes\n";
#my $cmd = "perl viral.fusion.pl --config vfs.conf --insertSIZE 350 HKCI5a CL7R_1_VFS.fq CL7R_2_VFS.fq";
#system ($cmd);
my $cmd = "perl viral.fusion.pl --config vfs.conf --insertSIZE 350 HKCI5a_gz CL7R_1_VFS.fq.gz CL7R_2_VFS.fq.gz";
system($cmd);
