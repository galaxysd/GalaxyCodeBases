#!/usr/bin/perl -w

use strict;
use File::Basename;

unless (@ARGV > 0) {
	print "perl $0 <SNP_dir> <qality> <Q20_dir>\n";
	exit 0;
}

my $filter_bin = './snpFilter.py';#"/nas/RD_09C/resequencing/soft/Pipeline/SNPcalling/subBin/snpFilter.py";

my $snp_dir = $ARGV[0];
my $q_val = $ARGV[1];
my $q18_dir = $ARGV[2];
mkdir "$q18_dir" unless (-d "$q18_dir");

my @snp = `find $snp_dir -name '*.snp'`;
chomp @snp;

foreach my $snp(@snp) {
	my $basename = basename $snp;
	open SH, ">$q18_dir/$basename.sh" || die "$!\n";
	print SH "#!/bin/sh
#\$ -N \"SNP$basename\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=50M,p=1
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash
python $filter_bin $snp $q_val 2 5 3 100 1 > $q18_dir/$basename.Q$q_val
";	# 2 5 1 40 1
	close SH;
	#`qsub -cwd -l vf=1g $q18_dir/$basename.sh`;
	#sleep 1;
}
