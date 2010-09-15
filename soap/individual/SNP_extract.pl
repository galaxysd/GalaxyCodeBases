#!/usr/bin/perl -w

use strict;
use File::Basename;

unless (@ARGV > 0) {
	print "perl $0 <CNS_dir> <SNP_dir>\n";
	exit 0;
}

my $cns_dir = $ARGV[0];
my $snp_dir = $ARGV[1];
mkdir "$snp_dir" unless (-d "$snp_dir");

my $bin = './snp_extract.py';#"/nas/RD_09C/resequencing/soft/Pipeline/SNPcalling/subBin/snp_extract.py";

my @cns = `find $cns_dir -name '*.cns'`;
chomp @cns;

foreach my $cns (@cns) {
	my $basename = basename $cns;
	if ($basename =~ /(.*)\.cns$/) {
		open SH, ">$snp_dir/$1.sh" || die "$!\n";
		print SH "#!/bin/sh
#\$ -N \"SNP$1\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=100M,p=1
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash
python $bin $cns > $snp_dir/$1.snp
";
		close SH;
		#`qsub -cwd -l vf=1G $snp_dir/$1.sh`;
		#sleep 1;
	}
}
