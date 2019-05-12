package main;
use strict;
use warnings;
use Data::Dump qw(ddx);

my $SgeQueue = 'fgi.q';
my $SgeProject = 'fgi';
my $SgeJobPrefix = 'nip'.(time() % 999);

sub Scutadapt($$$$) {
	my ($cwd,$cnt,$flst,$outP) = @_;
	my $SHcutadapt = <<"END_SH";
#!/bin/sh
#\$ -S /bin/bash
#\$ -q $SgeQueue -P $SgeProject
#\$ -N p0${SgeJobPrefix}cut
#\$ -l vf=600M,num_proc=2
#\$ -binding linear:3
#\$ -cwd -r y
#\$ -v PERL5LIB,PATH,LD_LIBRARY_PATH
#\$ -e /dev/null -o /dev/null

#\$ -t 1-$cnt

cd $cwd

INFILE=`sed -n "\${SGE_TASK_ID}p" $flst`
read -ra INDAT <<<"\$INFILE"

ADAPTERSTR='-g GAACGACATGGCTACGATCCGACTT -a AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -A AAGTCGGATCGTAGCCATGTCGTTC -G TTGTCTTCCTAAGACCGCTTGGCCTCCGACTT'

if [ -f "\${INDAT[1]}.fq.gz" ]; then
	cutadapt -j 2 -m 1 -O 4 -n 2 \$ADAPTERSTR --interleaved -o $outP/\${INDAT[0]}.fq.gz \${INDAT[1]}.fq.gz >$outP/\${INDAT[0]}.log
else
	cutadapt -j 2 -m 1 -O 4 -n 2 \$ADAPTERSTR --interleaved -o $outP/\${INDAT[0]}.fq.gz \${INDAT[1]}_1.fq.gz \${INDAT[1]}_2.fq.gz >$outP/\${INDAT[0]}.log
fi

END_SH
	return $SHcutadapt;
}

1;
