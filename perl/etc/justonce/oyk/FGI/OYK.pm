package main;
use strict;
use warnings;
use Data::Dump qw(ddx);

my $SgeQueue = 'fgi.q';
my $SgeProject = 'fgi';
my $SgeJobPrefix = 'nip'.(time() % 9999);

sub Scutadapt($$$$) {
	my ($cwd,$cnt,$flst,$outP) = @_;
	my $SHcutadapt = <<"END_SH";
#!/bin/bash
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

if [ -t 1 && "x\$SKIP" = 'x1' ]; then
	while read -u 9 THELINE; do
		read -ra INDAT <<<"\$THELINE"
		if [ -e "\${INDAT[1]}.fq.gz" ]; then
			ln -si \${INDAT[1]}.fq.gz $outP/\${INDAT[0]}.fq.gz
		elif [ -e "\${INDAT[1]}_1.fq.gz" ] && [ -e "\${INDAT[1]}_2.fq.gz" ]; then
			echo "[!]Cannot use SKIP mode for PE, using Read 1 only [\${INDAT[1]}_1.fq.gz] !"
			ln -si \${INDAT[1]}_1.fq.gz $outP/\${INDAT[0]}.fq.gz
		else
			echo "[x]Cannot find .fq.gz file(s) for [\${INDAT[1]}(_[12])?.fq.gz] !"
		fi
	done 9<$flst
	exit 0
fi

INFILE=`sed -n "\${SGE_TASK_ID}p" $flst`
read -ra INDAT <<<"\$INFILE"

ADAPTERSTR='-g GAACGACATGGCTACGATCCGACTT -a AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -A AAGTCGGATCGTAGCCATGTCGTTC -G TTGTCTTCCTAAGACCGCTTGGCCTCCGACTT'

if [ -e "\${INDAT[1]}.fq.gz" ]; then
	cutadapt -j 2 -m 1 -O 4 -n 2 \$ADAPTERSTR --interleaved -o $outP/\${INDAT[0]}.fq.gz \${INDAT[1]}.fq.gz >$outP/\${INDAT[0]}.log
elif [ -e "\${INDAT[1]}_1.fq.gz" ] && [ -e "\${INDAT[1]}_2.fq.gz" ]; then
	cutadapt -j 2 -m 1 -O 4 -n 2 \$ADAPTERSTR --interleaved -o $outP/\${INDAT[0]}.fq.gz \${INDAT[1]}_1.fq.gz \${INDAT[1]}_2.fq.gz >$outP/\${INDAT[0]}.log
else
	echo "[x]Cannot find .fq.gz file(s) for [\${INDAT[1]}(_[12])?.fq.gz] !"
	exit 1
fi
END_SH
	return $SHcutadapt;
}

sub Sbwamem($$$$$$) {
	my ($cwd,$cnt,$flst,$outP,$FQprefix,$pRef) = @_;
	my $SHbwa = <<"END_SH";
#!/bin/bash
#\$ -S /bin/bash
#\$ -q $SgeQueue -P $SgeProject
#\$ -N p1${SgeJobPrefix}bwa -hold_jid p0${SgeJobPrefix}cut
#\$ -l vf=6G,num_proc=7
#\$ -binding linear:12
#\$ -cwd -r y
#\$ -v PERL5LIB,PATH,LD_LIBRARY_PATH
#\$ -e /dev/null -o /dev/null

#\$ -t 1-$cnt

cd $cwd

INFILE=`sed -n "\${SGE_TASK_ID}p" $flst`
read -ra INDAT <<<"\$INFILE"

bwa mem -t 12 -Y $pRef -R "\@RG\\tID:\${INDAT[0]}\\tSM:\${INDAT[0]}" -p $FQprefix/\${INDAT[0]}.fq.gz 2>$outP/\${INDAT[0]}.log | samtools view -bS - | samtools sort -m 2097152 -T $outP/\${INDAT[0]}.tmp -o $outP/\${INDAT[0]}.bam

samtools index $outP/\${INDAT[0]}.bam

END_SH
	return $SHbwa;
}

1;
