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
#\$ -l vf=600M,num_proc=2,h=cngb-compute-f8-*
#\$ -binding linear:3
#\$ -cwd -r y
#\$ -v PERL5LIB,PATH,LD_LIBRARY_PATH
#\$ -e /dev/null -o /dev/null

#\$ -t 1-$cnt

cd $cwd
mkdir -p $outP

if [ -t 1 ] && [ "x\$SKIP" = 'x1' ]; then
	while read -u 9 THELINE; do
		read -ra INDAT <<<"\$THELINE"
		if [ -e "\${INDAT[1]}.fq.gz" ]; then
			ln -si \${INDAT[1]}.fq.gz $outP/\${INDAT[0]}.fq.gz
		elif [ -e "\${INDAT[1]}_1.fq.gz" ] && [ -e "\${INDAT[1]}_2.fq.gz" ]; then
			echo "[!]Cannot use SKIP mode for PE, using Read 1 only [\${INDAT[1]}_1.fq.gz] !"
			ln -si \${INDAT[1]}_1.fq.gz $outP/\${INDAT[0]}.fq.gz
		elif [ -e "\${INDAT[1]}_rawlib.basecaller.bam" ]; then
			echo "[!]PROTON/IonTorrent machine detected, assuming SE mode [\${INDAT[1]}_rawlib.basecaller.bam] !"
			ln -si \${INDAT[1]}_rawlib.basecaller.bam $outP/\${INDAT[0]}.ubam
		else
			echo "[x]Cannot find .fq.gz file(s) for [\${INDAT[1]}(_[12])?.fq.gz] !"
		fi
	done 9<$flst
	exit 0
fi

INFILE=`sed -n "\${SGE_TASK_ID}p" $flst`
read -ra INDAT <<<"\$INFILE"

SEADAPTERSTR='-g GAACGACATGGCTACGATCCGACTT -a AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA'
PEADAPTERSTR="\$SEADAPTERSTR -A AAGTCGGATCGTAGCCATGTCGTTC -G TTGTCTTCCTAAGACCGCTTGGCCTCCGACTT"
IONSEADAPTERSTR='-g GAACGACATGGCTACGATCCGACTT -a AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA'

if [ -e "\${INDAT[1]}.fq.gz" ]; then
	cutadapt -j 2 -m 1 -O 6 -n 2 \$SEADAPTERSTR -o $outP/\${INDAT[0]}.fq.gz \${INDAT[1]}.fq.gz >$outP/\${INDAT[0]}.log
elif [ -e "\${INDAT[1]}_1.fq.gz" ] && [ -e "\${INDAT[1]}_2.fq.gz" ]; then
	cutadapt -j 2 -m 1 -O 6 -n 2 \$PEADAPTERSTR --interleaved -o $outP/\${INDAT[0]}.fq.gz \${INDAT[1]}_1.fq.gz \${INDAT[1]}_2.fq.gz >$outP/\${INDAT[0]}.log
elif [ -e "\${INDAT[1]}_rawlib.basecaller.bam" ]; then
	samtools fastq \${INDAT[1]}_rawlib.basecaller.bam | cutadapt -j 2 -m 1 -O 6 -n 2 \$IONSEADAPTERSTR -o $outP/\${INDAT[0]}.fq.gz - >$outP/\${INDAT[0]}.log
else
	echo "[x]Cannot find .fq.gz file(s) for [\${INDAT[1]}(_[12])?.fq.gz] !"
	exit 1
fi

if [ ! -s $outP/\${INDAT[0]}.fq.gz ]; then
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
mkdir -p $outP

INFILE=`sed -n "\${SGE_TASK_ID}p" $flst`
read -ra INDAT <<<"\$INFILE"

if [ -e "$FQprefix/\${INDAT[0]}.fq.gz" ]; then
	bwa mem -t 12 -Y $pRef -R "\@RG\\tID:\${INDAT[0]}\\tSM:\${INDAT[0]}" -p $FQprefix/\${INDAT[0]}.fq.gz 2>$outP/\${INDAT[0]}.log | samtools view -bS - | samtools sort -m 2G -T $outP/\${INDAT[0]}.tmp -o $outP/\${INDAT[0]}.bam
elif [ -e "$FQprefix/\${INDAT[0]}.ubam" ]; then
	samtools fastq $FQprefix/\${INDAT[0]}.ubam | bwa mem -t 12 -Y $pRef -R "\@RG\\tID:\${INDAT[0]}\\tSM:\${INDAT[0]}" - 2>$outP/\${INDAT[0]}.log | samtools view -bS - | samtools sort -m 2G -T $outP/\${INDAT[0]}.tmp -o $outP/\${INDAT[0]}.bam
fi

samtools index $outP/\${INDAT[0]}.bam

if [ ! -s $outP/\${INDAT[0]}.bam ]; then
	exit 1
fi
END_SH
	return $SHbwa;
}

sub Smpileup($$$$$$$$$) {
	my ($cwd,$cnt,$flst,$outP,$BAMprefix,$LSTprefix,$pRef,$OYKprefix,$theMode) = @_;
	my $Smpileup = <<"END_SH";
#!/bin/bash
#\$ -S /bin/bash
#\$ -q $SgeQueue -P $SgeProject
#\$ -N p2${SgeJobPrefix}mplp -hold_jid p1${SgeJobPrefix}bwa
#\$ -l vf=500M,num_proc=1,h=cngb-compute-f8-*
#\$ -binding linear:2
#\$ -cwd -r y
#\$ -v PERL5LIB,PATH,LD_LIBRARY_PATH
#\$ -e /dev/null -o /dev/null

#\$ -t 1-$cnt

cd $cwd
mkdir -p $outP

INFILE=`sed -n "\${SGE_TASK_ID}p" $flst`
read -ra INDAT <<<"\$INFILE"

samtools mpileup -b $LSTprefix/p\${INDAT[2]}.bams.lst -d 4000 -Q 20 -f $pRef -v -t 'DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR' -p -o $outP/\${INDAT[2]}.vcf.gz
bcftools call -Oz -v -m $outP/\${INDAT[2]}.vcf.gz -o $outP/\${INDAT[2]}.snp.gz
bcftools index $outP/\${INDAT[2]}.vcf.gz &
bcftools index $outP/\${INDAT[2]}.snp.gz


bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\t%QUAL[\\t%TGT;%AD]\\n' -S $LSTprefix/p\${INDAT[2]}.M.lst -i'POS=501' $outP/\${INDAT[2]}.snp.gz >$OYKprefix/p\${INDAT[2]}.M.tsv
bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\t%QUAL[\\t%TGT;%AD]\\n' -S $LSTprefix/p\${INDAT[2]}.F.lst -i'POS=501' $outP/\${INDAT[2]}.snp.gz >$OYKprefix/p\${INDAT[2]}.F.tsv
bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\t%QUAL[\\t%TGT;%AD]\\n' -S $LSTprefix/p\${INDAT[2]}.C.lst -i'POS=501' $outP/\${INDAT[2]}.snp.gz >$OYKprefix/p\${INDAT[2]}.C.tsv

./oyka.pl $theMode $OYKprefix/p\${INDAT[2]}.M.tsv $OYKprefix/p\${INDAT[2]}.F.tsv $OYKprefix/p\${INDAT[2]}.C.tsv $OYKprefix/r\${INDAT[2]}

if [ ! -s $outP/\${INDAT[2]}.snp.gz ]; then
	exit 1
fi
END_SH
	return $Smpileup;
}


1;
