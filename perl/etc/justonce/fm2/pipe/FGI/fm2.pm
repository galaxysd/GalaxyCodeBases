package main;
use strict;
use warnings;

my $SgeQueue = 'fgi.q';
my $SgeProject = 'fgi';
my $SgeJobPrefix = 'nip'.(time() % 9999);
our $rRefn;

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
IONSEADAPTERSTR='-g CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT -a ATCnnnnnnnnnnCTGAGTCGGAGACACGCAGGGATGAGATGGTT'
IONPEADAPTERSTR="\$IONSEADAPTERSTR -G CCATCTCATCCCTGCGTGTCTCCGACTCAGnnnnnnnnnnGAT -A ATCACCGACTGCCCATAGAGAGGAAAGCGGAGGCGTAGTGGTT"
HISEQSEADAPTERSTR='-g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACnnnnnnnnATCTCGTATGCCGTCTTCTGCTTG'
HISEQPEADAPTERSTR="\$HISEQSEADAPTERSTR -G GATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -G CAAGCAGAAGACGGCATACGAGATnnnnnnnnGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"

if [ -e "\${INDAT[1]}.fq.gz" ]; then
	cutadapt -j 2 -m 1 -O 6 -n 2 \$HISEQSEADAPTERSTR -o $outP/\${INDAT[0]}.fq.gz \${INDAT[1]}.fq.gz >$outP/\${INDAT[0]}.log
elif [ -e "\${INDAT[1]}_1.fq.gz" ] && [ -e "\${INDAT[1]}_2.fq.gz" ]; then
	cutadapt -j 2 -m 1 -O 6 -n 2 \$HISEQPEADAPTERSTR --interleaved -o $outP/\${INDAT[0]}.fq.gz \${INDAT[1]}_1.fq.gz \${INDAT[1]}_2.fq.gz >$outP/\${INDAT[0]}.log
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

sub Sbwamem($$$$$$$) {
	my ($cwd,$cnt,$flst,$outP,$VCFprefix,$FQprefix,$OYKprefix) = @_;
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
	bwa mem -t 12 -Y -O 12,12 $RealBin/$rRefn -R "\@RG\\tID:\${INDAT[0]}\\tSM:\${INDAT[0]}" -p $FQprefix/\${INDAT[0]}.fq.gz 2>$outP/\${INDAT[0]}.log | samtools view -bS - | samtools sort -m 2G -T $outP/\${INDAT[0]}.tmp -o $outP/\${INDAT[0]}.bam
elif [ -e "$FQprefix/\${INDAT[0]}.ubam" ]; then
	samtools fastq $FQprefix/\${INDAT[0]}.ubam | bwa mem -t 12 -Y -O 12,12 $RealBin/$rRefn -R "\@RG\\tID:\${INDAT[0]}\\tSM:\${INDAT[0]}" - 2>$outP/\${INDAT[0]}.log | samtools view -bS - | samtools sort -m 2G -T $outP/\${INDAT[0]}.tmp -o $outP/\${INDAT[0]}.bam
fi

samtools index $outP/\${INDAT[0]}.bam

java -jar $RealBin/bin/picard.jar FilterSamReads I=$outP/\${INDAT[0]}.bam O=$VCFprefix/\${INDAT[0]}_M.bam JAVASCRIPT_FILE=$RealBin/bin/f200M.js FILTER=includeJavascript
samtools mpileup -l $RealBin/bin/LN-mid.bed $VCFprefix/\${INDAT[0]}_M.bam |awk '{print \$1"\\t"\$4}' > $VCFprefix/\${INDAT[0]}.0str
$RealBin/bin/fstr.pl $VCFprefix/\${INDAT[0]}.0str >$OYKprefix/\${INDAT[0]}.str

bcftools mpileup --threads 6 $outP/\${INDAT[0]}.bam -d 30000 -Q 30 -f $RealBin/$rRefn -p -Ob -o $VCFprefix/\${INDAT[0]}.bcf
bcftools call -Oz -A -m $VCFprefix/\${INDAT[0]}.bcf -o $VCFprefix/\${INDAT[0]}.snp.gz
bcftools index $VCFprefix/\${INDAT[0]}.snp.gz
bcftools query -f'%CHROM\\t[%DP\\t%QUAL\\t%TGT\\n]' -i 'POS==501' $VCFprefix/\${INDAT[0]}.snp.gz >$VCFprefix/\${INDAT[0]}.0snp
$RealBin/bin/fsnp.pl $VCFprefix/\${INDAT[0]}.0snp >$OYKprefix/\${INDAT[0]}.snp

cat $OYKprefix/\${INDAT[0]}.snp $OYKprefix/\${INDAT[0]}.str >$OYKprefix/\${INDAT[0]}.result

samtools stats -@ 12 $outP/\${INDAT[0]}.bam >$outP/\${INDAT[0]}.bam.stat
grep ^IS $outP/\${INDAT[0]}.bam.stat | cut -f 2- > $outP/\${INDAT[0]}.fragstats

if [ ! -s $outP/\${INDAT[0]}.bam ]; then
	exit 1
fi
END_SH
	return $SHbwa;
}



1;
