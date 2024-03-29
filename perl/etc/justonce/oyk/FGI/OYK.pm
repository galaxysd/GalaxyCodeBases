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
IONSEADAPTERSTR='-g CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT -a ATCnnnnnnnnnnCTGAGTCGGAGACACGCAGGGATGAGATGGTT'
IONPEADAPTERSTR="\$IONSEADAPTERSTR -G CCATCTCATCCCTGCGTGTCTCCGACTCAGnnnnnnnnnnGAT -A ATCACCGACTGCCCATAGAGAGGAAAGCGGAGGCGTAGTGGTT"
HISEQSEADAPTERSTR='-g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACnnnnnnnnATCTCGTATGCCGTCTTCTGCTTG'
HISEQPEADAPTERSTR="\$HISEQSEADAPTERSTR -G GATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -G CAAGCAGAAGACGGCATACGAGATnnnnnnnnGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"

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

sub Sbwamem($$$$$) {
	my ($cwd,$cnt,$flst,$outP,$FQprefix) = @_;
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
	bwa mem -t 12 -Y $RealBin/$rRefn -R "\@RG\\tID:\${INDAT[0]}\\tSM:\${INDAT[0]}" -p $FQprefix/\${INDAT[0]}.fq.gz 2>$outP/\${INDAT[0]}.log | samtools view -bS - | samtools sort -m 2G -T $outP/\${INDAT[0]}.tmp -o $outP/\${INDAT[0]}.bam
elif [ -e "$FQprefix/\${INDAT[0]}.ubam" ]; then
	samtools fastq $FQprefix/\${INDAT[0]}.ubam | bwa mem -t 12 -Y $RealBin/$rRefn -R "\@RG\\tID:\${INDAT[0]}\\tSM:\${INDAT[0]}" - 2>$outP/\${INDAT[0]}.log | samtools view -bS - | samtools sort -m 2G -T $outP/\${INDAT[0]}.tmp -o $outP/\${INDAT[0]}.bam
fi

samtools index $outP/\${INDAT[0]}.bam
samtools stats -@ 12 $outP/\${INDAT[0]}.bam >$outP/\${INDAT[0]}.bam.stat
grep ^IS $outP/\${INDAT[0]}.bam.stat | cut -f 2- > $outP/\${INDAT[0]}.fragstats

if [ ! -s $outP/\${INDAT[0]}.bam ]; then
	exit 1
fi
END_SH
	return $SHbwa;
}

sub Smpileup($$$$$$$$$) {
	my ($cwd,$cnt,$flst,$outP,$BAMprefix,$LSTprefix,$OYKprefix,$theMode,$theParentage) = @_;
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

samtools mpileup -b $LSTprefix/p\${INDAT[2]}.bams.lst -d 4000 -Q 30 -f $RealBin/$rRefn -v -t 'DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR' -p -o $outP/\${INDAT[2]}.vcf.gz
bcftools call -Oz -v -m $outP/\${INDAT[2]}.vcf.gz -o $outP/\${INDAT[2]}.snp.gz
bcftools index $outP/\${INDAT[2]}.vcf.gz &
bcftools index $outP/\${INDAT[2]}.snp.gz


bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\t%QUAL[\\t%TGT;%AD]\\n' -S $LSTprefix/p\${INDAT[2]}.M.lst -i'POS=501' $outP/\${INDAT[2]}.snp.gz >$OYKprefix/p\${INDAT[2]}.M.tsv
bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\t%QUAL[\\t%TGT;%AD]\\n' -S $LSTprefix/p\${INDAT[2]}.F.lst -i'POS=501' $outP/\${INDAT[2]}.snp.gz >$OYKprefix/p\${INDAT[2]}.F.tsv
bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\t%QUAL[\\t%TGT;%AD]\\n' -S $LSTprefix/p\${INDAT[2]}.C.lst -i'POS=501' $outP/\${INDAT[2]}.snp.gz >$OYKprefix/p\${INDAT[2]}.C.tsv

$RealBin/bin/oykn.pl $theMode $theParentage $OYKprefix/p\${INDAT[2]}.M.tsv $OYKprefix/p\${INDAT[2]}.F.tsv $OYKprefix/p\${INDAT[2]}.C.tsv $OYKprefix/r\${INDAT[2]}

$RealBin/bin/get_ChrNum.pl $OYKprefix/r\${INDAT[2]}.cpie $RealBin/db/nippt7274.tsv \${INDAT[1]} > $OYKprefix/r\${INDAT[2]}.F.txt

if [ ! -s $outP/\${INDAT[2]}.snp.gz ]; then
	exit 1
fi
END_SH
	return $Smpileup;
}

sub Sqc($$$$$$$$) {
	my ($cwd,$outP,$theMode,$theParentage,$theSam,$theFam,$theMachine,$theStore) = @_;
	my $Sqc = <<"END_SH";
#!/bin/bash
#\$ -S /bin/bash
#\$ -q $SgeQueue -P $SgeProject
#\$ -N p3${SgeJobPrefix}qc -hold_jid p2${SgeJobPrefix}mplp
#\$ -l vf=2G,num_proc=1,h=cngb-compute-f8-*
#\$ -binding linear:1
#\$ -cwd -r y
#\$ -v PERL5LIB,PATH,LD_LIBRARY_PATH
#\$ -e /dev/null -o /dev/null

cd $cwd
mkdir -p $outP
mkdir -p $outP/../6record

perl $RealBin/bin/contamination.F.pl $theMode $theParentage $outP $outP/contamination.F.txt
perl $RealBin/bin/contamination.FM.pl $outP $outP/contamination.FM.txt
perl $RealBin/bin/contamination.MC.pl $outP/../family.lst $outP/../4tsv $outP/contamination.MC.txt
perl $RealBin/bin/heterozygosity.pl $outP/../family.lst $outP/../4tsv $outP/heterozygosity.txt

perl $RealBin/bin/prepare.distribution.pl $outP/../family.lst $outP/../4tsv $outP/cffDNA.txt $outP/depth.txt
Rscript $RealBin/bin/ggplot_bar.cffdna.R $outP/cffDNA.txt $outP/cffDNA.pdf

if [[ "$theMode" == "CHIP" ]]; then
ls $outP/../2bam/*.fragstats > $outP/fragment.list
Rscript $RealBin/bin/ggplot_bar.frag.ReadList.R $outP/fragment.list $outP/fragment.pdf
	Rscript $RealBin/bin/ggplot_bar.depth.R $outP/depth.txt $outP/depth.pdf
elif [[ "$theMode" == "PCR" ]]; then
	Rscript $RealBin/bin/ggplot_bar.depth.forPCR.R $outP/depth.txt $outP/depth.pdf
fi

perl $RealBin/bin/base_qc.pl -sl $theSam -fl $theFam -mode $theMachine -s $theStore -r $outP/.. -o $outP
perl $RealBin/bin/record.pl $outP/../family.lst $outP/../4tsv $outP/../6record
END_SH
	return $Sqc;
}

1;
