# Usage

```
./virusinserts.pl <Host.fa[.gz]> <Virus.fa[.gz]> <Outprefix>

./virusinserts.pl /share/human/Chr1.fa /share/bsvir/HBV.AJ507799.2.fa test
```

This will read "Ref.fa" as **Host genome**, and "Virus.fa" as **Virus sequence**, save output file to "Outprefix.ref.fa","Outprefix.1.fq","Outprefix.2.fq".

For both input files, **only** the first sequence is loaded. All follow sequences, if any, are simply ignored.

# Description

* It will random pick `$SampleCnt = 100` fragments from **Host**, each with a length of twice of `$PEinsertLen = 200`.
* Then, `$SampleCnt = 100` samples from **Virus** with a length of [`$VirFragMin = 20` , `$VirFragMax = 500`] are selected.
* For each of the `$SampleCnt = 100` sets, *Virus fragment* is insereted into the middle of *host fragment* to form our **Simulation Reference**, which will be write to "Outprefix.ref.fa" in FASTA format. The direction of *Virus fragment* is random chosen.
* On each *Simulation Reference*, PE90 of inster size 200bp are picked follow uniform distribution. Reads are write to "Outprefix.1.fq" and "Outprefix.2.fq" in FASTQ format.
