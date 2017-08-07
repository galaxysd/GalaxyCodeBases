* From <https://bioinfo.uth.edu/VirusFinder/VirusFinder-2.0.tgz>.
* Cleanup:

````bash
rm VirusFinder2.0/icorn/iCORN-v0.97.tar.gz
rm VirusFinder2.0/icorn/pileup_v0.5b/ssaha_pileup/ssaha_pileup/ssaha_*[r,n,l,s,e,p,a]	# can be make
chmod -x VirusFinder2.0/icorn/pileup_v0.5b/ssaha_pileup/ssaha_pileup/*.c
gzip -9 VirusFinder2.0/icorn/pileup_v0.5b/ssaha_pileup/ssaha_pileup/pileup-old.tar
rm VirusFinder2.0/icorn/pileup_v0.5b/ssaha_pileup/ssaha_pileup/pileup.tar	# duplicate file
rm -r VirusFinder2.0/icorn/pileup_v0.5b/ssaha2/ssaha2_v2.1.2_x86_64/	# See ssaha2_v2.1.2_x86_64.tgz
cd VirusFinder2.0/icorn/pileup_v0.5b/ssaha2/
	tar -czvf ssaha2_v2.3.tgz ssaha2-2.3_*
	rm ssaha2-2.3_*[4,6]
# Also repack ssaha2_v2.1.2.tgz (ssaha2_v2.1.2_i686.tgz ssaha2_v2.1.2_ia64.tgz ssaha2_v2.1.2_x86_64.tgz)
rm VirusFinder2.0/bin/GenomeAnalysisTK.jar	# GATK version 2.4-9-g532efad
````

The *ssaha2* used is `VirusFinder2.0/icorn/pileup_v0.5b/ssaha2/ssaha2-2.3_x86_64`

## Patch

* Use `samtools` <= 0.1.20.
* Use **absolute path** for `$output_dir` and `$config_file`.
* `$blastn_index_human.fa` will be used as reference file name, `ln -s hg38.fa hg38.fa.fa` if you didnot `formatdb -n`.
