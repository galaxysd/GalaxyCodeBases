# My GATK4 WorkFlows

## seq-format-conversion
Workflows for converting between sequence data formats

<https://github.com/gatk-workflows/seq-format-conversion>

## Local CMD

```bash
time cromwell run paired-fastq-to-unmapped-bam.wdl -i paired-fastq-to-unmapped-bam.inputs.json >cromwell.log 2>cromwell.err &
```

### paired-fastq-to-unmapped-bam :
This WDL converts paired FASTQ to uBAM and adds read group information

*NOTE: paired-fastq-to-unmapped-bam-fc.wdl is a slightly modified version of the original to support users interested running on FireCloud.
As input this wdl takes a TSV with each row being a different readgroup and each column in the row being descriptors*

#### Requirements/expectations
- Pair-end sequencing data in FASTQ format (one file per orientation)
- The following metada descriptors per sample:
```
readgroup   fastq_pair1_file_path   fastq_pair2_file_path   sample_name   library_name   platform_unit   run_date   platform_name   sequecing_center
```

#### Outputs
- Set of unmapped BAMs, one per read group
- File containing a list of the generated unmapped BAMs

## gatk-somatic-with-preprocessing

This WDL pipeline implements data pre-processing and initial calling for somatic SNP,
Indel, and copy number variants in human whole-genome sequencing (WGS) data.

<https://github.com/gatk-workflows/gatk4-somatic-with-preprocessing>

Note: The gatk-somatic-with-preprocessing WDL is not used in any pipelines at the Broad Institute
and has been provided only as a convenience for the community.  Therefore, this WDL is unsupported.

## Local CMD

```bash
time cromwell run FullSomaticPipeline.wdl --imports FullSomaticPipeline.imports.zip -i FullSomaticPipeline.json >cromwell.`date '+%Y%m%d%H%M%S'`.log 2>cromwell.`date '+%Y%m%d%H%M%S'`.err &
```

### Requirements/expectations
 - Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
 - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
 - Input uBAM files must additionally comply with the following requirements:
 - - filenames all have the same suffix (we use ".unmapped.bam")
 - - files must pass validation by ValidateSamFile
 - - reads are provided in query-sorted order
 - - all reads must have an RG tag
 - Reference genome must be Hg38 with ALT contigs

