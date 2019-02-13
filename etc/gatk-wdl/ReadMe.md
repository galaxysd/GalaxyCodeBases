# My GATK4 WorkFlows

## seq-format-conversion
Workflows for converting between sequence data formats

<https://github.com/gatk-workflows/seq-format-conversion>

## Local CMD

```bash
time cromwell run paired-fastq-to-unmapped-bam.wdl -i paired-fastq-to-unmapped-bam.inputs.json >cromwell.uBAM.log &
```

    You can also use `date '+%Y%m%d%H%M%S'` for unique strings.

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



## gatk4-data-processing

<https://github.com/gatk-workflows/gatk4-data-processing>

## Local CMD

```bash
time cromwell run processing-for-variant-discovery-gatk4.wdl -i processing-for-variant-discovery-gatk4.hg38.wgs.D3B.inputs.json > cromwell.processing.D3B.log &
time cromwell run processing-for-variant-discovery-gatk4.wdl -i processing-for-variant-discovery-gatk4.hg38.wgs.Normal.inputs.json > cromwell.processing.Normal.log &
```

### Purpose :
Workflows for processing high-throughput sequencing data for variant discovery with GATK4 and related tools.

### processing-for-variant-discovery-gatk4 :
The processing-for-variant-discovery-gatk4 WDL pipeline implements data pre-processing according to the GATK Best Practices
(June 2016).

#### Requirements/expectations
- Pair-end sequencing data in unmapped BAM (uBAM) format
- One or more read groups, one per uBAM file, all belonging to a single sample (SM)
- Input uBAM files must additionally comply with the following requirements:
  - filenames all have the same suffix (we use ".unmapped.bam")
  - files must pass validation by ValidateSamFile
  - reads are provided in query-sorted order
  - all reads must have an RG tag

#### Outputs
- A clean BAM file and its index, suitable for variant discovery analyses.

### Software version requirements :
- GATK 4 or later
- Picard 2.x
- Samtools (see gotc docker)
- Python 2.7



## somatic-snvs-indels

<https://github.com/gatk-workflows/gatk4-somatic-snvs-indels>

### Purpose :
Workflows for somatic short variant analysis with GATK4.

### mutect2 :
Implements Somatic short variant discovery using [GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/workflow).
Note: Also provided in this repo is mutect2_nio which is a NIO supported version of the wdl.

#### Requirements/expectations
- Tumor bam and index
- Normal bam and index

#### Outputs
- unfiltered vcf
- unfiltered vcf index
- filtered vcf
- filtered vcf index

### mutect2_pon :
Creates a Panel of Norms to be implemented in somatic short variant discovery.

#### Requirements/expectations
- Normal bams and index

#### Outputs
- PON vcf and index
- Normal calls vcf and index

### mutect2-normal-normal :
Used to validate mutect2 workflow.

#### Requirements/expectations
- One analysis-ready BAM file (and its index) for each replicate

#### Outputs
- False Positive VCF files and its index with summary

### Software version requirements :
- GATK4 or later

Cromwell version support
- Successfully tested on v31


### Parameter descriptions :
#### mutect2 (single pair/sample)
- ``Mutect2.gatk4_jar`` -- Location *within the docker file* of the GATK4 jar file.  If you wish you to use a different jar file, such as one on your local filesystem or a google bucket, specify that location with ``Mutect2_Multi.gatk4_jar_override``.  This parameter is ignored if ``Mutect2_Multi.gatk4_jar_override`` is specified.
- ``Mutect2.intervals`` -- A file listing genomic intervals to search for somatic mutations.  This should be in the standard GATK4 format.
- ``Mutect2.ref_fasta`` 	-- reference fasta.  For Broad internal VM:  ``/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta``
- ``Mutect2.ref_fasta_index`` -- For Broad internal VM:  ``/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta.fai``
- ``Mutect2.ref_dict`` -- For Broad internal VM:  ``/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.dict``
- ``Mutect2.tumor_bam`` -- File path or storage location (depending on backend) of the tumor bam file.
- ``Mutect2.tumor_bam_index`` -- File path or storage location (depending on backend) of the tumor bam file index.
- ``Mutect2.normal_bam`` -- (optional) File path or storage location (depending on backend) of the normal bam file.
- ``Mutect2.normal_bam_index`` -- (optional, but required if ``Mutect2.normal_bam`` is specified)  File path or storage location (depending on backend) of the normal bam file index.
- ``Mutect2.pon`` -- (optional) Panel of normals VCF to use for false positive reduction.
- ``Mutect2.pon_index`` -- (optional, but required if ``Mutect2_Multi.pon`` is specified)  VCF index for the panel of normals.  Please see GATK4 tool ``IndexFeatureFile`` for creation of an index.
- ``Mutect2.scatter_count`` -- Number of executions to split the Mutect2 task into.  The more you put here, the faster Mutect2 will return results, but at a higher cost of resources.
- ``Mutect2.gnomad`` -- (optional)  gnomAD vcf containing population allele frequencies (AF) of common and rare alleles.  Download an exome or genome sites vcf [here](http://gnomad.broadinstitute.org/downloads).  Essential for determining possible germline variants in tumor
- ``Mutect2.gnomad_index`` -- (optional, but required if ``Mutect2_Multi.gnomad`` is specified)  VCF index for gnomAD.  Please see GATK4 tool ``IndexFeatureFile`` for creation of an index.
- ``Mutect2.variants_for_contamination`` -- (optional)  vcf containing population allele frequencies (AF) of common SNPs.  If omitted, cross-sample contamination will not be calculated and contamination filtering will not be applied.  This can be generated from a gnomAD vcf using the GATK4 tool ``SelectVariants`` with the argument ``--select "AF > 0.05"``.  For speed, one can get very good results using only SNPs on chromosome 1.  For example, ``java -jar $gatk SelectVariants -V gnomad.vcf -L 1 --select "AF > 0.05" -O variants_for_contamination.vcf``.
- ``Mutect2.variants_for_contamination_index`` -- (optional, but required if ``Mutect2_Multi.variants_for_contamination`` is specified)  VCF index for contamination variants.  Please see GATK4 tool ``IndexFeatureFile`` for creation of an index.
- ``Mutect2.is_run_orientation_bias_filter`` -- ``true``/``false`` whether the orientation bias filter should be run.
- ``Mutect2.is_run_oncotator`` -- ``true``/``false`` whether the command-line version of oncotator should be run.  If ``false``, ``Mutect2_Multi.oncotator_docker`` parameter is ignored.
- ``Mutect2.gatk_docker`` -- Docker image to use for Mutect2 tasks.  This is only used for backends configured to use docker.
- ``Mutect2.oncotator_docker`` -- (optional)  A GATK4 jar file to be used instead of the jar file in the docker image.  (See ``Mutect2_Multi.gatk4_jar``)  This can be very useful for developers.  Please note that you need to be careful that the docker image you use is compatible with the GATK4 jar file given here -- no automated checks are made.
- ``Mutect2.gatk4_jar_override`` -- (optional)  A GATK4 jar file to be used instead of the jar file in the docker image.  (See ``Mutect2_Multi.gatk4_jar``)  This can be very useful for developers.  Please note that you need to be careful that the docker image you use is compatible with the GATK4 jar file given here
- ``Mutect2.preemptible_attempts`` -- Number of times to attempt running a task on a preemptible VM.  This is only used for cloud backends in cromwell and is ignored for local and SGE backends.
- ``Mutect2.onco_ds_tar_gz`` -- (optional)  A tar.gz file of the oncotator datasources -- often quite large (>15GB).  This will be uncompressed as part of the oncotator task.  Depending on backend used, this can be specified as a path on the local filesystem of a cloud storage container (e.g. gs://...).  Typically the Oncotator default datasource can be downloaded at ``ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/oncotator/``.  Do not put the FTP URL into the json file.
- ``Mutect2.onco_ds_local_db_dir`` -- "(optional)  A direct path to the Oncotator datasource directory (uncompressed).  While this is the fastest approach, it cannot be used with docker unless your docker image already has the datasources in it.  For cromwell backends without docker, this can be a local filesystem path.  *This cannot be a cloud storage location*

 Note:  If neither ``Mutect2_Multi.onco_ds_tar_gz``, nor ``Mutect2_Multi.onco_ds_local_db_dir``, is specified, the Oncotator task will download and uncompress for each execution.

The following three parameters are useful for rendering TCGA MAFs using oncotator.  These parameters are ignored if ``is_run_oncotator`` is ``false``."
- ``Mutect2.artifact_modes`` -- List of artifact modes to search for in the orientation bias filter.  For example to filter the OxoG artifact, you would specify ``["G/T"]``.  For both the FFPE artifact and the OxoG artifact, specify ``["G/T", "C/T"]``.  If you do not wish to search for any artifacts, please set ``Mutect2_Multi.is_run_orientation_bias_filter`` to ``false``.
- ``Mutect2.picard_jar`` -- A direct path to a picard jar for using ``CollectSequencingArtifactMetrics``.  This parameter requirement will be eliminated in the future.
- ``Mutect2.m2_extra_args`` -- (optional) a string of additional command line arguments of the form "-argument1 value1 -argument2 value2" for Mutect 2.  Most users will not need this.
- ``Mutect2.m2_extra_filtering_args`` -- (optional) a string of additional command line arguments of the form "-argument1 value1 -argument2 value2" for Mutect 2.  Most users will not need this.
- ``Mutect2.sequencing_center`` 	-- (optional) center reporting this variant.     Please see ``https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4`` for more details.
- ``Mutect2.sequence_source`` -- (optional)  ``WGS`` or ``WXS`` for whole genome or whole exome sequencing, respectively.  Please note that the controlled vocabulary of the TCGA MAF spec is *not* enforced.  Please see ``https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4`` for more details.
- ``Mutect2.default_config_file`` -- "(optional)  A configuration file that can direct oncotator to use default values for unspecified annotations in the TCGA MAF.  This help prevents having MAF files with a lot of ""__UNKNOWN__"" values.  An usable example is given below.  Here is an example that should work for most users:

```
[manual_annotations]
override:NCBI_Build=37,Strand=+,status=Somatic,phase=Phase_I,sequencer=Illumina,Tumor_Validation_Allele1=,Tumor_Validation_Allele2=,Match_Norm_Validation_Allele1=,Match_Norm_Validation_Allele2=,Verification_Status=,Validation_Status=,Validation_Method=,Score=,BAM_file=,Match_Norm_Seq_Allele1=,Match_Norm_Seq_Allele2=
```
- ``Mutect2.filter_oncotator_maf`` -- (optional, default true) Whether Oncotator should remove filtered variants when rendering the MAF.  Ignored if `run_oncotator` is false.



## gatk4-germline-snps-indels

<https://github.com/gatk-workflows/gatk4-germline-snps-indels>

### Purpose : 
Workflows for germline short variant discovery with GATK4. 

### haplotypecaller-gvcf-gatk :
The haplotypecaller-gvcf-gatk4 workflow runs HaplotypeCaller 
from GATK4 in GVCF mode on a single sample according to the GATK Best Practices (June 2016), 
scattered across intervals.

#### Requirements/expectations
- One analysis-ready BAM file for a single sample (as identified in RG:SM)
- Set of variant calling intervals lists for the scatter, provided in a file
#### Outputs 
- One GVCF file and its index

### joint-discovery-gatk :
The second WDL implements the joint discovery and VQSR 
filtering portion of the GATK Best Practices (June 2016) for germline SNP and Indel 
discovery in human whole-genome sequencing (WGS) and exome sequencing data.

*NOTE: joint-discovery-gatk4-local.wdl is a slightly modified version of the original to support users interested in running the workflow locally.*

#### Requirements/expectations
- One or more GVCFs produced by HaplotypeCaller in GVCF mode
- Bare minimum 1 WGS sample or 30 Exome samples. Gene panels are not supported.
- When deteriming disk size in the json, use the guideline below
  - small_disk = (num_gvcfs / 10) + 10
  - medium_disk = (num_gvcfs * 15) + 10
  - huge_disk = num_gvcfs + 10

### Outputs 
- A VCF file and its index, filtered using variant quality score recalibration  
  (VQSR) with genotypes for all samples present in the input VCF. All sites that  
  are present in the input VCF are retained; filtered sites are annotated as such  
  in the FILTER field.

### Software version requirements :
- GATK 4 or later 
- Samtools (see gotc docker)
- Python 2.7

Cromwell version support 
- Successfully tested on v31
- Does not work on versions < v23 due to output syntax



---



## gatk-somatic-with-preprocessing

This WDL pipeline implements data pre-processing and initial calling for somatic SNP,
Indel, and copy number variants in human whole-genome sequencing (WGS) data.

<https://github.com/gatk-workflows/gatk4-somatic-with-preprocessing>

Note: The gatk-somatic-with-preprocessing WDL is not used in any pipelines at the Broad Institute
and has been provided only as a convenience for the community.  Therefore, this WDL is unsupported.

## Local CMD

```bash
time cromwell run FullSomaticPipeline.wdl --imports FullSomaticPipeline.imports.zip -i FullSomaticPipeline.json >cromwell.02.log &
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

---

## somatic-cnvs
Workflows for somatic copy number variant analysis

<https://github.com/gatk-workflows/gatk4-somatic-cnvs>

## Running the Somatic CNV WDL

### Which WDL should you use?

- Building a panel of normals (PoN): ``cnv_somatic_panel_workflow.wdl``
- Running a matched pair: ``cnv_somatic_pair_workflow.wdl``

#### Setting up parameter json file for a run

To get started, create the json template (using ``java -jar wdltool.jar inputs <workflow>``) for the workflow you wish to run and adjust parameters accordingly.

*Please note that there are optional workflow-level and task-level parameters that do not appear in the template file.  These are set to reasonable values by default, but can also be adjusted if desired.*

#### Required parameters in the somatic panel workflow

Important: The normal_bams samples in the json can be used test the wdl, they are NOT to be used to create a panel of normals for sequence analysis. For instructions on creating a proper PON please refer to user the documents https://software.broadinstitute.org/gatk/documentation/ .

The reference used must be the same between PoN and case samples.

- ``CNVSomaticPanelWorkflow.gatk_docker`` -- GATK Docker image (e.g., ``broadinstitute/gatk:latest``).
- ``CNVSomaticPanelWorkflow.intervals`` -- Picard or GATK-style interval list.  For WGS, this should typically only include the autosomal chromosomes.
- ``CNVSomaticPanelWorkflow.normal_bais`` -- List of BAI files.  This list must correspond to `normal_bams`.  For example, `["Sample1.bai", "Sample2.bai"]`.
- ``CNVSomaticPanelWorkflow.normal_bams`` -- List of BAM files.  This list must correspond to `normal_bais`.  For example, `["Sample1.bam", "Sample2.bam"]`.
- ``CNVSomaticPanelWorkflow.pon_entity_id`` -- Name of the final PoN file.
- ``CNVSomaticPanelWorkflow.ref_fasta_dict`` -- Path to reference dict file.
- ``CNVSomaticPanelWorkflow.ref_fasta_fai`` -- Path to reference fasta fai file.
- ``CNVSomaticPanelWorkflow.ref_fasta`` -- Path to reference fasta file.

In additional, there are optional workflow-level and task-level parameters that may be set by advanced users; for example:

- ``CNVSomaticPanelWorkflow.do_explicit_gc_correction`` -- (optional) If true, perform explicit GC-bias correction when creating PoN and in subsequent denoising of case samples.  If false, rely on PCA-based denoising to correct for GC bias.
- ``CNVSomaticPanelWorkflow.PreprocessIntervals.bin_length`` -- Size of bins (in bp) for coverage collection.  *This must be the same value used for all case samples.*
- ``CNVSomaticPanelWorkflow.PreprocessIntervals.padding`` -- Amount of padding (in bp) to add to both sides of targets for WES coverage collection.  *This must be the same value used for all case samples.*

Further explanation of other task-level parameters may be found by invoking the ``--help`` documentation available in the gatk.jar for each tool.

#### Required parameters in the somatic pair workflow

The reference and bins (if specified) must be the same between PoN and case samples.

- ``CNVSomaticPairWorkflow.common_sites`` -- Picard or GATK-style interval list of common sites to use for collecting allelic counts.
- ``CNVSomaticPairWorkflow.gatk_docker`` -- GATK Docker image (e.g., ``broadinstitute/gatk:latest``).
- ``CNVSomaticPairWorkflow.intervals`` -- Picard or GATK-style interval list.  For WGS, this should typically only include the autosomal chromosomes.
- ``CNVSomaticPairWorkflow.normal_bam`` -- Path to normal BAM file.
- ``CNVSomaticPairWorkflow.normal_bam_idx`` -- Path to normal BAM file index.
- ``CNVSomaticPairWorkflow.read_count_pon`` -- Path to read-count PoN created by the panel workflow. 
- ``CNVSomaticPairWorkflow.ref_fasta_dict`` -- Path to reference dict file.
- ``CNVSomaticPairWorkflow.ref_fasta_fai`` -- Path to reference fasta fai file.
- ``CNVSomaticPairWorkflow.ref_fasta`` -- Path to reference fasta file.
- ``CNVSomaticPairWorkflow.tumor_bam`` -- Path to tumor BAM file.
- ``CNVSomaticPairWorkflow.tumor_bam_idx`` -- Path to tumor BAM file index.

In additional, there are several task-level parameters that may be set by advanced users as above.

To invoke Oncotator on the called tumor copy-ratio segments:

- ``CNVSomaticPairWorkflow.is_run_oncotator`` -- (optional) If true, run Oncotator on the called copy-ratio segments.  This will generate both a simple TSV and a gene list.

---

