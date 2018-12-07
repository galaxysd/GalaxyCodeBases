## Copyright Broad Institute, 2017
##
## This WDL pipeline implements data pre-processing and initial calling for somatic SNP,
## Indel, and copy number variants in human whole-genome sequencing (WGS) data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.
##
## For documentation on the M2 and CNV parameters, please see the respective WDL files (imported below).
##

import "SomaticPairedSingleSampleWf.wdl" as PreProcess
import "mutect2.wdl" as M2
import "cnv_somatic_pair_workflow.wdl" as cnvSomaticPairWorkflow

workflow FullSomaticPipeline {

    ### Preprocessing parameters
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    File wgs_coverage_interval_list

    String tumor_base_file_name
    Array[File] tumor_flowcell_unmapped_bams

    String normal_base_file_name
    Array[File] normal_flowcell_unmapped_bams
    String unmapped_bam_suffix

    Int read_length = 250

    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_alt
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac

    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices

    Int preemptible_tries
    Int agg_preemptible_tries

    Float cutoff_for_large_rg_in_gb = 20.0

    # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
    Int? increase_disk_size

    Int compression_level = 2
    #########################################


    #### M2 parameters
    File? pon
    File? pon_index
    Int scatter_count
    File? gnomad
    File? gnomad_index
    File? variants_for_contamination
    File? variants_for_contamination_index
    Boolean is_run_orientation_bias_filter = true
    Boolean is_run_oncotator = true

    File? onco_ds_tar_gz
    String? onco_ds_local_db_dir
    Array[String] artifact_modes
    String? m2_extra_args
    String? m2_extra_filtering_args
    String? sequencing_center
    String? sequence_source
    File? default_config_file

    Int? preemptible_attempts
    String basic_bash_docker = "ubuntu:16.04"
    String oncotator_docker = "broadinstitute/oncotator:1.9.6.1"

    #####################################

    ### CNV parameters
    #File intervals
    File common_sites
    File read_count_pon

    String gatk_docker
    File? gatk4_jar_override
    Int? bin_length

    Int? mem_gb_for_model_segments

    call PreProcess.SomaticPairedEndSingleSampleWorkflow as PreProcessTumor {
        input:
            contamination_sites_ud = contamination_sites_ud,
            contamination_sites_bed = contamination_sites_bed,
            contamination_sites_mu = contamination_sites_mu,
            wgs_coverage_interval_list = wgs_coverage_interval_list,

            base_file_name = tumor_base_file_name,
            flowcell_unmapped_bams = tumor_flowcell_unmapped_bams,
            unmapped_bam_suffix = unmapped_bam_suffix,

            read_length = read_length,

            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            ref_alt = ref_alt,
            ref_bwt = ref_bwt,
            ref_sa = ref_sa,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_pac = ref_pac,

            dbSNP_vcf = dbSNP_vcf,
            dbSNP_vcf_index = dbSNP_vcf_index,
            known_indels_sites_VCFs = known_indels_sites_VCFs,
            known_indels_sites_indices = known_indels_sites_indices,

            preemptible_tries = preemptible_tries,
            agg_preemptible_tries = agg_preemptible_tries,

            cutoff_for_large_rg_in_gb = cutoff_for_large_rg_in_gb,

            increase_disk_size = increase_disk_size,

            compression_level = compression_level
    }

    call PreProcess.SomaticPairedEndSingleSampleWorkflow as PreProcessNormal {
        input:
            contamination_sites_ud = contamination_sites_ud,
            contamination_sites_bed = contamination_sites_bed,
            contamination_sites_mu = contamination_sites_mu,
            wgs_coverage_interval_list = wgs_coverage_interval_list,

            base_file_name = normal_base_file_name,
            flowcell_unmapped_bams = normal_flowcell_unmapped_bams,
            unmapped_bam_suffix = unmapped_bam_suffix,

            read_length = read_length,

            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            ref_alt = ref_alt,
            ref_bwt = ref_bwt,
            ref_sa = ref_sa,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_pac = ref_pac,

            dbSNP_vcf = dbSNP_vcf,
            dbSNP_vcf_index = dbSNP_vcf_index,
            known_indels_sites_VCFs = known_indels_sites_VCFs,
            known_indels_sites_indices = known_indels_sites_indices,

            preemptible_tries = preemptible_tries,
            agg_preemptible_tries = agg_preemptible_tries,

            cutoff_for_large_rg_in_gb = cutoff_for_large_rg_in_gb,

            increase_disk_size = increase_disk_size,

            compression_level = compression_level
    }

    call M2.Mutect2 as M2Pair {
        input:
            intervals = wgs_coverage_interval_list,
            tumor_bam = PreProcessTumor.output_bam,
            tumor_bai = PreProcessTumor.output_bam_index,
            normal_bam = PreProcessNormal.output_bam,
            normal_bai = PreProcessNormal.output_bam_index,
            pon = pon,
            pon_index = pon_index,
            scatter_count = scatter_count,
            gnomad = gnomad,
            gnomad_index = gnomad_index,
            variants_for_contamination = variants_for_contamination,
            variants_for_contamination_index = variants_for_contamination_index,
            run_orientation_bias_filter = is_run_orientation_bias_filter,
            run_oncotator = is_run_oncotator,

            gatk_override = gatk4_jar_override,
            onco_ds_tar_gz = onco_ds_tar_gz,
            onco_ds_local_db_dir = onco_ds_local_db_dir,
            artifact_modes = artifact_modes,
            m2_extra_args = m2_extra_args,
            m2_extra_filtering_args = m2_extra_filtering_args,
            sequencing_center = sequencing_center,
            sequence_source = sequence_source,
            default_config_file = default_config_file,

            preemptible_attempts = preemptible_attempts,
            gatk_docker = gatk_docker,
            basic_bash_docker = basic_bash_docker,
            oncotator_docker = oncotator_docker,

            ref_fasta = ref_fasta,
            ref_fai = ref_fasta_index,
            ref_dict = ref_dict,

            emergency_extra_disk = 20
    }

    call cnvSomaticPairWorkflow.CNVSomaticPairWorkflow as CNVPair {
        input:
            intervals = wgs_coverage_interval_list,
            common_sites = common_sites,
            tumor_bam = PreProcessTumor.output_bam,
            tumor_bam_idx = PreProcessTumor.output_bam_index,
            normal_bam = PreProcessNormal.output_bam,
            normal_bam_idx = PreProcessNormal.output_bam_index,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_dict,
            ref_fasta_fai = ref_fasta_index,
            read_count_pon = read_count_pon,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            is_run_oncotator = is_run_oncotator,
            bin_length = bin_length,
            mem_gb_for_model_segments = mem_gb_for_model_segments
    }
}


