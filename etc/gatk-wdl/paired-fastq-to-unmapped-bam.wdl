##Copyright Broad Institute, 2018
## 
## This WDL converts paired FASTQ to uBAM and adds read group information 
##
## Requirements/expectations :
## - Pair-end sequencing data in FASTQ format (one file per orientation)
## - The following metada descriptors per sample:
## ```readgroup   fastq_pair1_file_path   fastq_pair2_file_path   sample_name   library_name   platform_unit   run_date   platform_name   sequecing_center``` 
##
## Outputs :
## - Set of unmapped BAMs, one per read group
## - File of a list of the generated unmapped BAMs
##
## Cromwell version support 
## - Successfully tested on v32
## - Does not work on versions < v23 due to output syntax
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

# WORKFLOW DEFINITION
workflow ConvertPairedFastQsToUnmappedBamWf {

  Array[String] sample_name 
  Array[String] fastq_1 
  Array[String] fastq_2 
  Array[String] readgroup_name 
  Array[String] library_name 
  Array[String] platform_unit 
  Array[String] run_date 
  Array[String] platform_name 
  Array[String] sequencing_center 

  String ubam_list_name

  String? gatk_docker_override
  String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:latest"])
  String? gatk_path_override
  String gatk_path = select_first([gatk_path_override, "gatk"])
  Int? preemptible_attempts

  # Convert multiple pairs of input fastqs in parallel
  scatter (i in range(length(readgroup_name))) {

    # Convert pair of FASTQs to uBAM
    call PairedFastQsToUnmappedBAM {
      input:
        sample_name = sample_name[i],
        fastq_1 = fastq_1[i],
        fastq_2 = fastq_2[i],
        readgroup_name = readgroup_name[i],
        library_name = library_name[i],
        platform_unit = platform_unit[i],
        run_date = run_date[i],
        platform_name = platform_name[i],
        sequencing_center = sequencing_center[i],
        gatk_path = gatk_path,
        docker = gatk_docker,
        preemptible_attempts = preemptible_attempts
    }
  }

  #Create a file with a list of the generated ubams
  call CreateFoFN {
    input:
      array_of_files = PairedFastQsToUnmappedBAM.output_bam,
      fofn_name = ubam_list_name,
      docker = gatk_docker
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] output_bams = PairedFastQsToUnmappedBAM.output_bam
    File unmapped_bam_list = CreateFoFN.fofn_list
  }
}

# TASK DEFINITIONS

# Convert a pair of FASTQs to uBAM
task PairedFastQsToUnmappedBAM {
  # Command parameters
  String sample_name
  File fastq_1
  File fastq_2
  String readgroup_name
  String library_name
  String platform_unit
  String run_date
  String platform_name
  String sequencing_center

  # Runtime parameters
  Int? disk_space_gb
  Int? machine_mem_gb
  Int? preemptible_attempts
  String docker
  String gatk_path

  command {
    ${gatk_path} --java-options "-Xmx3000m" \
    FastqToSam \
    --FASTQ ${fastq_1} \
    --FASTQ2 ${fastq_2} \
    --OUTPUT ${readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME ${readgroup_name} \
    --SAMPLE_NAME ${sample_name} \
    --LIBRARY_NAME ${library_name} \
    --PLATFORM_UNIT ${platform_unit} \
    --RUN_DATE ${run_date} \
    --PLATFORM ${platform_name} \
    --SEQUENCING_CENTER ${sequencing_center} 
  }
  runtime {
    #docker: docker
    memory: select_first([machine_mem_gb, 10]) + " GB"
    cpu: "1"
    disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_bam = "${readgroup_name}.unmapped.bam"
  }
}

task CreateFoFN {
  # Command parameters
  Array[String] array_of_files
  String fofn_name
  
  # Runtime parameters
  String docker
  
  command {
    mv ${write_lines(array_of_files)}  ${fofn_name}.list
  }
  output {
    File fofn_list = "${fofn_name}.list"
  }
  runtime {
    #docker: docker
    preemptible: 3
  }
}

