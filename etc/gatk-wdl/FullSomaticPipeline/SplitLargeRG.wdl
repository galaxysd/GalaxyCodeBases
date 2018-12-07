import "CommonTasks.wdl" as tasks

workflow SplitLargeRG {
  File input_bam

  String bwa_commandline
  String bwa_version
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative".
  File ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  Int additional_disk
  Int compression_level
  Int preemptible_tries
  Int reads_per_file = 48000000

  Float bwa_ref_size
  Float disk_multiplier

  Float unmapped_bam_size

  call SamSplitter {
    input :
      input_bam = input_bam,
      n_reads = reads_per_file,
      # Since the output bams are less compressed than the input bam we need a disk multiplier
      # that's larger than 2.
      disk_size = ceil(disk_multiplier * unmapped_bam_size + additional_disk),
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }

  scatter(unmapped_bam in SamSplitter.split_bams) {
    Float current_unmapped_bam_size = size(unmapped_bam, "GB")
    String current_name = basename(unmapped_bam, ".bam")

    call tasks.SamToFastqAndBwaMemAndMba as Alignment {
      input:
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = current_name,
        ref_fasta = ref_fasta,
         ref_fasta_index = ref_fasta_index,
         ref_dict = ref_dict,
         ref_alt = ref_alt,
         ref_bwt = ref_bwt,
         ref_amb = ref_amb,
         ref_ann = ref_ann,
         ref_pac = ref_pac,
         ref_sa = ref_sa,
         bwa_version = bwa_version,
         # The merged bam can be bigger than only the aligned bam,
         # so account for the output size by multiplying the input size by 2.75.
         disk_size = current_unmapped_bam_size + bwa_ref_size + (disk_multiplier * current_unmapped_bam_size) + additional_disk,
         compression_level = compression_level,
         preemptible_tries = preemptible_tries
    }

    Float current_mapped_size = size(Alignment.output_bam, "GB")
  }

  call tasks.SumFloats as SumSplitAlignedSizes {
    input:
      sizes = current_mapped_size,
      preemptible_tries = preemptible_tries
  }

  call GatherBamFiles {
    input:
      input_bams = Alignment.output_bam,
      disk_size = ceil((2 * SumSplitAlignedSizes.total_size) + additional_disk),
      output_bam_basename = output_bam_basename,
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }

  output {
    File aligned_bam = GatherBamFiles.output_bam
  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  Array[File] input_bams
  String output_bam_basename
  Int disk_size
  Int compression_level
  Int preemptible_tries

  command {
    java -Dsamjdk.compression_level=${compression_level} -Xms3000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=${sep=' INPUT=' input_bams} \
      OUTPUT=${output_bam_basename}.bam 
  }
  runtime {
    preemptible: preemptible_tries
    memory: "3.75 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}

task SamSplitter {
  File input_bam
  Int n_reads
  Int disk_size
  Int preemptible_tries
  Int compression_level

  command {
    set -e
    mkdir output_dir

    total_reads=$(samtools view -c ${input_bam})

    java -Dsamjdk.compression_level=${compression_level} -Xms3000m -jar /usr/gitc/picard.jar SplitSamByNumberOfReads \
      INPUT=${input_bam} \
      OUTPUT=output_dir \
      SPLIT_TO_N_READS=${n_reads} \
      TOTAL_READS_IN_INPUT=$total_reads
  }
  output {
    Array[File] split_bams = glob("output_dir/*.bam")
  }
  runtime {
    preemptible: preemptible_tries
    memory: "3.75 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
}

