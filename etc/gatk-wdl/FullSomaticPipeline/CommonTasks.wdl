# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
task SamToFastqAndBwaMemAndMba {
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
  Float disk_size
  Int compression_level
  Int preemptible_tries

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    # if ref_alt has data in it,
    if [ -s ${ref_alt} ]; then
      java -Xms5000m -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT=${input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      /usr/gitc/${bwa_commandline} /dev/stdin - 2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) | \
      java -Dsamjdk.compression_level=${compression_level} -Xms3000m -jar /usr/gitc/picard.jar \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=${input_bam} \
        OUTPUT=${output_bam_basename}.bam \
        REFERENCE_SEQUENCE=${ref_fasta} \
        PAIRED_RUN=true \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="bwamem" \
        PROGRAM_GROUP_VERSION="${bwa_version}" \
        PROGRAM_GROUP_COMMAND_LINE="${bwa_commandline}" \
        PROGRAM_GROUP_NAME="bwamem" \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false

      grep -m1 "read .* ALT contigs" ${output_bam_basename}.bwa.stderr.log | \
      grep -v "read 0 ALT contigs"

    # else ref_alt is empty or could not be found
    else
      exit 1;
    fi
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "14 GB"
    cpu: "16"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
  }
}

# Calculates sum of a list of floats
task SumFloats {
  Array[Float] sizes
  Int preemptible_tries

  command <<<
  python -c "print ${sep="+" sizes}"
  >>>
  output {
    Float total_size = read_float(stdout())
  }
  runtime {
    docker: "python:2.7"
    preemptible: preemptible_tries
  }
}

