version 1.0

import "tasks/sample.wdl" as sampleWf
import "tasks/common.wdl" as common
import "tasks/structs.wdl" as structs

workflow fm2 {
	input {
		File sampleConfigFile
		String outputDir = "."
		Int? scatterSize
		Int bwaThreads = 4
	}

	#Array[Array[File]] inputSamples = read_tsv(SampleFqTSV)
	SampleConfig sampleConfig = read_json(sampleConfigFile)

	scatter (sample in sampleConfig.samples) {
		String sampleIds = sample.id
		String sampleDir = outputDir + "/samples/" + sample.id + "/"
		call sampleWf.SampleWorkflow as sampleWorkflow {
			input:
				sampleDir = sampleDir,
				sample = sample,
				bwaThreads = bwaThreads
				#sampleName = sample.id,
				#gender = select_first([sample.gender, "unknown"]),
		}
	}

	output {
		File multiqcReport = sampleConfigFile
	}
}













task Cutadapter {
	input {
		String sampleID
		String sampleName
		String sampleFQ1
		String sampleFQ2
	}

	command {
		mkdir -p fq/${sampleName}
		echo "${sampleID} ${sampleName} ${sampleFQ1} ${sampleFQ2}" > fq/${sampleName}/${sampleID}.fq
	}

	output {
		File cfq = "fq/${sampleName}/${sampleID}.fq"
	}
}

task mergeFQ {
	input{
		Array[File] cFQs
	}

	command {
		mkdir -p fq
		echo "${sep=", " cFQs}" >fq/merge.txt
	}
	output {
		File mfq = "fq/merge.txt"
	}
}