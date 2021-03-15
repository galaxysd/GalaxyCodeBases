workflow fm2 {
	String SampleFqTSV
	Array[Array[File]] inputSamples = read_tsv(SampleFqTSV)
	
	scatter (sample in inputSamples) {
		call Cutadapter {
			input: sampleID=sample[0],
			sampleName=sample[1], 
			sampleFQ1=sample[2],
			sampleFQ2=sample[3]
		}
	}
	call mergeFQ {
		input: cFQs=Cutadapter.cfq
	}
}

task Cutadapter {
	String sampleID
	String sampleName
	String sampleFQ1
	String sampleFQ2

	command {
		mkdir -p fq/${sampleName}
		echo "${sampleID} ${sampleName} ${sampleFQ1} ${sampleFQ2}" > fq/${sampleName}/${sampleID}.fq
	}
	output {
		File cfq = "fq/${sampleName}/${sampleID}.fq"
	}
}

task mergeFQ {
	Array[File] cFQs

	command {
		mkdir -p fq
		echo "${sep=", " cFQs}" >fq/merge.txt
	}
	output {
		File mfq = "fq/merge.txt"
	}
}