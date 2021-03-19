version 1.0

import "tasks/sample.wdl" as sampleWf
import "tasks/common.wdl" as common
import "tasks/structs.wdl" as structs

workflow fm2 {
	input {
		File sampleConfigFile
		String outputDir = "."
		String platform = "illumina"
		Boolean useBwaKit = false
		Int scatterSizeMillions = 1000
		Int? minBWAmatchLen = 200
		BwaIndex bwaIndex
		GatkIndex GatkIndex
		String? adapterForward = "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"  # Illumina universal adapter.
		String? adapterReverse = "AAGTCGGATCGTAGCCATGTCGTTC"  # Illumina universal adapter.
		Int? scatterSize
		Int bwaThreads = 4
		File dockerImagesFile
	}

	# Parse docker Tags configuration and sample sheet.
	call common.YamlToJson as convertDockerImagesFile {
		input:
			yaml = dockerImagesFile,
			outputJson = "dockerImages.json"
	}

	Map[String, String] dockerImages = read_json(convertDockerImagesFile.json)
	#Array[Array[File]] inputSamples = read_tsv(SampleFqTSV)
	SampleConfig sampleConfig = read_json(sampleConfigFile)

	scatter (sample in sampleConfig.samples) {
		String sampleIds = sample.id
		String sampleDir = outputDir + "/samples/" + sample.id + "/"
		call sampleWf.SampleWorkflow as sampleWorkflow {
			input:
				sampleDir = sampleDir,
				sample = sample,
				bwaThreads = bwaThreads,
				bwaIndex = bwaIndex,
				#bwaMem2Index = bwaMem2Index,
				GatkIndex = GatkIndex,
				adapterForward = adapterForward,
				adapterReverse = adapterReverse,
				useBwaKit = useBwaKit,
				dockerImages = dockerImages,
				#scatters = scatterList.scatters,
				bwaThreads = bwaThreads,
				platform = platform,
				minBWAmatchLen = minBWAmatchLen
				#sampleName = sample.id,
				#gender = select_first([sample.gender, "unknown"]),
		}
	}

	Array[File] allReports = flatten([
		flatten(sampleWorkflow.reports)
		#, flatten(singleSampleCalling.reports)
	])
	
	output {
		File dockerImagesList = convertDockerImagesFile.json
		#File multiqcReport = multiqcTask.multiqcReport
		Array[File] reports = allReports
		#Array[File] recalibratedBams = sampleWorkflow.recalibratedBam
		#Array[File] recalibratedBamIndexes = sampleWorkflow.recalibratedBamIndex
		Array[File] markdupBams = sampleWorkflow.markdupBam
		Array[File] markdupBamIndexes = sampleWorkflow.markdupBamIndex
		Array[File] filteredBam = sampleWorkflow.filteredBam
		Array[File] filteredBamIndex = sampleWorkflow.filteredBamIndex
		Array[File] outSNPtxts = sampleWorkflow.outSNPtxt
		Array[File] outSTRtxts = sampleWorkflow.outSTRtxt
		#Array[File?] mantaVCFs = svCalling.mantaVcf
		#Array[File?] dellyVCFs = svCalling.dellyVcf
		#Array[File?] survivorVCFs = svCalling.survivorVcf
		#Array[Array[File]?] modifiedVcfs = svCalling.modifiedVcfs
	}

	parameter_meta {
		# inputs
		sampleConfigFile: {description: "The samplesheet, including sample ids, library ids, readgroup ids and fastq file locations.", category: "required"}
		outputDir: {description: "The directory the output should be written to.", category: "common"}
		platform: {description: "The platform used for sequencing.", category: "advanced"}
		useBwaKit: {description: "Whether or not BWA kit should be used. If false BWA mem will be used.", category: "advanced"}
		scatterSizeMillions:{description: "Same as scatterSize, but is multiplied by 1000000 to get scatterSize. This allows for setting larger values more easily.", category: "advanced"}
		bwaIndex: {description: "The BWA index files. When these are provided BWA will be used.", category: "common"}
		adapterForward: {description: "The adapter to be removed from the reads first or single end reads.", category: "common"}
		adapterReverse: {description: "The adapter to be removed from the reads second end reads.", category: "common"}
		scatterSize: {description: "The size of the scattered regions in bases for the GATK subworkflows. Scattering is used to speed up certain processes. The genome will be seperated into multiple chunks (scatters) which will be processed in their own job, allowing for parallel processing. Higher values will result in a lower number of jobs. The optimal value here will depend on the available resources.", category: "advanced"}
		bwaThreads: {description: "The amount of threads for the alignment process.", category: "advanced"}
		dockerImagesFile: {description: "A YAML file describing the docker image used for the tasks. The dockerImages.yml provided with the pipeline is recommended.", category: "advanced"}

		# outputs
		dockerImagesList: {description: "Json file describing the docker images used by the pipeline."}
		reports: {description: ""}
		recalibratedBams: {description: ""}
		recalibratedBamIndexes: {description: ""}
		markdupBams: {description: ""}
		markdupBamIndexes: {description: ""}
	}
}
