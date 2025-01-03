version 1.0

# Copyright (c) 2018 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#import "BamMetrics/bammetrics.wdl" as bammetrics
#import "gatk-preprocess/gatk-preprocess.wdl" as preprocess
import "structs.wdl" as structs
import "bwa.wdl" as bwa
import "fgi.wdl" as fgifm2
import "sambamba.wdl" as sambamba
#import "tasks/picard.wdl" as picard
import "QC.wdl" as qc
#import "tasks/umi-tools.wdl" as umiTools

workflow SampleWorkflow {
	input {
		Sample sample
		String sampleDir
		String platform = "illumina"
		Boolean useBwaKit = false

		BwaIndex bwaIndex
		GatkIndex GatkIndex
		String? adapterForward
		String? adapterReverse
		Int? minBWAmatchLen = 200

		Int bwaThreads = 4
		Map[String, String] dockerImages
	}

	scatter (readgroup in sample.readgroups) {
		String readgroupDir = sampleDir + "/lib_" + readgroup.lib_id + "--rg_" + readgroup.id
		call qc.QC as qualityControl {
			input:
				outputDir = readgroupDir,
				read1 = readgroup.R1,
				read2 = readgroup.R2,
				adapterForward = adapterForward,
				adapterReverse = adapterReverse,
		}
		call bwa.Mem as bwaMem {
			input:
				read1 = qualityControl.qcRead1,
				read2 = qualityControl.qcRead2,
				outputPrefix = readgroupDir + "/" + sample.id + "-" + readgroup.lib_id + "-" + readgroup.id,
				readgroup = "@RG\\tID:~{sample.id}-~{readgroup.lib_id}-~{readgroup.id}\\tLB:~{readgroup.lib_id}\\tSM:~{sample.id}\\tPL:~{platform}",
				bwaIndex = select_first([bwaIndex]),
				threads = bwaThreads,
				usePostalt = useBwaKit,
				dockerImage = dockerImages["bwakit+samtools"]
		}
		Boolean paired = defined(readgroup.R2)

	}
	call sambamba.Markdup as markdup {
		input:
			inputBams = select_all(bwaMem.outputBam),
			outputPath = sampleDir + "/" + sample.id + ".markdup.bam",
			dockerImage = dockerImages["sambamba"]
	}

	call fgifm2.FilterSamReads as FilterSam {
		input:
			inputBam = markdup.outputBam,
			inputBamIndex = markdup.outputBamIndex,
			outputBamPath = sampleDir + "/" + sample.id + ".fM.bam",
			minMatchLen = minBWAmatchLen
	}
	call fgifm2.callSNP as fm2callSNP {
		input:
			inputBam = markdup.outputBam,
			inputBamIndex = markdup.outputBamIndex,
			GatkIndex = GatkIndex,
			outputPath = sampleDir + "/SNP/",
			helperPl = "bin/fsnp.pl"
	}
	call fgifm2.callSTR as fm2callSTR {
		input:
			inputBam = FilterSam.outputBam,
			inputBamIndex = FilterSam.outputBamIndex,
			outputPath = sampleDir + "/STR/",
			helperBED = "bin/LN-mid.bed",
			helperPl = "bin/fstr.pl"
	}

	output {
		File markdupBam = markdup.outputBam
		File markdupBamIndex = markdup.outputBamIndex
		File filteredBam = FilterSam.outputBam
		File filteredBamIndex = FilterSam.outputBamIndex
		File outSNP0txt = fm2callSNP.outSNP0txt
		File outSTR0txt = fm2callSTR.outSTR0txt
		File outSNP1txt = fm2callSNP.outSNP1txt
		File outSNPtxt = fm2callSNP.outSNPtxt
		File outSTRtxt = fm2callSTR.outSTRtxt
		Array[File] reports = flatten(qualityControl.reports)
		File outbcfFile = fm2callSNP.outbcfFile
		File outsnpFile = fm2callSNP.outsnpFile
		File outsnpIndexFile = fm2callSNP.outsnpIndexFile
		File outsnp0File = fm2callSNP.outsnp0File
		File outsnp0IndexFile = fm2callSNP.outsnp0IndexFile
	}

	parameter_meta {
		# inputs
		sample: {description: "The sample information: sample id, readgroups, etc.", category: "required"}
		sampleDir: {description: "The directory the output should be written to.", category: "required"}
		platform: {description: "The platform used for sequencing.", category: "advanced"}
		useBwaKit: {description: "Whether or not BWA kit should be used. If false BWA mem will be used.", category: "advanced"}
		bwaIndex: {description: "The BWA index files. These or the bwaMem2Index should be provided.", category: "common"}
		adapterForward: {description: "The adapter to be removed from the reads first or single end reads.", category: "common"}
		adapterReverse: {description: "The adapter to be removed from the reads second end reads.", category: "common"}
		bwaThreads: {description: "The amount of threads for the alignment process.", category: "advanced"}
		dockerImages: {description: "The docker images used.", category: "required"}

		# outputs
		markdupBam: {description: ""}
		markdupBamIndex: {description: ""}
		reports: {description: ""}
	}
}

