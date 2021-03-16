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
#import "tasks/bwa.wdl" as bwa
#import "tasks/bwa-mem2.wdl" as bwamem2
#import "tasks/sambamba.wdl" as sambamba
#import "tasks/picard.wdl" as picard
import "QC.wdl" as qc
#import "tasks/umi-tools.wdl" as umiTools

workflow SampleWorkflow {
	input {
		Sample sample
		String sampleDir
		String platform = "illumina"

		String? adapterForward
		String? adapterReverse

		Int bwaThreads = 4
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
	}
	output {
		Array[File] reports = flatten([flatten(qualityControl.reports),
		])
	}
}

