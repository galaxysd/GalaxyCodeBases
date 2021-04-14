version 1.0

task FilterSamReads {
	input {
		File inputBam
		File inputBamIndex
		String outputBamPath
		Int minMatchLen = 200
		Boolean? createMd5File = false

		Int? compressionLevel

		Int? javaXmxMb = 1024
		Int? memoryMb = javaXmxMb + 512
		# One minute per input gigabyte.
		Int? timeMinutes = 1 + ceil(size(inputBam, "G") * 1)
		String? dockerImage = "quay.io/biocontainers/picard:2.23.8--0"
	}
	Array[String] FilterScriptContents = [
		'function accept(rec) {',
		'  if (rec.getReadUnmappedFlag()) return false;',
		'  var cigar = rec.getCigar();',
		'  if (cigar == null) return false;',
		'  var readMatch = 0;',
		'  for (var i=0;i < cigar.numCigarElements();++i) {',
		'    var ce = cigar.getCigarElement(i);',
		'    if (ce.getOperator().name() == "M") readMatch += ce.length;',
		'  }',
		'  if (readMatch > '+ minMatchLen +') return true;',
		'}',
		'accept(record);'
	]

	command {
		set -e
		mkdir -p "$(dirname ~{outputBamPath})"
		JAVA_OPTS="-Xmx~{javaXmxMb}M -XX:ParallelGCThreads=1" picard \
		FilterSamReads \
		INPUT=~{sep=' INPUT=' inputBam} \
		OUTPUT=~{outputBamPath} \
		~{"COMPRESSION_LEVEL=" + compressionLevel} \
		JAVASCRIPT_FILE=${write_lines(FilterScriptContents)} \
		FILTER=includeJavascript \
		CREATE_INDEX=true \
		CREATE_MD5_FILE=~{true="true" false="false" createMd5File}
	}

	output {
		File outputBam = outputBamPath
		File outputBamIndex = sub(outputBamPath, "\.bam$", ".bai")
		File? outputBamMd5 = outputBamPath + ".md5"
	}

	runtime {
		memory: "~{memoryMb}M"
		time_minutes: timeMinutes
		#docker: dockerImage
	}
}

task callSNP {
	input {
		File inputBam
		File inputBamIndex
		GatkIndex GatkIndex
		String outputPath
		File helperPl
	}
	File referenceFasta = GatkIndex.fastaFile
	File dbsnpVCF = GatkIndex.dbsnpVCF
	String bcfFile = outputPath + "/mpileup.bcf"
	String snp0File = outputPath + "/snp0.gz"
	String snpFile = outputPath + "/snp.gz"
	File SNPosFile = "bin/snpos.lst"
	command {
		set -e
		mkdir -p "~{outputPath}"
		bcftools mpileup --threads 6 ~{inputBam} -d 30000 -Q 30 -f ~{referenceFasta} -p -Ob \
		-a FORMAT/AD,FORMAT/SCR,FORMAT/ADF,FORMAT/ADR \
		--ff UNMAP,SECONDARY,QCFAIL -B \
		-o ~{bcfFile}
		bcftools call -Oz -A -c -p 0.9 ~{bcfFile} -o ~{snp0File}
		bcftools index ~{snp0File}
		bcftools annotate -a ~{dbsnpVCF} ~{snp0File} -c ID --collapse all -R ~{SNPosFile} -Oz -o ~{snpFile}
		bcftools index ~{snpFile}
		bcftools query -f'%CHROM\t[%DP\t%QUAL\t%TGT\n]' -i 'POS==501' ~{snpFile} > ~{outputPath + "/snp0.txt"}
		bcftools query -f'%ID\t[%DP\t%QUAL\t%TGT]\t%CHROM:%POS\n' ~{snpFile} -o ~{outputPath + "/../snpG.txt"}
		perl ~{helperPl} ~{outputPath + "/snp0.txt"} > ~{outputPath + "/../snp.txt"}
	}
	# See <https://github.com/samtools/bcftools/issues/658> for `-c -p 0.9`. This fix low recalculated BaseQ next to INDEL.

	output {
		File outSNP0txt = outputPath + "/snp0.txt"
		File outSNP1txt = outputPath + "/../snpG.txt"
		File outSNPtxt = outputPath + "/../snp.txt"
		File outbcfFile = bcfFile
		File outsnpFile = snpFile
		File outsnpIndexFile = snpFile + ".csi"
	}
}

task callSTR {
	input {
		File inputBam
		File inputBamIndex
		String outputPath
		File helperBED
		File helperPl
	}
	String awkStr = "'\{print $1\"\\t\"$4\}'"
	command {
		set -e
		mkdir -p "~{outputPath}"
		samtools mpileup -l ~{helperBED} ~{inputBam} |awk ~{awkStr} > ~{outputPath + "/str0.txt"}
		perl ~{helperPl} ~{outputPath + "/str0.txt"} > ~{outputPath + "/../str.txt"}
	}
	output {
		File outSTR0txt = outputPath + "/str0.txt"
		File outSTRtxt = outputPath + "/../str.txt"
	}
}

struct GatkIndex {
	File fastaFile
	File dbsnpVCF
	Array[File] indexFiles
}
