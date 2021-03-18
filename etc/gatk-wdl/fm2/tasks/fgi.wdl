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