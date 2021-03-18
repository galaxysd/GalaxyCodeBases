version 1.0

task FilterSamReads {
	input {
		File inputBam
		File inputBamIndex
		String outputBamPath
		Boolean? createMd5File = false

		Int? compressionLevel

		Int? javaXmxMb = 1024
		Int? memoryMb = javaXmxMb + 512
		# One minute per input gigabyte.
		Int? timeMinutes = 1 + ceil(size(inputBam, "G") * 1)
		String? dockerImage = "quay.io/biocontainers/picard:2.23.8--0"
	}
	String FilterScriptContent = "function accept(e){if(e.getReadUnmappedFlag())return!1;var r=e.getCigar();if(null==r)return!1;for(var t=0,a=0;a<r.numCigarElements();++a){var n=r.getCigarElement(a);\\"M\\"==n.getOperator().name()&&(t+=n.length)}return 200<t||void 0}accept(record);"

	command {
		set -e
		mkdir -p "$(dirname ~{outputBamPath})"
		picard -Xmx~{javaXmxMb}M -XX:ParallelGCThreads=1 \
		FilterSamReads \
		INPUT=~{sep=' INPUT=' inputBam} \
		OUTPUT=~{outputBamPath} \
		~{"COMPRESSION_LEVEL=" + compressionLevel} \
		JAVASCRIPT_FILE=${write_lines([FilterScriptContent])} \
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