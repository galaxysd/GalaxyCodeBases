version 1.0

task markdown {
	input {
		String ReportID
		Int libNum
		File makerPL
		File renderPY
		File template
		File addStyle
		File snptxt
		File bamStats
		String outputPath
		Array[File] quals
		Array[File] datas
	}

	String tmpFile = outputPath + "/" + ReportID + "_include.txt"
	String mdFile = outputPath + "/" + ReportID + ".md"
	String htmlFile = outputPath + "/" + ReportID + ".html"
	String pdfFile = outputPath + "/" + ReportID + ".pdf"
	command {
		set -e
		mkdir -p "~{outputPath}"

		echo "sample id	~{ReportID}" > ~{tmpFile}
		echo "library number	~{libNum}" >> ~{tmpFile}
		echo "snp table	~{snptxt}" >> ~{tmpFile}
		echo "align stats	~{bamStats}" >> ~{tmpFile}
		echo "quality control	~{sep=',' datas}" >> ~{tmpFile}
		echo "qual plot	~{sep=',' quals}" >> ~{tmpFile}
		perl ~{makerPL} ~{template} ~{tmpFile} > ~{mdFile}
		python3 ~{renderPY} -i ~{mdFile} -o ~{htmlFile} -s ~{addStyle}
		weasyprint ~{htmlFile} ~{pdfFile}
	}

	output {
		File mdReport = mdFile
		File htmlReport = htmlFile
		File pdfReport = pdfFile
	}
}
