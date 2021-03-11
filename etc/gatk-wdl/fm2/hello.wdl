task hello {
  String pattern
  File in

  command {
    egrep '${pattern}' '${in}'
  }

  runtime {
    #docker: "broadinstitute/my_image"
  }

  output {
    Array[String] matches = read_lines(stdout())
  }
}

workflow wf {
  call hello
}

