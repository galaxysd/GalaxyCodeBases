#include <stdio.h>
#include <stdlib.h> //EXIT_FAILURE
#include <err.h>
#include "cfgparser.h"

void write_example_cfg(const char * const filename) {
    fputs(".\b", stderr);
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) err(EXIT_FAILURE, "Cannot write to [%s]", filename);
    fputs("\
# Estimated Genome Size in M bp\n\
Genome_Size=3000\n\
\n\
# Heterozygosity rate\n\
Het=0.01\n\
\n\
# DNA sequence error rate\n\
Seq_Err=0.02\n\
\n\
# Sequence Depth\n\
Seq_Depth=20\n\
\n\
[Sequence Files]\n\
./fullpath/to/fasta.fa\n\
./fullpath/to/fastaq.fq.any.filename.is.OK\n\
./fullpath/to/fasta.or.fastaq.gz\n\
", fp);
    if ( fclose(fp) )
        err(EXIT_FAILURE, "Error writing to [%s]", filename);
    fputs("\n", stderr);
}
