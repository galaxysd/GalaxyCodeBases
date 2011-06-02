#include <stdio.h>
#include <stdlib.h> //EXIT_FAILURE
#include <err.h>
#include <string.h>
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

SDLConfig *read_SDL_cfg(const char * const filename){
    SDLConfig *psdlcfg = malloc(sizeof(SDLConfig));
    psdlcfg->fp = fopen(filename, "r");
    psdlcfg->seqfilename=NULL;
    psdlcfg->seqfilename_length=0u;
    if (psdlcfg->fp == NULL) err(EXIT_FAILURE, "Cannot write to [%s]", filename);
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    while ((read = getline(&line, &len, psdlcfg->fp)) != -1) {
        if (line[0]=='[' && line[read-2]==']') {
            *(line+(--read))='\0';  //chop
            break;
        }
    	if (*line=='\n' || line[0]=='#') continue;
        if (*(line+read-1)=='\n') *(line+(--read))='\0';    // line[read-1] = '\0';
    }
    if (! strcmp(line,"[Sequence Files]") ) {
        err(EXIT_FAILURE, "The first section is \"%s\", not \"[Sequence Files]\". See example input_config.", line);
    }
    return psdlcfg;
}

int_fast8_t get_next_seqfile(SDLConfig * const psdlcfg){

}

void destory_seqfile(SDLConfig * psdlcfg){
    fclose(psdlcfg->fp);
    free(psdlcfg->seqfilename);
    free(psdlcfg);
}
