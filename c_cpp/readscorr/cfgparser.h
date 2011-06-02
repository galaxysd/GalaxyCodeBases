#ifndef _G_CFGPARSER_H
#define _G_CFGPARSER_H
#include <stdio.h>
#include <stdint.h>

#define ERR_KMER_CORRECT_RATIO 0.95

void write_example_cfg(const char * const filename);

typedef struct __SDLConfig {
    double GenomeSize,HetRatio,SeqErrRate,SeqDepth;
    char * seqfilename;
    size_t seqfilename_length;
    FILE *fp;
} SDLConfig;

SDLConfig *read_SDL_cfg(const char * const filename);
int_fast8_t get_next_seqfile(SDLConfig * const psdlcfg);
void destory_seqfile(SDLConfig * psdlcfg);

#endif /* cfgparser.h */
