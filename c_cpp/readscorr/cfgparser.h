#ifndef _G_CFGPARSER_H
#define _G_CFGPARSER_H
#include <stdio.h>
#include <stdint.h>

void write_example_cfg(const char * const filename);

#define SDLConfig_Count 5u
typedef struct __SDLConfig {
    double GenomeSize,HetRatio,SeqErrRate,SeqDepth,MutEffect;
    char * seqfilename;
    size_t seqfilename_length;
    FILE *fp;
} SDLConfig;

SDLConfig *read_SDL_cfg(const char * const filename);
ssize_t get_next_seqfile(SDLConfig * const psdlcfg);
void destory_seqfile(SDLConfig * psdlcfg);

#endif /* cfgparser.h */
