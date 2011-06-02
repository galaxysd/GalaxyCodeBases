#ifndef _G_CFGPARSER_H
#define _G_CFGPARSER_H

#define ERR_KMER_CORRECT_RATIO 0.95

void write_example_cfg(const char * const filename);

typedef struct __sCFG {
    double GenomeSize,HetRatio,SeqErrRate,SeqDepth;
} sConfig;



#endif /* cfgparser.h */
