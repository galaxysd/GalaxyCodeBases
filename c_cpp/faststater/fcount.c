#define _GNU_SOURCE
#include <stdlib.h> //EXIT_FAILURE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <err.h>
#include <argp.h>
#include <math.h>
#include "gFileIO.h"
#include "chrseq.h"
#include "timer.h"

#define MAXREADLEN (64ul*1024*1024 - 1)
uint64_t ReadsLenArr[MAXREADLEN + 1];

const char *argp_program_version =
    "fcounter 0.2 @"__TIME__ "," __DATE__;
const char *argp_program_bug_address =
    "<galaxy001@gmail.com>";

/* Program documentation. */
static char doc[] =
    "Fasta/q Counter"
#ifdef DEBUG
    " (Debug Version)"
#endif
#ifdef TEST
    " (Test Version)"
#endif
;

/* A description of the arguments we accept. */
static char args_doc[] = "input_files (FASTA or FASTQ)";

/* The options we understand. */
static struct argp_option options[] = {
    {"outfile",     'o', "./out.stat",0,  "Output file" },
    { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments {
    char **args;
    char *outfile;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state) {
/* Get the input argument from argp_parse, which we
  know is a pointer to our arguments structure. */
    struct arguments *arguments = state->input;
    
    switch (key) {
        case 'o':
            arguments->outfile = arg;
            break;
        
        case ARGP_KEY_ARG:
            arguments->args[state->arg_num] = arg;
            break;
        
        case ARGP_KEY_END:
            if (state->arg_num < 1)   /* Not enough arguments. */
               //argp_usage (state);
               argp_state_help (state,stderr,ARGP_HELP_STD_HELP);
            break;
        
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

int main (int argc, char **argv) {
    struct arguments arguments;
    arguments.args=calloc(sizeof(size_t),argc);
    
    // Default values.
    arguments.outfile = "./out.stat";
    
    // Parse our arguments; every option seen by parse_opt will be reflected in arguments.
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
    
    uint64_t allbases=0;
    uint64_t allreads=0;
    uint64_t maxReadLen=0;
    uint64_t minReadLen=MAXREADLEN+1;
    float128 SS=0.0;
    double SStd;

    char *const*line = arguments.args-1;    // so that we can use *(++line)

    G_TIMER_START;

    fputs("\nParsing Sequence Files:\n", stderr);
    while ( *(++line) ) {
        ssize_t readlength;
        fprintf(stderr, " <%s> ...", *line);
        SeqFileObj *seqobj = inSeqFinit(*line,GFIOCHRBASE);
        if (seqobj) {
        	while ( (readlength = (*seqobj->getNextSeq)(seqobj)) >= 0 ) {
        	#ifdef DEBUG
	        	printf("-ID:[%s,%s] %zu %zu %zd\nSeq:[%s]\nQ:[%s] *%zx,%u\n",
        			seqobj->name,seqobj->comment,seqobj->readlength,seqobj->binMallocedQQWord,readlength,
				    seqobj->seq,seqobj->qual,(size_t)seqobj->seq,seqobj->type);
			#endif
                allbases += readlength;
                SS += readlength*readlength;
                ++allreads;
                if (maxReadLen < readlength) maxReadLen = readlength;
                if (minReadLen > readlength) minReadLen = readlength;
                if (readlength < MAXREADLEN) {
                    ++ReadsLenArr[readlength];
                } else {
                    ++ReadsLenArr[MAXREADLEN];
                }
        	}
        } else continue;
        fputs("\b\b\b\b, done !\n", stderr);
        inSeqFdestroy(seqobj);
    }
    free(arguments.args);
    fputs("\nCount Done!\n", stderr);

    SStd = sqrtl(( SS - (float128)allbases*(float128)allbases/(float128)allreads ) / (long double)(allreads -1));

    FILE *fp = fopen(arguments.outfile, "w");
    fprintf(fp,"#Total_Bases: %llu\n#Total_Reads: %llu\n#Avg_Read_Len: %.1f\tStd: %.3f\n"
        "#Read_Len_Range: [%llu,%llu]\n"
        "#Overflow: %llu\n"
        "\n#Read_Len\tCount\tRatio\n",
        allbases,allreads,
        0.05+((double)allbases/(double)allreads), SStd,
        minReadLen,maxReadLen,ReadsLenArr[0]);
    for (uint64_t i=minReadLen; i<=maxReadLen; i++) {
        if (ReadsLenArr[i])
            fprintf(fp,"%llu\t%llu\t%g\n",i,ReadsLenArr[i],(double)ReadsLenArr[i]/(double)allreads);
    }
    if (ReadsLenArr[MAXREADLEN]) {
        fprintf(fp,"#>=%lu\t%llu\t%g\n",MAXREADLEN,ReadsLenArr[MAXREADLEN],(double)ReadsLenArr[MAXREADLEN]/(double)allreads);
    }
    fclose(fp);

    G_TIMER_END;
    G_TIMER_PRINT;

    exit(EXIT_SUCCESS);
}

