#define _GNU_SOURCE
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h> //EXIT_FAILURE
#include "gkseq.h"

#include <stdint.h>
#include <string.h>
#include <err.h>
#include <argp.h>
#include <math.h>
#include "timer.h"

const char *argp_program_version =
    "readsdedump simple 0.1 @"__TIME__ "," __DATE__;
const char *argp_program_bug_address =
    "<huxuesong@genomics.org.cn>";

/* Program documentation. */
static char doc[] =
    "Simple readsdedump"
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

    double SS=0.0;
    double SStd;

    char *const*line = arguments.args-1;    // so that we can use *(++line)

    fputs("\nParsing Sequence Files:\n", stderr);
    G_TIMER_START;

	gzFile fpi;
	kseq_t *kseq;
	int_fast8_t read_type;

    while ( *(++line) ) {
        ssize_t readlength;
        fprintf(stderr, " <%s> ...", *line);
        fpi = gzopen(*line, "r");
        kseq = kseq_init(fpi);
        if (kseq) {
        	while ( (read_type = kseq_read(kseq)) >= 0 ) {
        	readlength = kseq->seq.l;
        	#ifdef DEBUG
	        	printf("-ID:[%s,%s] %zu %zu %zd\nSeq:[%s]\nQ:[%s] *%zx,%d\n",
        			kseq->name.s,kseq->comment.s,kseq->seq.l,kseq->qual.l,readlength,
				    kseq->seq.s,kseq->qual.s,(size_t)&(kseq->seq.s),read_type);
			#endif
                allbases += readlength;
                SS += readlength*readlength;
                ++allreads;

        	}
        } else continue;
        fputs("\b\b\b\b, done !\n", stderr);
	    kseq_destroy(kseq);
    	gzclose(fpi);
    }
    free(arguments.args);
    fputs("\nCount Done!\n", stderr);

    FILE *fp = fopen(arguments.outfile, "w");
    fprintf(fp,"#Total_Bases: %lu\n#Total_Reads: %lu\n#Avg_Read_Len: %.1f\tStd: %.3f\n"
        "#Read_Len_Range: [%lu,%lu]\n"
        "#Overflow: %lu\n"
        "\n#Read_Len\tCount\tRatio\n",
        allbases,allreads,
        0.05+((double)allbases/(double)allreads), SStd,
        0,maxReadLen,0);

    fclose(fp);

    G_TIMER_END;
    G_TIMER_PRINT;

    exit(EXIT_SUCCESS);
}

