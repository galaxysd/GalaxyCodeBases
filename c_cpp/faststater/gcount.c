#define _GNU_SOURCE
#include <stdlib.h> //EXIT_FAILURE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <err.h>
#include <argp.h>
#include "gFileIO.h"
#include "chrseq.h"
#include "timer.h"

#define MAXHOMOPOLYLEN (128ul*1024)
#define BaseA 1
#define BaseT 3
#define BaseC 2
#define BaseG 4
#define BaseATCG 5
#define BaseN 0
uint64_t HomoPoly[6][MAXHOMOPOLYLEN];
uint64_t maxHomoPoly[6],BaseCounter[5];
unsigned char theBase;

const char *argp_program_version =
    "HomoPolymer Counter 0.1 @"__TIME__ "," __DATE__;
const char *argp_program_bug_address =
    "<huxuesong@genomics.org.cn>";

/* Program documentation. */
static char doc[] =
    "HomoPolymer Counter"
#ifdef DEBUG
    " (Debug Version)"
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
    
    char lastbase;
    ssize_t lastpos,polymerLen;
    uint64_t allbases=0;

    char *const*line = arguments.args-1;

    G_TIMER_START;

    fputs("\nParsing Sequence Files:\n", stderr);
    while ( *(++line) ) {
        ssize_t readlength;
        fprintf(stderr, " <%s>:\n", *line);
        SeqFileObj *seqobj = inSeqFinit(*line,GFIOCHRBASE);
        if (seqobj) {
        	while ( (readlength = (*seqobj->getNextSeq)(seqobj)) >= 0 ) {
        	#ifdef DEBUG
	        	printf("-ID:[%s,%s] %zu %zu %zd\nSeq:[%s]\nQ:[%s] *%zx,%u\n",
        			seqobj->name,seqobj->comment,seqobj->readlength,seqobj->binMallocedQQWord,readlength,
				    seqobj->seq,seqobj->qual,(size_t)seqobj->seq,seqobj->type);
			#endif
			    fprintf(stderr, "  >%s %s,%zu ...",seqobj->name,seqobj->comment,seqobj->readlength);
                lastbase=seqobj->seq[0];
                lastpos=0;
                allbases += readlength;
                for (ssize_t i=1;i<=readlength;i++) {   // *(readlength+1)=='\0'
                    if (seqobj->seq[i] != lastbase) {
                        polymerLen = i-lastpos;
                        switch (lastbase) {
                        case 'A':
                            theBase=BaseA;
                            break;
                        case 'T':
                            theBase=BaseT;
                            break;
                        case 'C':
                            theBase=BaseC;
                            break;
                        case 'G':
                            theBase=BaseG;
                            break;
                        default:
                            theBase=BaseN;
                            break;
                        }
                        if (polymerLen<MAXHOMOPOLYLEN) {
                            ++HomoPoly[theBase][polymerLen];
                            if (lastbase != 'N') ++HomoPoly[BaseATCG][polymerLen];
                            if (maxHomoPoly[theBase]<polymerLen) {
                                maxHomoPoly[theBase]=polymerLen;
                                if (lastbase != 'N' && maxHomoPoly[BaseATCG]<polymerLen)
                                    maxHomoPoly[BaseATCG]=polymerLen;
                            }
                        } else {
                            ++HomoPoly[theBase][0];
                            if (lastbase != 'N') ++HomoPoly[BaseATCG][0];
                        }
                        BaseCounter[theBase] += polymerLen;
                        lastbase=seqobj->seq[i];
                        lastpos=i;
                    }
                }
                fputs("\b\b\b\b, done !\n", stderr);
        	}
        } else continue;
        inSeqFdestroy(seqobj);
    }
    free(arguments.args);
    fputs("\nStat Done!\n", stderr);

    FILE *fp = fopen(arguments.outfile, "w");
    fprintf(fp,"# Total_Bases: %llu\n"
        "# A: %llu\tT: %llu\tC: %llu\tG: %llu\tN: %llu\n"
        "# A%%: %.3f\tT%%: %.3f\tC%%: %.3f\tG%%: %.3f\tN%%: %.3f\n"
        "\n#Base\tRepeat_Count\tTimes\tBaseRatio\n",
        allbases,
        BaseCounter[BaseA],BaseCounter[BaseT],BaseCounter[BaseC],BaseCounter[BaseG],BaseCounter[BaseN],
        100*(double)BaseCounter[BaseA]/allbases,100*(double)BaseCounter[BaseT]/allbases,
        100*(double)BaseCounter[BaseC]/allbases,100*(double)BaseCounter[BaseG]/allbases,
        100*(double)BaseCounter[BaseN]/allbases);
    for (unsigned char base=0;base<6;base++) {
        switch (base) {
        case BaseA:
            theBase='A';
            break;
        case BaseT:
            theBase='T';
            break;
        case BaseC:
            theBase='C';
            break;
        case BaseG:
            theBase='G';
            break;
        case BaseATCG:
            theBase='H';
            break;
        default:
            theBase='N';
            break;
        }
        for (size_t polymerLen=1;polymerLen<=maxHomoPoly[base];polymerLen++) {
            if (HomoPoly[base][polymerLen])
                fprintf(fp,"%c\t%lu\t%llu\t%g\n",theBase,polymerLen,HomoPoly[base][polymerLen],
                (double)HomoPoly[base][polymerLen]*polymerLen/allbases);
        }
        if (HomoPoly[base][0]) {
            fprintf(fp,"#%c\t>=%lu\t%llu\n",theBase,MAXHOMOPOLYLEN,HomoPoly[base][0]);
        }
        fprintf(fp,"#----------------------------------------------------------------\n");
    }
    fclose(fp);

    G_TIMER_END;
    G_TIMER_PRINT;

    exit(EXIT_SUCCESS);
}

