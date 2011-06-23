#define _GNU_SOURCE
#include <stdlib.h> //EXIT_FAILURE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
//#include <errno.h>
#include <err.h>
#include <argp.h>
#include <math.h>
#include "getch.h"
#include "gFileIO.h"
#include "chrseq.h"
#include "timer.h"
#include "gtools.h"

#define MAXREADLEN (8ul*1024*1024)
uint64_t ReadsLenArr[MAXREADLEN];

const char *argp_program_version =
    "fcounter 0.1 @"__TIME__ "," __DATE__;
const char *argp_program_bug_address =
    "<huxuesong@genomics.org.cn>";

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
static char args_doc[] = "input_list (of FASTA or FASTQ files)";

/* The options we understand. */
static struct argp_option options[] = {
    {"interactive", 'i', 0,   0,  "Prompt arguments before procedure"},
    {"outfile",     'o', "./out.stat",0,  "Output file" },
//    {"countbase",   'c', 0,   0,  "Also count bases" },
    { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments {
    char *args[1];                /* arg1 */
    uint_fast8_t countbase, interactive;   //_Bool is the same under amd64, as typedef unsigned char uint_fast8_t;
    char *outfile;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state) {
/* Get the input argument from argp_parse, which we
  know is a pointer to our arguments structure. */
    struct arguments *arguments = state->input;
    
    switch (key) {
        case 'i':
            arguments->interactive = 1;
            break;
        case 'c':
            arguments->countbase = 1;
            break;
        case 'o':
            arguments->outfile = arg;
            break;
        
        case ARGP_KEY_ARG:
            if (state->arg_num >= 1)
             /* Too many arguments. */
             argp_usage (state);
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
    
    // Default values.
    arguments.countbase = 0;
    arguments.interactive = 0;
    arguments.outfile = "./out.stat";
    
    // Parse our arguments; every option seen by parse_opt will be reflected in arguments.
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
    
    if (arguments.interactive) {
      printf ("input_list = %s\nOutput = %s\nCountBases = %s\n",
           arguments.args[0],
           arguments.outfile,
           arguments.countbase ? "yes" : "no"
           );
      pressAnyKey();
    }
//printf("Out[%s][%s][%s]\n",outStat,outDat,outLog);

    uint64_t allbases=0;
    uint64_t allreads=0;
    uint64_t maxReadLen=0;
    float128 SS=0.0;
    double SStd;

    FILE *fp;
    ssize_t read;
    char *line = NULL;
    size_t len = 0;
    fp = fopen(arguments.args[0], "r");
    if (fp == NULL) err(EXIT_FAILURE, "Cannot open input_list [%s]", arguments.args[0]); //exit(EXIT_FAILURE);

    G_TIMER_START;

    fputs("\nParsing Sequence Files:\n", stderr);
     while ((read = getline(&line, &len, fp)) != -1) {
        ssize_t readlength;
        if (*line=='\n') continue;      // skip empty lines, which is definitely impossible
        if (*(line+read-1)=='\n') *(line+(--read))='\0';    // line[read-1] = '\0';
        fprintf(stderr, " <%s> ...", line);
        //sleep(1);    // the Call ...
        SeqFileObj *seqobj = inSeqFinit(line,GFIOCHRBASE);
        if (seqobj) {
        	while ( (readlength = (*seqobj->getNextSeq)(seqobj)) >= 0 ) {
        	#ifdef DEBUG
        		puts(line);
	        	printf("-ID:[%s,%s] %zu %zu %zd\nSeq:[%s]\nQ:[%s] *%zx,%u\n",
        			seqobj->name,seqobj->comment,seqobj->readlength,seqobj->binMallocedQQWord,readlength,
				    seqobj->seq,seqobj->qual,(size_t)seqobj->seq,seqobj->type);
				/*char *tmpseq;
			    //NormalizeChrSeq((char*)seqobj->seq);
			    printf("rc[%s]\n",tmpseq=ChrSeqRevComp(seqobj->seq,seqobj->readlength));
			    free(tmpseq);*/
			#endif
                allbases += readlength;
                SS += readlength*readlength;
                ++allreads;
                if (maxReadLen < readlength) maxReadLen = readlength;
                if (readlength < MAXREADLEN) {
                    ++ReadsLenArr[readlength];
                } else {
                    ++ReadsLenArr[0];
                }
        	}
        } else continue;
        fputs("\b\b\b\b, done !\n", stderr);
        inSeqFdestroy(seqobj);
    }

    fputs("\nCount Done!\n", stderr);
    free(line);
    fclose(fp);

    SStd = ( SS - (float128)allbases*(float128)allbases/(float128)allreads )
            / (double)(allreads -1);

    fp = fopen(arguments.outfile, "w");
    fprintf(fp,"#Total_Bases: %lu\n#Total_Reads: %lu\n#Avg_Read_Len: %.1f\tStd: %.3f\n"
        "#Max_Read_Len: %lu\n"
        "#Overflow: %lu\n"
        "\n#Read_Len\tCount\tRatio\n",
        allbases,allreads,
        0.05+((double)allbases/(double)allreads), SStd,
        maxReadLen,ReadsLenArr[0]);
    for (uint64_t i=1; i<=maxReadLen; i++) {
        fprintf(fp,"%lu\t%lu\t%g\n",i,ReadsLenArr[i],(double)ReadsLenArr[i]/(double)allreads);
        //printf("%lu\t%lu\n",i,ReadsLenArr[i]);
    }
    if (ReadsLenArr[0]) {
        fprintf(fp,"#>=%lu\t%lu\t%g\n",MAXREADLEN,ReadsLenArr[0],(double)ReadsLenArr[0]/(double)allreads);
    }
    fclose(fp);

    G_TIMER_END;
    G_TIMER_PRINT;

    exit(EXIT_SUCCESS);
}

