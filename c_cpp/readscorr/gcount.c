#define _GNU_SOURCE
#include <stdlib.h> //EXIT_FAILURE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
//#include <errno.h>
#include <err.h>
#include <argp.h>
#include "getch.h"
#include "gFileIO.h"
#include "chrseq.h"
#include "timer.h"

#define MAXHOMOPOLYLEN (2l*1024*1024)
#define BaseA 1
#define BaseT 3
#define BaseC 2
#define BaseG 4
#define BaseN 0
uint64_t HomoPoly[5][MAXHOMOPOLYLEN];
uint64_t maxHomoPoly[5];
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
      printf ("input_list = %s\nOutput = %s\n",
           arguments.args[0],
           arguments.outfile
           );
      pressAnyKey();
    }
//printf("Out[%s][%s][%s]\n",outStat,outDat,outLog);

    char lastbase;
    ssize_t lastpos,polymerLen;
    uint64_t allbases=0;

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
                            if (maxHomoPoly[theBase]<polymerLen) maxHomoPoly[theBase]=polymerLen;
                        } else {
                            ++HomoPoly[theBase][0];
                        }
                        lastbase=seqobj->seq[i];
                        lastpos=i;
                    }
                }
        	}
        } else continue;
        fputs("\b\b\b\b, done !\n", stderr);
        inSeqFdestroy(seqobj);
    }

    fputs("\nStat Done!\n", stderr);
    free(line);
    fclose(fp);

    fp = fopen(arguments.outfile, "w");
    fprintf(fp,"#Total_Bases: %lu\n"
        "\n#Base\tRepeat_Count\tTimes\n",allbases);
    for (unsigned char base=0;base<5;base++) {
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
        default:
            theBase='N';
            break;
        }
        for (size_t polymerLen=1;polymerLen<=maxHomoPoly[base];polymerLen++) {
            fprintf(fp,"%c\t%lu\t%lu\n",theBase,polymerLen,HomoPoly[base][polymerLen]);
        }
        if (HomoPoly[base][0]) {
            fprintf(fp,"#%c\t>=%lu\t%lu\n",theBase,MAXHOMOPOLYLEN,HomoPoly[base][0]);
        }
        fprintf(fp,"#----------------------------------------------------------------\n");
    }
    fclose(fp);

    G_TIMER_END;
    G_TIMER_PRINT;

    exit(EXIT_SUCCESS);
}

