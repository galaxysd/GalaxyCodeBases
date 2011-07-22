#define _GNU_SOURCE
#include <stdlib.h> //EXIT_FAILURE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
//#include <errno.h>
#include <err.h>
#include <argp.h>
#include <math.h>
#include "MurmurHash3.h"
#include "getch.h"
#include "timer.h"

#define SUBARRAY_MAX 8192u

const char *argp_program_version =
    "depthstator 0.1 @"__TIME__ "," __DATE__;
const char *argp_program_bug_address =
    "<huxuesong@genomics.org.cn>";

/* Program documentation. */
static char doc[] =
    "a ???"
#ifdef DEBUG
    " (Debug Version)"
#endif
#ifdef TEST
    " (Test Version)"
#endif
;

/* A description of the arguments we accept. */
static char args_doc[] = "input_config (of nfo & FASTA or FASTQ files)";

/* The options we understand. */
static struct argp_option options[] = {
    {"interactive", 'i', 0,   0,  "Prompt arguments before procedure"},
    {"outprefix",'o', "./out",0,  "Output to [./out.{dat,stat,log}]" },
    {"ref", 'r', NULL,   0,  "Reference fasta file" },
    { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments {
    char **args;
    uint_fast8_t interactive;
    char *outprefix, *refile;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state) {
/* Get the input argument from argp_parse, which we
  know is a pointer to our arguments structure. */
    struct arguments *arguments = state->input;
    
    int tmpArgValue;
    switch (key) {
        case 'i':
            arguments->interactive = 1;
            break;
        case 'r':
            arguments->refile = arg;
            break;
        case 'o':
            arguments->outprefix = arg;
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

char *strlinker(char *main, char *suffix) {
    char *outstr = malloc(strlen(main)+strlen(suffix)+1);
    strcpy(outstr,main);
    strcat(outstr,suffix);
    return outstr;
}

int main (int argc, char **argv) {
    struct arguments arguments;
    arguments.args=calloc(sizeof(size_t),argc);

    // Default values.
    arguments.interactive = 0;
    arguments.refile  = NULL;
    arguments.outprefix = "./out";
    
    // Parse our arguments; every option seen by parse_opt will be reflected in arguments.
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
    
    if (arguments.interactive) {
      printf ("input_config = %s\nOutputPrefix = %s\n"
           "SDLA { CountBit:%u. rBit:%u(%u), ArraySize:%lu(k), SubItemCount:%u(%u) }\n"
           "kmersize = %i, SDLAsize = %.1f MiB (%s Bytes)\n"
           "HashBit:%.3g, ErrRate:%.12g\nMaxCapacity: %s\n",
           arguments.args[0],
           arguments.outprefix,
           arguments.CountBit,rBit,arguments.rBit,arguments.ArraySizeK,SubItemCount,arguments.SubItemCount,
           arguments.kmersize,(double)SDLsize/1048576.0,commaprint(SDLsize),
           HashBit,pow(0.5,rBit)/ArraySize,commaprint(Capacity)
           );
      pressAnyKey();
    }
    char *outStat = strlinker(arguments.outprefix, ".stat");
    char *outDat  = strlinker(arguments.outprefix, ".dat");
    char *outLog  = strlinker(arguments.outprefix, ".log");
//printf("Out[%s][%s][%s]\n",outStat,outDat,outLog);
    SDLConfig *psdlcfg=read_SDL_cfg((double)arguments.kmersize,arguments.args[0]);
    ssize_t read;
    char *line;
    uint64_t insertedCount=0;

    SDLeftArray_t *dleftp = dleft_arrayinit(arguments.CountBit,arguments.rBit,arguments.ArraySizeK*1024,arguments.SubItemCount); // (9,31,5000000,48);
    fputs("SDLA nfo: ", stderr);
    fprintSDLAnfo(stderr,dleftp);

    G_TIMER_START;

    fputs("\nParsing Sequence Files:\n", stderr);
    while ((read = get_next_seqfile(psdlcfg)) != -1) {
    	ssize_t readlength;
    	line=psdlcfg->seqfilename;
        fprintf(stderr, " <%s> ...", line);
        //sleep(1);    // the Call ...
        SeqFileObj *seqobj = inSeqFinit(line,GFIOCHRBASE);
        if (seqobj) {
        	while ( (readlength = (*seqobj->getNextSeq)(seqobj)) >= 0 ) {
        	#ifdef DEBUG
        		puts(line);
	        	printf("-ID:[%s,%s] %zu %zu\nSeq:[%s]\nQ:[%s] *%zx,%u\n",
        			seqobj->name,seqobj->comment,seqobj->readlength,seqobj->binMallocedQQWord,
				    seqobj->seq,seqobj->qual,(size_t)seqobj->seq,seqobj->type);
				/*char *tmpseq;
			    //NormalizeChrSeq((char*)seqobj->seq);
			    printf("rc[%s]\n",tmpseq=ChrSeqRevComp(seqobj->seq,seqobj->readlength));
			    free(tmpseq);*/
			#endif
                size_t insertedKmer = dleft_insert_read(arguments.kmersize,seqobj->seq,seqobj->readlength,dleftp);
                insertedCount += insertedKmer;
            #ifdef DEBUG
                printf("Inserted:[%zu of %lu]\n",insertedKmer,insertedCount);
            #endif
        	}
        } else continue;
        fputs("\b\b\b\b, done !\n", stderr);
        fprintf(stderr, "[!]Accumulated Inserted Kmer:[%lu] times\n", insertedCount);
        inSeqFdestroy(seqobj);
    }
/*
    rewind(fp);
    fputs("\nSecond pass:\n", stderr);
    while ((read = getline(&line, &len, fp)) != -1) {
        if (*(line+read-1)=='\n') *(line+(--read))='\0';
        printf("[%zu]+%zu <%s>\n", read, len, line);
    }
*/
    fputs("\nHashing Done!\n", stderr);
    //__clock_diff = clock() - __clock_start;
    //fprintf(stderr,"\nHashing Done with CPU time: %zd.%06zd s\n", __clock_diff/CLOCKS_PER_SEC, __clock_diff%CLOCKS_PER_SEC);

    destory_seqfile(psdlcfg);
    fputs("SDLA nfo: ", stderr);
    fprintSDLAnfo(stderr,dleftp);
    FILE *fpstat, *fpdat;
    fpstat = fopen(outStat, "w");
    fpdat = fopen(outDat, "w");
    fprintf(fpstat,"#KmerSize: %d\n#Kmer_theory_count: %.34Lg\n",arguments.kmersize,powl(4.0,(long double)arguments.kmersize));   // 34 from http://en.wikipedia.org/wiki/Quadruple_precision

#ifdef TEST    /* in TEST mode, "out.stat" changed to subArray filling count */
    SDLeftStat_t *SDLeftStat = dleft_stat(dleftp,fpstat,fpdat);
#else           /* Normal out.stat  */
    SDLeftStat_t *SDLeftStat = dleft_stat(dleftp,fpstat);

    SDLdumpHead *DatHead = malloc(sizeof(SDLdumpHead));
    DatHead->kmersize = arguments.kmersize;
    DatHead->HistMaxCntVal = SDLeftStat->HistMaxCntVal;
    DatHead->HistMaxHistVal = SDLeftStat->HistMaxHistVal;
    DatHead->HistMean = SDLeftStat->HistMean;
    DatHead->HistSStd = SDLeftStat->HistSStd;
    //free(SDLeftStat);
    dleft_dump(dleftp,DatHead,fpdat);
    free(DatHead);
#endif
    fclose(fpstat);
    fclose(fpdat);

    free(SDLeftStat);
    dleft_arraydestroy(dleftp);
    free(outStat); free(outDat); free(outLog);  // just filenames

    G_TIMER_END;
    G_TIMER_PRINT;

    exit(EXIT_SUCCESS);
}
