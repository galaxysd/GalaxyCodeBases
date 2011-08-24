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

#include "prb.h"
#include "prb-tree.h"

#define CAT0(x) #x
#define CAT(x) CAT0(x)
// 1048576 + 128k
//#define READSCOUNT_INIT (1048576)
//#define READSCOUNT_INC  (128*1024)
#define READSCOUNT_INIT (16)
#define READSCOUNT_INC  (4)

const char *argp_program_version =
    "readsdedump simple 0.1 @"__TIME__ "," __DATE__;
const char *argp_program_bug_address =
    "<huxuesong@genomics.org.cn>";

/* Program documentation. */
static char doc[] =
    "Seeding readsdedump"
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
    {"maxmismatch", 'm', "10"           ,0, "Max mismatch to be PCR duplicate"},
    {"seedlength",  's', "7"            ,0, "Seed length"},
    {"outfile",     'o', "./out.stat"   ,0, "Output file" },
    { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments {
    uint_fast16_t mismatch,seedLen;
    char **args;
    char *outfile;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state) {
/* Get the input argument from argp_parse, which we
  know is a pointer to our arguments structure. */
    struct arguments *arguments = state->input;
    int tmpArgValue;
    switch (key) {
        case 'o':
            arguments->outfile = arg;
            break;
        case 'm':
            tmpArgValue = atoi(arg);
            if (tmpArgValue>=0 && tmpArgValue <= UINT16_MAX) {
               arguments->mismatch = tmpArgValue;
            } else {
               errx(2,"-m \"%s\"=%i is not between [0,%d] !",arg,tmpArgValue,UINT16_MAX);
            }
            break;
        case 's':
            tmpArgValue = atoi(arg);
            if (tmpArgValue>=0 && tmpArgValue <= UINT16_MAX) {
               arguments->seedLen = tmpArgValue;
            } else {
               errx(2,"-s \"%s\"=%i is not between [0,%d] !",arg,tmpArgValue,UINT16_MAX);
            }
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

static inline uint_fast16_t compseq(char const* strA, char const* strB, uint_fast16_t maxMismatch) {
    uint_fast16_t mismatch=0;
    #ifdef DEBUG
    printf("[%s]-[%s]:\n",strA,strB);
    #endif
    while (*strA && *strB) {
        #ifdef DEBUG
        printf(" %.1s,%.1s:",strA,strB);
        #endif
        mismatch += (*strA++ != *strB++);
        #ifdef DEBUG
        printf("%d ",(int)mismatch);
        #endif
        if (mismatch > maxMismatch) {
            #ifdef DEBUG
            puts("");
            #endif
            return 0;
        }
    }
    #ifdef DEBUG
    puts("");
    #endif
    return mismatch;
}

struct prb_table *MainKmerPosTree, *KmerCacheTree;

void insertSeeds(char const* seq, struct prb_table * tree, uint_fast16_t SeedLen) {
    size_t ReadLength = strlen(seq);
    #ifdef DEBUG
    if (ReadLength<SeedLen) {
        fprintf(stderr,"[!]Read_Length(%d) < Seed_Length(%d) !\n",
            (unsigned int)ReadLength,(unsigned int)SeedLen);
        exit(1);
    }
    #endif
    for(uint_fast16_t i=0;i<=ReadLength-SeedLen;++i){
        currentPos=i;
    }
}

int main (int argc, char **argv) {
    struct arguments arguments;
    arguments.args=calloc(sizeof(size_t),argc);
    
    // Default values.
    arguments.outfile  = "./out.stat";
    arguments.mismatch = 10;
    arguments.seedLen  = 7;
    char *const*line = arguments.args-1;    // so that we can use *(++line)
    
    // Parse our arguments; every option seen by parse_opt will be reflected in arguments.
    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    fprintf(stderr,"Arguments: max_mismatch=%u, seed_length=%u, out_file=%s\nInput_files:",
        (unsigned int)arguments.mismatch,(unsigned int)arguments.seedLen,arguments.outfile);
    while ( *(++line) ) {
        fprintf(stderr, " [%s]", *line);
    }
    fputs("\n", stderr);
    line = arguments.args-1;
    
    //uint8_t* pSeqMismatchArray = calloc(READSCOUNT_INIT,sizeof(uint8_t));
    char** pSeqArray = malloc(READSCOUNT_INIT*sizeof(size_t));
    size_t SeqArrayLength = READSCOUNT_INIT;
    uint64_t* pMismatchCount = calloc(arguments.mismatch+1,sizeof(uint64_t));
    
    uint64_t BasesCount=0;
    uint64_t ReadsCount=0;
    uint64_t MaxReadLength=0;
    int_fast8_t diffRL=0;
    int ReadsTooLong=0;
    int ReadsTooShort=0;

    fputs("\nParsing Sequence Files:\n", stderr);
    G_TIMER_START;

	kseq_t *kseq;
	int_fast8_t ReadType;

    while ( *(++line) ) {
        fprintf(stderr, " <%s>", *line);
        kseq = kseq_open(*line);
        uint64_t RealSize,RealOffset;
        size_t ReadLength;
        if (kseq) {
            RealSize = kseq->f->RealSize;
            fprintf(stderr," RealSize:%ld ...",RealSize);
        	while ( (ReadType = kseq_read(kseq)) > 0 ) {
            	ReadLength = kseq->seq.l;
            	RealOffset = kseq->f->RealOffset;
            	#ifdef DEBUG
	            	printf("\n-ID:[%s,%s] %zu %zu %zd %ld\nSeq:[%s]\nQ:[%s] *%zx,%d    \n",
            			kseq->name.s,kseq->comment.s,kseq->seq.l,kseq->qual.l,ReadLength,kseq->f->RealOffset,
				        kseq->seq.s,kseq->qual.s,(size_t)&(kseq->seq.s),ReadType);
			    #endif
                if (ReadLength>0) {
                    if (ReadLength>UINT16_MAX) {
                        ReadLength=UINT16_MAX;
                        ReadsTooLong = 1;
                    } else if (ReadLength<arguments.seedLen) {
                        ReadsTooShort = 1;
                        continue;
                    }
                    if (ReadLength>MaxReadLength) MaxReadLength=ReadLength;
                    if (ReadsCount+1>SeqArrayLength) {
                        pSeqArray = realloc(pSeqArray,(SeqArrayLength+READSCOUNT_INC)*sizeof(size_t));
                        //memset(pSeqSeqArray+SeqMismatchArrayLength,0,READSCOUNT_INC*sizeof(size_t));
                        SeqArrayLength += READSCOUNT_INC;
                    }
                    pSeqArray[ReadsCount]=malloc(ReadLength+1);
                    strncpy(pSeqArray[ReadsCount],kseq->seq.s,ReadLength);
                    pSeqArray[ReadsCount][ReadLength]='\0';
                    ++ReadsCount;
                    BasesCount += ReadLength;
                    if (MaxReadLength != ReadLength) diffRL=1;
                }
        	}
        	//printf(" lasttype:%d   lastpos:%ld ------",ReadType,kseq->f->RealOffset);
        } else {
            fputs(": File not found !\n", stderr);
            continue;
        }
        fputs("\b\b\b\b, done !\n", stderr);
	    kseq_destroy(kseq);
    }
    uint8_t* pSeqMismatchArray = calloc(ReadsCount,sizeof(uint8_t));

    fprintf(stderr,"[!]Total_Reads: %ld\tMean_Read_Length: %f\n",ReadsCount,(double)BasesCount/(double)ReadsCount);

    KmerCacheTree = prb_create(compare_fixed_strings, &arguments.seedLen, NULL);
    MainKmerPosTree = prb_create(compare_fixed_strings, &arguments.seedLen, NULL);
    for(uint_fast32_t i=0;i<ReadsCount;++i) {
        insertSeeds(pSeqArray[i],KmerCacheTree,arguments.seedLen);
    }

    //uint64_t CmpCount=0;  // CmpCount == ReadsCount*(ReadsCount-1)/2
    uint64_t MisCount=0;
    uint64_t MisChgCount=0;
    for (size_t i=0;i<ReadsCount;++i) {
        for (size_t j=i+1;j<ReadsCount;++j) {
            //++CmpCount;
            uint_fast16_t miscount=compseq(pSeqArray[i],pSeqArray[j],arguments.mismatch);
            if (miscount) {
                ++MisCount;
                if (pSeqMismatchArray[i]==0 || pSeqMismatchArray[i]>miscount) {
                    pSeqMismatchArray[i]=miscount;
                    ++MisChgCount;
                }
                if (pSeqMismatchArray[j]==0 || pSeqMismatchArray[j]>miscount) {
                    pSeqMismatchArray[j]=miscount;
                    ++MisChgCount;
                }
            }
        }
    }

    for (size_t i=0;i<ReadsCount;++i) {
        #ifdef DEBUG
        printf("%zd/%zd\t[%s] %d\n",i,ReadsCount,pSeqArray[i],(int)pSeqMismatchArray[i]);
        #endif
        ++pMismatchCount[pSeqMismatchArray[i]];
    }

    fprintf(stderr,"[!]CmpTimes:%ld, mismatch_OK:%ld, mismatch_Update:%ld\n",
        ReadsCount*(ReadsCount-1)/2,MisCount,MisChgCount);
    fputs("\nCount Done!\n\nMismatch\tReadsCount\n", stderr);
    if (diffRL) fprintf(stderr,"[!]Different Read Length Found. Mismatch may be in correct.\n");

    FILE *fp = fopen(arguments.outfile, "w");
    fprintf(fp,"#Total_Reads: %ld\n#Mean_Read_Length: %f\n#Max_mismatch: %d\n\n#Mismatch\tReadsCount\n",
        ReadsCount,(double)BasesCount/(double)ReadsCount,(int)arguments.mismatch);
    for (uint_fast16_t i=0;i<=arguments.mismatch;++i) {
        if (pMismatchCount[i]) {
            fprintf(fp,"%d\t%ld\n",(int)i,pMismatchCount[i]);
            fprintf(stderr,"%d\t%ld\n",(int)i,pMismatchCount[i]);
        }
    }
    fprintf(fp,"\n# Mismatch=0 means unique.\n");
    if (ReadsTooLong) {
        fputs("#There are reads longer than "CAT(UINT16_MAX)".\n", fp);
        fputs("[!]There are reads longer than "CAT(UINT16_MAX)".\n", stderr);
    }
    if (ReadsTooShort) {
        fprintf(fp, "#There are reads shorter than Seed_Length:%u.\n", (unsigned int)arguments.seedLen);
        fprintf(stderr, "#There are reads shorter than Seed_Length:%u.\n", (unsigned int)arguments.seedLen);
    }
    fclose(fp);

    for (size_t i=0;i<ReadsCount;++i) {
        free(pSeqArray[i]);
    }
    free(pSeqArray);
    free(pSeqMismatchArray);
    free(arguments.args);

    G_TIMER_END;
    G_TIMER_PRINT;

    exit(EXIT_SUCCESS);
}

