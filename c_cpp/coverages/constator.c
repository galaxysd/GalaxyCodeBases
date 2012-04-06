#define _GNU_SOURCE
#include <stdlib.h> //EXIT_FAILURE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <err.h>
#include <argp.h>
#include <math.h>
#include <bam/sam.h>
#include "uthash/utarray.h"
#include "timer.h"

//#define MAXREADLEN (8ul*1024*1024)
//uint64_t ReadsLenArr[MAXREADLEN];

const char *argp_program_version =
    "Continuous Stator 0.1 (SAM/BAM) @"__TIME__ "," __DATE__;
const char *argp_program_bug_address =
    "<galaxy001@gmail.com>";

/* Program documentation. */
static char doc[] =
    "Continuous Stator"
#ifdef DEBUG
    " (Debug Version)"
#endif
#ifdef TEST
    " (Test Version)"
#endif
;

/* A description of the arguments we accept. */
static char args_doc[] = "input_files (SAM or BAM)";

/* The options we understand. */
static struct argp_option options[] = {
    {"interactive",      'i', 0,        0,  "Prompt arguments before procedure"},
    {"verbose",          'v', 0,        0,  "Produce verbose output" },
    {"overlap_length",   'l', "25",     0,  "min overlap to connect reads" },
    {"min_depths",       'd', "1,2,3,5",0,  "list of min acceptable depth" },
    {"out_prefix",       'o', "./ctg",  0,  "prefix of output files" },
    { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments {
    uint16_t overlap;
    char *deplstStr;
    UT_array *deplst;
    char **args;
    char *outfile;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state) {
/* Get the input argument from argp_parse, which we
  know is a pointer to our arguments structure. */
    struct arguments *arguments = state->input;
    utarray_new(arguments->deplst,&ut_int_icd);
    int tmpArgValue;
    switch (key) {
        case 'l':
            tmpArgValue = atoi(arg);
            if (tmpArgValue>=1 && tmpArgValue <= UINT16_MAX) {
               arguments->overlap = tmpArgValue;
            } else {
               errx(2,"-%c \"%s\"=%i is not a integer of [1,%u] !",key,arg,tmpArgValue,UINT16_MAX);
            }
            break;
        case 'd':
            arguments->deplstStr = arg;
            char *tmpStr, *tmpStrN;
            if ( (tmpStr = tmpStrN = malloc(strlen(arg)+1)) == NULL )
                err(3,NULL);
            //strcpy(tmpStr,arg);
            //printf("[%s] %s\n",tmpStr,arg);
            for ( strcpy(tmpStrN,arg); ; tmpStrN = NULL) {
                char *token = strtok(tmpStrN, ",; .");
                if (token == NULL) break;
                tmpArgValue = atoi(token);
                if (tmpArgValue>=1 && tmpArgValue <= UINT8_MAX) {
                    utarray_push_back(arguments->deplst,&tmpArgValue);
                    //printf("[%s] %s %d\n",tmpStr,token,tmpArgValue);
                } else {
                    errx(2,"In -%c \"%s\":[%i] is not a integer of [1,%u] !",key,arg,tmpArgValue,UINT8_MAX);
                }
            }
            free(tmpStr);
            break;
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
    uint64_t minReadLen=1;
    double SStd;

    char *const*line = arguments.args-1;    // so that we can use *(++line)

    G_TIMER_START;

    fputs("\nParsing Sequence Files:\n", stderr);
    while ( *(++line) ) {
        fprintf(stderr, " <%s> ...", *line);
        int *seqobj = NULL;
        if (seqobj) {
        } else continue;
        fputs("\b\b\b\b, done !\n", stderr);
    }
    free(arguments.args);
    fputs("\nCount Done!\n", stderr);

    FILE *fp = fopen(arguments.outfile, "w");
    fprintf(fp,"#Total_Bases: %lu\n#Total_Reads: %lu\n#Avg_Read_Len: %.1f\tStd: %.3f\n"
        "#Read_Len_Range: [%lu,%lu]\n"
        "#Overflow: %u\n"
        "\n#Read_Len\tCount\tRatio\n",
        allbases,allreads,
        0.05+((double)allbases/(double)allreads), SStd,
        minReadLen,maxReadLen,0);
    fclose(fp);

    G_TIMER_END;
    G_TIMER_PRINT;

    exit(EXIT_SUCCESS);
}

