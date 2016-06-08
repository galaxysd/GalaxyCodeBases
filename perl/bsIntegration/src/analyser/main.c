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
#include "getch.h"
#include "timer.h"

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

int UT_array_intsort(const void *a,const void*b) {
    int _a = *(int*)a;
    int _b = *(int*)b;
    return _a - _b;
}
//struct ChrData_hash_struct *ChrData = NULL;    /* important! initialize to NULL */

/* A description of the arguments we accept. */
static char args_doc[] = "input_files (SAM or BAM)";

/* The options we understand. */
static struct argp_option options[] = {
    {"interactive",      'i', 0,        0,  "pause before procedure"},
    {"sam",              's', 0,        0,  "input is SAM" },
    {"overlap_length",   'l', "25",     0,  "min overlap to connect reads" },
    {"min_depths",       'd', "1,2,3,5",0,  "list of min acceptable depth" },
    {"out_prefix",       'o', "./ctg",  0,  "prefix of output files" },
    { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments {
    uint_fast8_t isSAM, interactive;
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
    //utarray_new(arguments->deplst,&ut_int_icd);
    int tmpArgValue;
    switch (key) {
        case 'l':
            tmpArgValue = atoi(arg);
            if (tmpArgValue>=1 && tmpArgValue < UINT16_MAX) {   // k=l+1, so not == UINT16_MAX
               arguments->overlap = tmpArgValue;
            } else {
               errx(2,"-%c \"%s\"=%i is not a integer of [1,%u] !",key,arg,tmpArgValue,UINT16_MAX-1);
            }
            break;
        case 'd':
            ;
            char *tmpStr, *tmpStrN;
            if ( (tmpStrN = tmpStr = malloc(strlen(arg)+2)) == NULL )   // +1 for '\0' and +1 for extra ','
                err(3,NULL);
            //strcpy(tmpStr,arg);
            //printf("[%s] %s\n",tmpStr,arg);
            utarray_init(arguments->deplst,&ut_int_icd);
            for ( strcpy(tmpStrN,arg);
                   ;
                   tmpStrN = NULL) {
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
            utarray_sort(arguments->deplst,&UT_array_intsort);
            //free(tmpStr);
            *tmpStr = '\0';
            if ( (tmpStrN = malloc(strlen(arg)+2)) == NULL ) err(3,NULL);
            for ( int *p=(int*)utarray_front(arguments->deplst);
                  p!=NULL;
                  p=(int*)utarray_next(arguments->deplst,p) ) {
                sprintf(tmpStrN,"%i,",*p);
                strcat(tmpStr,tmpStrN);
            }
            free(tmpStrN);
            arguments->deplstStr = tmpStr;
            tmpStr += strlen(tmpStr)-1;
            *tmpStr = '\0';
            //printf("[%s] %s\n",arguments->deplstStr,arg);
            break;
        case 'o':
            arguments->outfile = arg;
            break;
        case 's':
            arguments->isSAM = 1;
            break;
        case 'i':
            arguments->interactive = 1;
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
    arguments.isSAM = 0;
    arguments.interactive = 0;
    arguments.deplstStr = "1,2,3,5";
    int tmpArgValue;
    utarray_new(arguments.deplst,&ut_int_icd);
    tmpArgValue = 1; utarray_push_back(arguments.deplst,&tmpArgValue);
    tmpArgValue = 2; utarray_push_back(arguments.deplst,&tmpArgValue);
    tmpArgValue = 3; utarray_push_back(arguments.deplst,&tmpArgValue);
    tmpArgValue = 5; utarray_push_back(arguments.deplst,&tmpArgValue);
    arguments.overlap = 25;
    arguments.outfile = "./ctg";

    // Parse our arguments; every option seen by parse_opt will be reflected in arguments.
    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    printf("Options:\noverlap_length = %u,\nmin_depths = [%s],\nout_prefix = [%s]\nInput file(s) type = %s .\n",
        arguments.overlap,arguments.deplstStr,arguments.outfile,arguments.isSAM?"SAM":"BAM");
    if (arguments.interactive) {
      pressAnyKey();
    }

    char *const*line = arguments.args-1;    // so that we can use *(++line)

    G_TIMER_START;

    samfile_t *samfp;
    fputs("\nParsing SAM/BAM Files:\n", stderr);
    while ( *(++line) ) {
        fprintf(stderr, " <%s> ... ", *line);
        if (arguments.isSAM) {
            samfp = samopen(*line, "r", 0);
        } else {
            samfp = samopen(*line, "rb", 0);
        }
        if (samfp) {
            bam_header_t *samhead = samfp->header;
            if ( samhead == NULL ) errx(3,"File Header Error.");
            //ChrDat_init(samhead);

            bam1_t *balignd = bam_init1();
            //while (samread(samfp, balignd) >= 0) do_stat(balignd, arguments.overlap, &Data);
            bam_destroy1(balignd);

            fputs("done !\n", stderr);
            samclose(samfp);
            continue;
        }
        warn("failed to open file");
        //fputs("failed !\n", stderr);
    }
    free(arguments.args);
    fputs("\nReading Done!\n\nBegin Stat:\n", stderr);

    for ( int *p=(int*)utarray_front(arguments.deplst);
          p!=NULL;
          p=(int*)utarray_next(arguments.deplst,p) ) {
        fprintf(stderr,"minDepth = %i: ",*p);
        //do_contig(*p, &Data, arguments.outfile);
        fputs("done !\n", stderr);
    }

    FILE *fp = fopen(arguments.outfile, "w");
    fprintf(fp,"#Main\n");
    fclose(fp);

    G_TIMER_END;
    G_TIMER_PRINT;

    exit(EXIT_SUCCESS);
}
