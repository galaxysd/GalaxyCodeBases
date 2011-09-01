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
	   ,
           arguments.args[0],
           arguments.outprefix
           );
      pressAnyKey();
    }

//printf("Out[%s][%s][%s]\n",outStat,outDat,outLog);

    ssize_t read;
    char *line;
    uint64_t insertedCount=0;

    fputs("SDLA nfo: ", stderr);


    G_TIMER_START;

    fputs("\nParsing Sequence Files:\n", stderr);

    G_TIMER_END;
    G_TIMER_PRINT;

    exit(EXIT_SUCCESS);
}
