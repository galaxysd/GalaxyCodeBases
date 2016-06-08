#define _GNU_SOURCE
#include <stdlib.h> //EXIT_FAILURE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <err.h>
#include <argp.h>
#include <math.h>
//#include <bam/sam.h>
#include "uthash/utarray.h"
#include "ini.h"
#include "getch.h"
#include "timer.h"

const char *argp_program_version =
	"bsuit Analyser 0.1 @"__TIME__ "," __DATE__;
const char *argp_program_bug_address =
	"<galaxy001@gmail.com>";

/* Program documentation. */
static char doc[] =
	"BS-viral-inte Analyser"
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
//struct ChrData_hash_struct *ChrData = NULL;	/* important! initialize to NULL */

/* A description of the arguments we accept. */
static char args_doc[] = "<Grep INI file>";

/* The options we understand. */
static struct argp_option options[] = {
	{"interactive",	  'i', 0,		0,  "pause before procedure"},
//	{"sam",			  's', 0,		0,  "input is SAM" },
//	{"overlap_length",   'l', "25",	 0,  "min overlap to connect reads" },
//	{"min_depths",	   'd', "1,2,3,5",0,  "list of min acceptable depth" },
	{"programme",		'p', "grep",	   0,  "grep / analyse" },
	{"verbose", 'v', 0, 0, "Produce verbose output"},
	{ 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments {
	uint_fast8_t isSAM, interactive;
	uint16_t overlap;
	char *deplstStr;
	UT_array *deplst;
	char **args;
	char const* programme;
	char const* infile;
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
		case 'p':
			if ( arg[0] == 'g' || arg[0] == 'G' ) {
				arguments->programme = "grep";
			} else if ( arg[0] == 'a' || arg[0] == 'A' ) {
				arguments->programme = "analyse";
			} else {
				errx(2,"-%c \"%s\" must be either \"grep\" or \"analyse\" !",key,arg);
			}
			//arguments->programme = arg;
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
			if (state->arg_num != 1) {
				errx(2,"There can be only one input file, found [%d] !",state->arg_num);
				//argp_usage (state);
				argp_state_help (state,stderr,ARGP_HELP_STD_HELP);
		   } else
				arguments->infile = arguments->args[0];
			if (arguments->programme == NULL)
				errx(2,"-p must be specified as either \"grep\" or \"analyse\" !");
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
	arguments.overlap = 25;

	// Parse our arguments; every option seen by parse_opt will be reflected in arguments.
	argp_parse (&argp, argc, argv, 0, 0, &arguments);

	printf("[!]Runing %s [%s] on [%s].\n",argv[0],arguments.programme,arguments.infile);
	if (arguments.interactive) {
		pressAnyKey();
	}

	G_TIMER_START;

	free(arguments.args);

	G_TIMER_END;
	G_TIMER_PRINT;

	exit(EXIT_SUCCESS);
}
