#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include "MurmurHash3.h"
#include <stdio.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

const char *argp_program_version =
"readscorr 0.1";
const char *argp_program_bug_address =
"<huxuesong@genomics.org.cn>";

/* Program documentation. */
static char doc[] =
"a Solexa Reads Corrector using Bloom filter"
#ifdef DEBUG
" (Debug Version)"
#endif
;

/* A description of the arguments we accept. */
static char args_doc[] = "input_list (of FASTA or FASTQ files, gzipped is OK)";

/* The options we understand. */
static struct argp_option options[] = {
{"verbose",  'v', 0,      0,  "Produce verbose output" },
{"quiet",    'q', 0,      0,  "Don't produce any output" },
{"silent",   's', 0,      OPTION_ALIAS },
{"outprefix",   'o', "./out", 0, "Output to [./out.{dat,stat,log}]" },
{ 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
char *args[1];                /* arg1 */
int silent, verbose;
char *outprefix;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
/* Get the input argument from argp_parse, which we
  know is a pointer to our arguments structure. */
struct arguments *arguments = state->input;

switch (key)
 {
 case 'q': case 's':
   arguments->silent = 1;
   break;
 case 'v':
   arguments->verbose = 1;
   break;
 case 'o':
   arguments->outprefix = arg;
   break;

 case ARGP_KEY_ARG:
   if (state->arg_num >= 1)
     /* Too many arguments. */
     argp_usage (state);

   arguments->args[state->arg_num] = arg;

   break;

 case ARGP_KEY_END:
   if (state->arg_num < 1)
     /* Not enough arguments. */
     argp_usage (state);
   break;

 default:
   return ARGP_ERR_UNKNOWN;
 }
return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

int main (int argc, char **argv)
{
uint64_t tout[2];
MurmurHash3_x64_128("test",4,123,tout);

struct arguments arguments;

/* Default values. */
arguments.silent = 0;
arguments.verbose = 0;
arguments.outprefix = "./out";

/* Parse our arguments; every option seen by parse_opt will
  be reflected in arguments. */
argp_parse (&argp, argc, argv, 0, 0, &arguments);

printf ("ARG1 = %s\nOutputPrefix = %s\n"
       "VERBOSE = %s\nSILENT = %s\n",
       arguments.args[0],
       arguments.outprefix,
       arguments.verbose ? "yes" : "no",
       arguments.silent ? "yes" : "no");

char * str=(char *) arguments.outprefix;
size_t len=strlen(str);
uint32_t key=atoi(arguments.args[0]);
MurmurHash3_x64_128(str,len,key,tout);
printf("Str:[%s] Seed:[%i] -> [%016lx %016lx]\n",str,key,tout[0],tout[1]);
//free(tout);

	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen(arguments.args[0], "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		printf("name: %s\n", seq->name.s);
		if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
		printf("seq: %s\n", seq->seq.s);
		if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
	}
	printf("return value: %d\n", l);
	kseq_destroy(seq);
	gzclose(fp);


exit (0);
}

