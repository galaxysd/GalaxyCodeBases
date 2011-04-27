#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <err.h>
#include <argp.h>
#include "MurmurHash3.h"
//#include <zlib.h>
//#include "kseq.h"
//KSEQ_INIT(gzFile, gzread)

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
    {"interactive", 'i', 0,   0,  "Prompt arguments before procedure"},
    {"verbose",  'v', 0,      0,  "Produce verbose output" },
    {"quiet",    'q', 0,      0,  "Don't produce any output" },
    //{"silent",   's', 0,      OPTION_ALIAS },
    {"outprefix",'o', "./out",0,  "Output to [./out.{dat,stat,log}]" },
    {"kmersize", 'k', "21",   0,  "K-mer size, must be odd number" },
    {"bloomsize",'b', "256",  0,  "Size in MiB for Bloom Filter"  },
    { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments {
    char *args[1];                /* arg1 */
    uint_fast8_t silent, verbose, interactive;  //_Bool is the same under amd64, as typedef unsigned char uint_fast8_t;
    int bloomsize, kmersize;
    char *outprefix;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state) {
/* Get the input argument from argp_parse, which we
  know is a pointer to our arguments structure. */
    struct arguments *arguments = state->input;
    
    int tmpArgValue;
    switch (key) {
        case 'q':// case 's':
            arguments->silent = 1;
            break;
        case 'v':
            arguments->verbose = 1;
            break;
        case 'i':
            arguments->interactive = 1;
            break;
        case 'o':
            arguments->outprefix = arg;
            break;
        case 'k':
            tmpArgValue = atoi(arg);
            if (tmpArgValue>2 && tmpArgValue%2) {    // odd numbers only
               arguments->kmersize = tmpArgValue;
            } else {
               errx(2,"-k \"%s\"=%i is not a positive odd number !",arg,tmpArgValue);
            }
            break;
        case 'b':
            tmpArgValue = atoi(arg);
            if (tmpArgValue>0)
               arguments->bloomsize = tmpArgValue;
            else
               errx(2,"-b \"%s\"=%i is not a positive number !",arg,tmpArgValue);
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

/*
int main (int argc, char **argv)                                          
{                                                                         
uint64_t tout[2];                                                         
MurmurHash3_x64_128("test",4,123,tout);                                   
                                                                          
struct arguments arguments;                                               
                                                                          
// Default values.              
arguments.silent = 0;                                                     
arguments.verbose = 0;                                                    
arguments.outprefix = "./out";                                            
                                                                          
// Parse our arguments; every option seen by parse_opt will be reflected in arguments.                           
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
*/

int main (int argc, char **argv) {
    struct arguments arguments;
    
    // Default values.
    arguments.silent = 0;
    arguments.verbose = 0;
    arguments.interactive = 0;
    arguments.outprefix = "./out";
    arguments.kmersize = 21;
    arguments.bloomsize = 256;
    
    // Parse our arguments; every option seen by parse_opt will be reflected in arguments.
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
    
    size_t bloomLen = 1024*1024*arguments.bloomsize;
    printf ("ARG1 = %s\nOutputPrefix = %s\n"
           "VERBOSE = %s\nSILENT = %s\n"
            "kmersize = %i, bloomsize = %i\n",
           arguments.args[0],
           arguments.outprefix,
           arguments.verbose ? "yes" : "no",
           arguments.silent ? "yes" : "no",
           arguments.kmersize,arguments.bloomsize);

	FILE *fp;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	fp = fopen(arguments.args[0], "r");
	if (fp == NULL) err(EXIT_FAILURE, "Cannot open input_list [%s]", arguments.args[0]); //exit(EXIT_FAILURE);

    fputs("\nFirst pass:\n", stderr);
	while ((read = getline(&line, &len, fp)) != -1) {
		if (*(line+read-1)=='\n') *(line+(--read))='\0';	// line[read-1] = '\0';
		fprintf(stderr, " <%s> ...", line);
		sleep(1);	// the Call ...
		fputs("\b\b\b\b, done !\n", stderr);
	}

    rewind(fp);
	fputs("\nSecond pass:\n", stderr);
	while ((read = getline(&line, &len, fp)) != -1) {
		if (*(line+read-1)=='\n') *(line+(--read))='\0';
		printf("[%zu]+%zu <%s>\n", read, len, line);
	}
	fputs("\nHashing Done!\n", stderr);
	free(line);

    exit(EXIT_SUCCESS);
}
