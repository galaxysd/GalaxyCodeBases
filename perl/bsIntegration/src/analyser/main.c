#define _GNU_SOURCE
#include <err.h>
#include <argp.h>
//#include <bam/sam.h>
#include "klib/khash.h"
#include "klib/kvec.h"
//#include "uthash/utarray.h"
#include "ini.h"
#include "getch.h"
#include "timer.h"
#include "functions.h"

kvec_t(char*) aRefChrIDs, aVirusChrIDs;

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
//	UT_array *deplst;
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
			   errx(2,"[x]-%c \"%s\"=%i is not a integer of [1,%u] !",key,arg,tmpArgValue,UINT16_MAX-1);
			}
			break;
		case 'p':
			if ( arg[0] == 'g' || arg[0] == 'G' ) {
				arguments->programme = "grep";
			} else if ( arg[0] == 'a' || arg[0] == 'A' ) {
				arguments->programme = "analyse";
			} else {
				errx(2,"[x]-%c \"%s\" must be either \"grep\" or \"analyse\" !",key,arg);
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
			if (arguments->programme == NULL)
				errx(2,"[x]-p must be specified as either \"grep\" or \"analyse\" !");
			if (state->arg_num != 1) {
				errx(2,"[x]There can be only one input file, found [%d] !",state->arg_num);
				//argp_usage (state);
				argp_state_help (state,stderr,ARGP_HELP_STD_HELP);
			} else
				arguments->infile = arguments->args[0];
			break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

#ifdef DEBUGa
static int dumper(void* user, const char* section, const char* name, const char* value) {
    static char prev_section[50] = "";
    if (strcmp(section, prev_section)) {
        printf("%s[%s]\n", (prev_section[0] ? "\n" : ""), section);
        strncpy(prev_section, section, sizeof(prev_section));
        prev_section[sizeof(prev_section) - 1] = '\0';
    }
    printf("%s = %s\n", name, value);
    return 1;
}
#endif

static int ReadGrepINI(void* user, const char* section, const char* name, const char* value) {
	khiter_t ki;
	ChrInfo_t tmp; tmp.isHum = 9;	// We'll check whether T/F cover all items later, baka ⑨.
	BamInfo_t tbam, *pbam;
	int absent;
	char *word, *strtok_lasts;
	char *sep = ", ;";
	//#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
	if (strcmp(section, "Ref") == 0) {
		if (strcmp(name, "RefChrIDs") == 0) {
			while((word = strsep((char**)&value, sep)) != NULL) {
				if (*word) {
#ifdef DEBUGa
					printf("-[%s]-",word);
#endif
					kv_push(char*, aRefChrIDs, strdup(word));
				}
			}
			/*
			putchar('\n');
			for (size_t i=0; i<kv_size(aRefChrIDs);++i) {
				printf("=[%s]=",kv_A(aRefChrIDs, i));
			}
			putchar('\n');
			//printf("<%s>\n",value);
			*/
		} else if (strcmp(name, "VirusChrIDs") == 0) {
			for (word = strtok_r( (char*)value, sep, &strtok_lasts);
				word;	// We need to alter this 'const char *' of `value`.
				word = strtok_r(NULL, sep, &strtok_lasts))
				kv_push(char*, aVirusChrIDs, strdup(word));
		} else if (strcmp(name, "Refilename") == 0) {
			myConfig.RefileName = strdup(value);
		} else {
			tmp.ChrLen = atol(value);
			ki = kh_put(chrNFO, chrNFOp, name, &absent);
 			if (absent) kh_key(chrNFOp, ki) = strdup(name);
			kh_value(chrNFOp, ki) = tmp;
	   	}
	} else if (strcmp(section, "BamFiles") == 0) {
		ki = kh_put(bamNFO, bamNFOp, name, &absent);
		if (absent) {
			kh_key(bamNFOp, ki) = strdup(name);
			tbam.fileName = strdup(value);
			kh_value(bamNFOp, ki) = tbam;
		} else {
			pbam = &kh_value(bamNFOp, ki);
			pbam->fileName = strdup(value);
		}
	} else if (strcmp(section, "InsertSizes") == 0) {
		char * id = strdup(name);
		BamInfo_t *pt;
		size_t offset;
		char *ptr = strstr(id,".SD");
		if (ptr) {
			*ptr = '\0';
			offset = offsetof(BamInfo_t, SD);
		} else {	// http://stackoverflow.com/questions/37757902/in-c-how-to-point-to-different-members-with-offset-at-run-time
			offset = offsetof(BamInfo_t, insertSize);
		}
		ki = kh_put(bamNFO, bamNFOp, id, &absent);
		if (absent) {
			kh_key(bamNFOp, ki) = id;
			pt = &tbam;
			*(uint16_t*)((char*)pt+offset) = atol(value);	// http://stackoverflow.com/questions/3523145/pointer-arithmetic-for-void-pointer-in-c
			kh_value(bamNFOp, ki) = tbam;
		} else {
			free(id);
			pt = &kh_value(bamNFOp, ki);
			*(uint16_t*)((char*)pt+offset) = atol(value);	// thus use (char*) instead of (void*)
		}
	} else if (strcmp(section, "Output") == 0) {
		if (strcmp(name, "ProjectID") == 0)
			myConfig.ProjectID = strdup(value);
		else if (strcmp(name, "WorkDir") == 0)
			myConfig.WorkDir = strdup(value);
	} else {
		return 0;
	}
	return 1;
}

int main (int argc, char **argv) {
	struct arguments arguments;
	arguments.args=calloc(sizeof(size_t),argc);

	// Default values.
	arguments.isSAM = 0;
	arguments.interactive = 0;
	arguments.overlap = 25;
	myConfig.minGrepSlen = 5;

	// Parse our arguments; every option seen by parse_opt will be reflected in arguments.
	argp_parse (&argp, argc, argv, 0, 0, &arguments);

	printf("[!]Runing %s [%s] on [%s].\n",argv[0],arguments.programme,arguments.infile);
	if (arguments.interactive) {
		pressAnyKey();
	}

	G_TIMER_START;
	free(arguments.args);

	chrNFOp = kh_init(chrNFO);
	bamNFOp = kh_init(bamNFO);
	kv_init(aRefChrIDs);
	kv_init(aVirusChrIDs);
	int error;
#ifdef DEBUGa
	error = ini_parse(arguments.infile, dumper, NULL);
#endif
	error = ini_parse(arguments.infile, ReadGrepINI, NULL);
	if (error < 0) {
	    printf("[x]Can't read '%s'!\n", arguments.infile);
	    return 2;
	}
	else if (error) {
	    printf("[x]Bad config file (first error on line %d)!\n", error);
	    return 3;
	}
#ifdef DEBUGa
	printf("------------\n");
#endif
	ChrInfo_t * tmp;
	BamInfo_t * pbam;
	kh_cstr_t ChrID, BamID;
	khiter_t ki;
	for (size_t i=0; i<kv_size(aRefChrIDs);++i) {
		ChrID = kv_A(aRefChrIDs, i);
		ki = kh_get(chrNFO, chrNFOp, ChrID);
		if (ki == kh_end(chrNFOp)) {
			errx(3,"Cannot find Human ChrLen for [%s] !",ChrID);
		}
		tmp = &kh_value(chrNFOp, ki);
		tmp->isHum = true;
		free((char*)ChrID);
	}
	for (size_t i=0; i<kv_size(aVirusChrIDs);++i) {
		ChrID = kv_A(aVirusChrIDs, i);
		ki = kh_get(chrNFO, chrNFOp, ChrID);
		if (ki == kh_end(chrNFOp)) {
			errx(3,"[x]Cannot find Virus ChrLen for [%s] !",ChrID);
		}
		tmp = &kh_value(chrNFOp, ki);
		tmp->isHum = false;
		free((char*)ChrID);
	}
	kv_destroy(aRefChrIDs);
	kv_destroy(aVirusChrIDs);
	for (ki = kh_begin(chrNFOp); ki != kh_end(chrNFOp); ++ki)
		if (kh_exist(chrNFOp, ki)) {
			tmp = &kh_value(chrNFOp, ki);
			ChrID = kh_key(chrNFOp, ki);
			if (tmp->isHum == 9) {
				warnx("[!]Extra ChrLen exists [%s] !",ChrID);
			}
		}
	if (strcmp(arguments.programme,"grep")==0) {
		error = do_grep();
	} else if (strcmp(arguments.programme,"analyse")==0) {
		error = do_analyse();
	}
#ifdef DEBUGa
	printf("[!]ProjectID:[%s], WorkDir:[%s]\nRefileName:[%s]\n",myConfig.ProjectID,myConfig.WorkDir,myConfig.RefileName);
#endif
	for (ki = kh_begin(chrNFOp); ki != kh_end(chrNFOp); ++ki) {
		if (kh_exist(chrNFOp, ki)) {
			ChrID = kh_key(chrNFOp, ki);
#ifdef DEBUGa
			tmp = &kh_value(chrNFOp, ki);
			printf("%u [%s]=%d %u\n",ki,ChrID,tmp->ChrLen,tmp->isHum);
#endif
			free((char*)ChrID);
		}
	}
	for (ki = kh_begin(bamNFOp); ki != kh_end(bamNFOp); ++ki) {
		if (kh_exist(bamNFOp, ki)) {
			BamID = kh_key(bamNFOp, ki);
			pbam = &kh_value(bamNFOp, ki);
#ifdef DEBUGa
			printf("%u [%s]=%s\t%u %u\n",ki,BamID,pbam->fileName,pbam->insertSize,pbam->SD);
#endif
			free((char*)BamID);
			free((char*)pbam->fileName);
		}
	}
	kh_destroy(chrNFO, chrNFOp);
	kh_destroy(bamNFO, bamNFOp);
	free((char*)myConfig.ProjectID);
	free((char*)myConfig.WorkDir);
	free((char*)myConfig.RefileName);	// not const anymore
	G_TIMER_END;
	G_TIMER_PRINT;

	exit(EXIT_SUCCESS);
}
