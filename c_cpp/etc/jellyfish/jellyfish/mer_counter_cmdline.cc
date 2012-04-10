/*
  File autogenerated by gengetopt version 2.22.4
  generated with the following command:
  gengetopt --show-required --default-option -c cc -H hpp -F mer_counter_cmdline -f mer_counter_cmdline -a mer_counter_args --unamed-opts=file.f[aq]

  The developers of gengetopt consider the fixed text that goes in all
  gengetopt output files to be in the public domain:
  we make no copyright claims on it.
*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

#include <getopt.h>

#include "mer_counter_cmdline.hpp"

const char *mer_counter_args_purpose = "Count k-mers or qmers in fasta or fastq files";

const char *mer_counter_args_usage = "Usage: jellyfish count [OPTIONS]... [file.f[aq]]...";

const char *mer_counter_args_description = "";

const char *mer_counter_args_full_help[] = {
  "  -h, --help                    Print help and exit",
  "      --full-help               Print help, including hidden options, and exit",
  "  -V, --version                 Print version and exit",
  "  -m, --mer-len=INT             Length of mer (mandatory)",
  "  -s, --size=LONG               Hash size (mandatory)",
  "  -t, --threads=INT             Number of threads  (default=`1')",
  "  -o, --output=STRING           Output prefix  (default=`mer_counts')",
  "  -c, --counter-len=Length in bits\n                                Length of counting field  (default=`7')",
  "      --out-counter-len=Length in bytes\n                                Length of counter field in output  \n                                  (default=`4')",
  "  -C, --both-strands            Count both strand, canonical representation  \n                                  (default=off)",
  "  -p, --reprobes=INT            Maximum number of reprobes  (default=`62')",
  "  -r, --raw                     Write raw database  (default=off)",
  "  -q, --quake                   Quake compatibility mode  (default=off)",
  "      --quality-start=INT       Starting ASCII for quality values  \n                                  (default=`64')",
  "      --min-quality=INT         Minimum quality. A base with lesser quality \n                                  becomes an N  (default=`0')",
  "  -L, --lower-count=LONG        Don't output k-mer with count < lower-count",
  "  -U, --upper-count=LONG        Don't output k-mer with count > upper-count",
  "      --matrix=Matrix file      Hash function binary matrix",
  "      --timing=Timing file      Print timing information",
  "  -w, --no-write                Don't write database  (default=off)",
  "  -u, --measure                 Write usage statistics  (default=off)",
  "      --buffers=LONG            Number of buffers per thread",
  "      --buffer-size=LONG        Size of buffers  (default=`8192')",
  "      --out-buffer-size=LONG    Size of output buffer per thread  \n                                  (default=`20000000')",
    0
};

static void
init_help_array(void)
{
  mer_counter_args_help[0] = mer_counter_args_full_help[0];
  mer_counter_args_help[1] = mer_counter_args_full_help[1];
  mer_counter_args_help[2] = mer_counter_args_full_help[2];
  mer_counter_args_help[3] = mer_counter_args_full_help[3];
  mer_counter_args_help[4] = mer_counter_args_full_help[4];
  mer_counter_args_help[5] = mer_counter_args_full_help[5];
  mer_counter_args_help[6] = mer_counter_args_full_help[6];
  mer_counter_args_help[7] = mer_counter_args_full_help[7];
  mer_counter_args_help[8] = mer_counter_args_full_help[8];
  mer_counter_args_help[9] = mer_counter_args_full_help[9];
  mer_counter_args_help[10] = mer_counter_args_full_help[10];
  mer_counter_args_help[11] = mer_counter_args_full_help[11];
  mer_counter_args_help[12] = mer_counter_args_full_help[12];
  mer_counter_args_help[13] = mer_counter_args_full_help[13];
  mer_counter_args_help[14] = mer_counter_args_full_help[14];
  mer_counter_args_help[15] = mer_counter_args_full_help[15];
  mer_counter_args_help[16] = mer_counter_args_full_help[16];
  mer_counter_args_help[17] = mer_counter_args_full_help[17];
  mer_counter_args_help[18] = mer_counter_args_full_help[18];
  mer_counter_args_help[19] = 0; 
  
}

const char *mer_counter_args_help[20];

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
  , ARG_LONG
} mer_counter_cmdline_arg_type;

static
void clear_given (struct mer_counter_args *args_info);
static
void clear_args (struct mer_counter_args *args_info);

static int
mer_counter_cmdline_internal (int argc, char **argv, struct mer_counter_args *args_info,
                        struct mer_counter_cmdline_params *params, const char *additional_error);

static int
mer_counter_cmdline_required2 (struct mer_counter_args *args_info, const char *prog_name, const char *additional_error);

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct mer_counter_args *args_info)
{
  args_info->help_given = 0 ;
  args_info->full_help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->mer_len_given = 0 ;
  args_info->size_given = 0 ;
  args_info->threads_given = 0 ;
  args_info->output_given = 0 ;
  args_info->counter_len_given = 0 ;
  args_info->out_counter_len_given = 0 ;
  args_info->both_strands_given = 0 ;
  args_info->reprobes_given = 0 ;
  args_info->raw_given = 0 ;
  args_info->quake_given = 0 ;
  args_info->quality_start_given = 0 ;
  args_info->min_quality_given = 0 ;
  args_info->lower_count_given = 0 ;
  args_info->upper_count_given = 0 ;
  args_info->matrix_given = 0 ;
  args_info->timing_given = 0 ;
  args_info->no_write_given = 0 ;
  args_info->measure_given = 0 ;
  args_info->buffers_given = 0 ;
  args_info->buffer_size_given = 0 ;
  args_info->out_buffer_size_given = 0 ;
}

static
void clear_args (struct mer_counter_args *args_info)
{
  FIX_UNUSED (args_info);
  args_info->mer_len_orig = NULL;
  args_info->size_orig = NULL;
  args_info->threads_arg = 1;
  args_info->threads_orig = NULL;
  args_info->output_arg = gengetopt_strdup ("mer_counts");
  args_info->output_orig = NULL;
  args_info->counter_len_arg = 7;
  args_info->counter_len_orig = NULL;
  args_info->out_counter_len_arg = 4;
  args_info->out_counter_len_orig = NULL;
  args_info->both_strands_flag = 0;
  args_info->reprobes_arg = 62;
  args_info->reprobes_orig = NULL;
  args_info->raw_flag = 0;
  args_info->quake_flag = 0;
  args_info->quality_start_arg = 64;
  args_info->quality_start_orig = NULL;
  args_info->min_quality_arg = 0;
  args_info->min_quality_orig = NULL;
  args_info->lower_count_orig = NULL;
  args_info->upper_count_orig = NULL;
  args_info->matrix_arg = NULL;
  args_info->matrix_orig = NULL;
  args_info->timing_arg = NULL;
  args_info->timing_orig = NULL;
  args_info->no_write_flag = 0;
  args_info->measure_flag = 0;
  args_info->buffers_orig = NULL;
  args_info->buffer_size_arg = 8192;
  args_info->buffer_size_orig = NULL;
  args_info->out_buffer_size_arg = 20000000;
  args_info->out_buffer_size_orig = NULL;
  
}

static
void init_args_info(struct mer_counter_args *args_info)
{

  init_help_array(); 
  args_info->help_help = mer_counter_args_full_help[0] ;
  args_info->full_help_help = mer_counter_args_full_help[1] ;
  args_info->version_help = mer_counter_args_full_help[2] ;
  args_info->mer_len_help = mer_counter_args_full_help[3] ;
  args_info->size_help = mer_counter_args_full_help[4] ;
  args_info->threads_help = mer_counter_args_full_help[5] ;
  args_info->output_help = mer_counter_args_full_help[6] ;
  args_info->counter_len_help = mer_counter_args_full_help[7] ;
  args_info->out_counter_len_help = mer_counter_args_full_help[8] ;
  args_info->both_strands_help = mer_counter_args_full_help[9] ;
  args_info->reprobes_help = mer_counter_args_full_help[10] ;
  args_info->raw_help = mer_counter_args_full_help[11] ;
  args_info->quake_help = mer_counter_args_full_help[12] ;
  args_info->quality_start_help = mer_counter_args_full_help[13] ;
  args_info->min_quality_help = mer_counter_args_full_help[14] ;
  args_info->lower_count_help = mer_counter_args_full_help[15] ;
  args_info->upper_count_help = mer_counter_args_full_help[16] ;
  args_info->matrix_help = mer_counter_args_full_help[17] ;
  args_info->timing_help = mer_counter_args_full_help[18] ;
  args_info->no_write_help = mer_counter_args_full_help[19] ;
  args_info->measure_help = mer_counter_args_full_help[20] ;
  args_info->buffers_help = mer_counter_args_full_help[21] ;
  args_info->buffer_size_help = mer_counter_args_full_help[22] ;
  args_info->out_buffer_size_help = mer_counter_args_full_help[23] ;
  
}

void
mer_counter_cmdline_print_version (void)
{
  printf ("%s %s\n",
     (strlen(MER_COUNTER_CMDLINE_PACKAGE_NAME) ? MER_COUNTER_CMDLINE_PACKAGE_NAME : MER_COUNTER_CMDLINE_PACKAGE),
     MER_COUNTER_CMDLINE_VERSION);
}

static void print_help_common(void) {
  mer_counter_cmdline_print_version ();

  if (strlen(mer_counter_args_purpose) > 0)
    printf("\n%s\n", mer_counter_args_purpose);

  if (strlen(mer_counter_args_usage) > 0)
    printf("\n%s\n", mer_counter_args_usage);

  printf("\n");

  if (strlen(mer_counter_args_description) > 0)
    printf("%s\n\n", mer_counter_args_description);
}

void
mer_counter_cmdline_print_help (void)
{
  int i = 0;
  print_help_common();
  while (mer_counter_args_help[i])
    printf("%s\n", mer_counter_args_help[i++]);
}

void
mer_counter_cmdline_print_full_help (void)
{
  int i = 0;
  print_help_common();
  while (mer_counter_args_full_help[i])
    printf("%s\n", mer_counter_args_full_help[i++]);
}

void
mer_counter_cmdline_init (struct mer_counter_args *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);

  args_info->inputs = 0;
  args_info->inputs_num = 0;
}

void
mer_counter_cmdline_params_init(struct mer_counter_cmdline_params *params)
{
  if (params)
    { 
      params->override = 0;
      params->initialize = 1;
      params->check_required = 1;
      params->check_ambiguity = 0;
      params->print_errors = 1;
    }
}

struct mer_counter_cmdline_params *
mer_counter_cmdline_params_create(void)
{
  struct mer_counter_cmdline_params *params = 
    (struct mer_counter_cmdline_params *)malloc(sizeof(struct mer_counter_cmdline_params));
  mer_counter_cmdline_params_init(params);  
  return params;
}

static void
free_string_field (char **s)
{
  if (*s)
    {
      free (*s);
      *s = 0;
    }
}


static void
mer_counter_cmdline_release (struct mer_counter_args *args_info)
{
  unsigned int i;
  free_string_field (&(args_info->mer_len_orig));
  free_string_field (&(args_info->size_orig));
  free_string_field (&(args_info->threads_orig));
  free_string_field (&(args_info->output_arg));
  free_string_field (&(args_info->output_orig));
  free_string_field (&(args_info->counter_len_orig));
  free_string_field (&(args_info->out_counter_len_orig));
  free_string_field (&(args_info->reprobes_orig));
  free_string_field (&(args_info->quality_start_orig));
  free_string_field (&(args_info->min_quality_orig));
  free_string_field (&(args_info->lower_count_orig));
  free_string_field (&(args_info->upper_count_orig));
  free_string_field (&(args_info->matrix_arg));
  free_string_field (&(args_info->matrix_orig));
  free_string_field (&(args_info->timing_arg));
  free_string_field (&(args_info->timing_orig));
  free_string_field (&(args_info->buffers_orig));
  free_string_field (&(args_info->buffer_size_orig));
  free_string_field (&(args_info->out_buffer_size_orig));
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, const char *values[])
{
  FIX_UNUSED (values);
  if (arg) {
    fprintf(outfile, "%s=\"%s\"\n", opt, arg);
  } else {
    fprintf(outfile, "%s\n", opt);
  }
}


int
mer_counter_cmdline_dump(FILE *outfile, struct mer_counter_args *args_info)
{
  int i = 0;

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot dump options to stream\n", MER_COUNTER_CMDLINE_PACKAGE);
      return EXIT_FAILURE;
    }

  if (args_info->help_given)
    write_into_file(outfile, "help", 0, 0 );
  if (args_info->full_help_given)
    write_into_file(outfile, "full-help", 0, 0 );
  if (args_info->version_given)
    write_into_file(outfile, "version", 0, 0 );
  if (args_info->mer_len_given)
    write_into_file(outfile, "mer-len", args_info->mer_len_orig, 0);
  if (args_info->size_given)
    write_into_file(outfile, "size", args_info->size_orig, 0);
  if (args_info->threads_given)
    write_into_file(outfile, "threads", args_info->threads_orig, 0);
  if (args_info->output_given)
    write_into_file(outfile, "output", args_info->output_orig, 0);
  if (args_info->counter_len_given)
    write_into_file(outfile, "counter-len", args_info->counter_len_orig, 0);
  if (args_info->out_counter_len_given)
    write_into_file(outfile, "out-counter-len", args_info->out_counter_len_orig, 0);
  if (args_info->both_strands_given)
    write_into_file(outfile, "both-strands", 0, 0 );
  if (args_info->reprobes_given)
    write_into_file(outfile, "reprobes", args_info->reprobes_orig, 0);
  if (args_info->raw_given)
    write_into_file(outfile, "raw", 0, 0 );
  if (args_info->quake_given)
    write_into_file(outfile, "quake", 0, 0 );
  if (args_info->quality_start_given)
    write_into_file(outfile, "quality-start", args_info->quality_start_orig, 0);
  if (args_info->min_quality_given)
    write_into_file(outfile, "min-quality", args_info->min_quality_orig, 0);
  if (args_info->lower_count_given)
    write_into_file(outfile, "lower-count", args_info->lower_count_orig, 0);
  if (args_info->upper_count_given)
    write_into_file(outfile, "upper-count", args_info->upper_count_orig, 0);
  if (args_info->matrix_given)
    write_into_file(outfile, "matrix", args_info->matrix_orig, 0);
  if (args_info->timing_given)
    write_into_file(outfile, "timing", args_info->timing_orig, 0);
  if (args_info->no_write_given)
    write_into_file(outfile, "no-write", 0, 0 );
  if (args_info->measure_given)
    write_into_file(outfile, "measure", 0, 0 );
  if (args_info->buffers_given)
    write_into_file(outfile, "buffers", args_info->buffers_orig, 0);
  if (args_info->buffer_size_given)
    write_into_file(outfile, "buffer-size", args_info->buffer_size_orig, 0);
  if (args_info->out_buffer_size_given)
    write_into_file(outfile, "out-buffer-size", args_info->out_buffer_size_orig, 0);
  

  i = EXIT_SUCCESS;
  return i;
}

int
mer_counter_cmdline_file_save(const char *filename, struct mer_counter_args *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", MER_COUNTER_CMDLINE_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  i = mer_counter_cmdline_dump(outfile, args_info);
  fclose (outfile);

  return i;
}

void
mer_counter_cmdline_free (struct mer_counter_args *args_info)
{
  mer_counter_cmdline_release (args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = 0;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
mer_counter_cmdline (int argc, char **argv, struct mer_counter_args *args_info)
{
  return mer_counter_cmdline2 (argc, argv, args_info, 0, 1, 1);
}

int
mer_counter_cmdline_ext (int argc, char **argv, struct mer_counter_args *args_info,
                   struct mer_counter_cmdline_params *params)
{
  int result;
  result = mer_counter_cmdline_internal (argc, argv, args_info, params, 0);

  if (result == EXIT_FAILURE)
    {
      mer_counter_cmdline_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
mer_counter_cmdline2 (int argc, char **argv, struct mer_counter_args *args_info, int override, int initialize, int check_required)
{
  int result;
  struct mer_counter_cmdline_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = mer_counter_cmdline_internal (argc, argv, args_info, &params, 0);

  if (result == EXIT_FAILURE)
    {
      mer_counter_cmdline_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
mer_counter_cmdline_required (struct mer_counter_args *args_info, const char *prog_name)
{
  int result = EXIT_SUCCESS;

  if (mer_counter_cmdline_required2(args_info, prog_name, 0) > 0)
    result = EXIT_FAILURE;

  if (result == EXIT_FAILURE)
    {
      mer_counter_cmdline_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
mer_counter_cmdline_required2 (struct mer_counter_args *args_info, const char *prog_name, const char *additional_error)
{
  int error = 0;
  FIX_UNUSED (additional_error);

  /* checks for required options */
  if (! args_info->mer_len_given)
    {
      fprintf (stderr, "%s: '--mer-len' ('-m') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->size_given)
    {
      fprintf (stderr, "%s: '--size' ('-s') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  
  /* checks for dependences among options */

  return error;
}


static char *package_name = 0;

/**
 * @brief updates an option
 * @param field the generic pointer to the field to update
 * @param orig_field the pointer to the orig field
 * @param field_given the pointer to the number of occurrence of this option
 * @param prev_given the pointer to the number of occurrence already seen
 * @param value the argument for this option (if null no arg was specified)
 * @param possible_values the possible values for this option (if specified)
 * @param default_value the default value (in case the option only accepts fixed values)
 * @param arg_type the type of this option
 * @param check_ambiguity @see mer_counter_cmdline_params.check_ambiguity
 * @param override @see mer_counter_cmdline_params.override
 * @param no_free whether to free a possible previous value
 * @param multiple_option whether this is a multiple option
 * @param long_opt the corresponding long option
 * @param short_opt the corresponding short option (or '-' if none)
 * @param additional_error possible further error specification
 */
static
int update_arg(void *field, char **orig_field,
               unsigned int *field_given, unsigned int *prev_given, 
               char *value, const char *possible_values[],
               const char *default_value,
               mer_counter_cmdline_arg_type arg_type,
               int check_ambiguity, int override,
               int no_free, int multiple_option,
               const char *long_opt, char short_opt,
               const char *additional_error)
{
  char *stop_char = 0;
  const char *val = value;
  int found;
  char **string_field;
  FIX_UNUSED (field);

  stop_char = 0;
  found = 0;

  if (!multiple_option && prev_given && (*prev_given || (check_ambiguity && *field_given)))
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: `--%s' (`-%c') option given more than once%s\n", 
               package_name, long_opt, short_opt,
               (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: `--%s' option given more than once%s\n", 
               package_name, long_opt,
               (additional_error ? additional_error : ""));
      return 1; /* failure */
    }

  FIX_UNUSED (default_value);
    
  if (field_given && *field_given && ! override)
    return 0;
  if (prev_given)
    (*prev_given)++;
  if (field_given)
    (*field_given)++;
  if (possible_values)
    val = possible_values[found];

  switch(arg_type) {
  case ARG_FLAG:
    *((int *)field) = !*((int *)field);
    break;
  case ARG_INT:
    if (val) *((int *)field) = strtol (val, &stop_char, 0);
    break;
  case ARG_LONG:
    if (val) *((long *)field) = (long)strtol (val, &stop_char, 0);
    break;
  case ARG_STRING:
    if (val) {
      string_field = (char **)field;
      if (!no_free && *string_field)
        free (*string_field); /* free previous string */
      *string_field = gengetopt_strdup (val);
    }
    break;
  default:
    break;
  };

  /* check numeric conversion */
  switch(arg_type) {
  case ARG_INT:
  case ARG_LONG:
    if (val && !(stop_char && *stop_char == '\0')) {
      fprintf(stderr, "%s: invalid numeric value: %s\n", package_name, val);
      return 1; /* failure */
    }
    break;
  default:
    ;
  };

  /* store the original value */
  switch(arg_type) {
  case ARG_NO:
  case ARG_FLAG:
    break;
  default:
    if (value && orig_field) {
      if (no_free) {
        *orig_field = value;
      } else {
        if (*orig_field)
          free (*orig_field); /* free previous string */
        *orig_field = gengetopt_strdup (value);
      }
    }
  };

  return 0; /* OK */
}


int
mer_counter_cmdline_internal (
  int argc, char **argv, struct mer_counter_args *args_info,
                        struct mer_counter_cmdline_params *params, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error = 0;
  struct mer_counter_args local_args_info;
  
  int override;
  int initialize;
  int check_required;
  int check_ambiguity;
  
  package_name = argv[0];
  
  override = params->override;
  initialize = params->initialize;
  check_required = params->check_required;
  check_ambiguity = params->check_ambiguity;

  if (initialize)
    mer_counter_cmdline_init (args_info);

  mer_counter_cmdline_init (&local_args_info);

  optarg = 0;
  optind = 0;
  opterr = params->print_errors;
  optopt = '?';

  while (1)
    {
      int option_index = 0;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "full-help",	0, NULL, 0 },
        { "version",	0, NULL, 'V' },
        { "mer-len",	1, NULL, 'm' },
        { "size",	1, NULL, 's' },
        { "threads",	1, NULL, 't' },
        { "output",	1, NULL, 'o' },
        { "counter-len",	1, NULL, 'c' },
        { "out-counter-len",	1, NULL, 0 },
        { "both-strands",	0, NULL, 'C' },
        { "reprobes",	1, NULL, 'p' },
        { "raw",	0, NULL, 'r' },
        { "quake",	0, NULL, 'q' },
        { "quality-start",	1, NULL, 0 },
        { "min-quality",	1, NULL, 0 },
        { "lower-count",	1, NULL, 'L' },
        { "upper-count",	1, NULL, 'U' },
        { "matrix",	1, NULL, 0 },
        { "timing",	1, NULL, 0 },
        { "no-write",	0, NULL, 'w' },
        { "measure",	0, NULL, 'u' },
        { "buffers",	1, NULL, 0 },
        { "buffer-size",	1, NULL, 0 },
        { "out-buffer-size",	1, NULL, 0 },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVm:s:t:o:c:Cp:rqL:U:wu", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          mer_counter_cmdline_print_help ();
          mer_counter_cmdline_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
          mer_counter_cmdline_print_version ();
          mer_counter_cmdline_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'm':	/* Length of mer.  */
        
        
          if (update_arg( (void *)&(args_info->mer_len_arg), 
               &(args_info->mer_len_orig), &(args_info->mer_len_given),
              &(local_args_info.mer_len_given), optarg, 0, 0, ARG_INT,
              check_ambiguity, override, 0, 0,
              "mer-len", 'm',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* Hash size.  */
        
        
          if (update_arg( (void *)&(args_info->size_arg), 
               &(args_info->size_orig), &(args_info->size_given),
              &(local_args_info.size_given), optarg, 0, 0, ARG_LONG,
              check_ambiguity, override, 0, 0,
              "size", 's',
              additional_error))
            goto failure;
        
          break;
        case 't':	/* Number of threads.  */
        
        
          if (update_arg( (void *)&(args_info->threads_arg), 
               &(args_info->threads_orig), &(args_info->threads_given),
              &(local_args_info.threads_given), optarg, 0, "1", ARG_INT,
              check_ambiguity, override, 0, 0,
              "threads", 't',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* Output prefix.  */
        
        
          if (update_arg( (void *)&(args_info->output_arg), 
               &(args_info->output_orig), &(args_info->output_given),
              &(local_args_info.output_given), optarg, 0, "mer_counts", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output", 'o',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Length of counting field.  */
        
        
          if (update_arg( (void *)&(args_info->counter_len_arg), 
               &(args_info->counter_len_orig), &(args_info->counter_len_given),
              &(local_args_info.counter_len_given), optarg, 0, "7", ARG_INT,
              check_ambiguity, override, 0, 0,
              "counter-len", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'C':	/* Count both strand, canonical representation.  */
        
        
          if (update_arg((void *)&(args_info->both_strands_flag), 0, &(args_info->both_strands_given),
              &(local_args_info.both_strands_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "both-strands", 'C',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Maximum number of reprobes.  */
        
        
          if (update_arg( (void *)&(args_info->reprobes_arg), 
               &(args_info->reprobes_orig), &(args_info->reprobes_given),
              &(local_args_info.reprobes_given), optarg, 0, "62", ARG_INT,
              check_ambiguity, override, 0, 0,
              "reprobes", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'r':	/* Write raw database.  */
        
        
          if (update_arg((void *)&(args_info->raw_flag), 0, &(args_info->raw_given),
              &(local_args_info.raw_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "raw", 'r',
              additional_error))
            goto failure;
        
          break;
        case 'q':	/* Quake compatibility mode.  */
        
        
          if (update_arg((void *)&(args_info->quake_flag), 0, &(args_info->quake_given),
              &(local_args_info.quake_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "quake", 'q',
              additional_error))
            goto failure;
        
          break;
        case 'L':	/* Don't output k-mer with count < lower-count.  */
        
        
          if (update_arg( (void *)&(args_info->lower_count_arg), 
               &(args_info->lower_count_orig), &(args_info->lower_count_given),
              &(local_args_info.lower_count_given), optarg, 0, 0, ARG_LONG,
              check_ambiguity, override, 0, 0,
              "lower-count", 'L',
              additional_error))
            goto failure;
        
          break;
        case 'U':	/* Don't output k-mer with count > upper-count.  */
        
        
          if (update_arg( (void *)&(args_info->upper_count_arg), 
               &(args_info->upper_count_orig), &(args_info->upper_count_given),
              &(local_args_info.upper_count_given), optarg, 0, 0, ARG_LONG,
              check_ambiguity, override, 0, 0,
              "upper-count", 'U',
              additional_error))
            goto failure;
        
          break;
        case 'w':	/* Don't write database.  */
        
        
          if (update_arg((void *)&(args_info->no_write_flag), 0, &(args_info->no_write_given),
              &(local_args_info.no_write_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "no-write", 'w',
              additional_error))
            goto failure;
        
          break;
        case 'u':	/* Write usage statistics.  */
        
        
          if (update_arg((void *)&(args_info->measure_flag), 0, &(args_info->measure_given),
              &(local_args_info.measure_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "measure", 'u',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
          if (strcmp (long_options[option_index].name, "full-help") == 0) {
            mer_counter_cmdline_print_full_help ();
            mer_counter_cmdline_free (&local_args_info);
            exit (EXIT_SUCCESS);
          }

          /* Length of counter field in output.  */
          if (strcmp (long_options[option_index].name, "out-counter-len") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->out_counter_len_arg), 
                 &(args_info->out_counter_len_orig), &(args_info->out_counter_len_given),
                &(local_args_info.out_counter_len_given), optarg, 0, "4", ARG_INT,
                check_ambiguity, override, 0, 0,
                "out-counter-len", '-',
                additional_error))
              goto failure;
          
          }
          /* Starting ASCII for quality values.  */
          else if (strcmp (long_options[option_index].name, "quality-start") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->quality_start_arg), 
                 &(args_info->quality_start_orig), &(args_info->quality_start_given),
                &(local_args_info.quality_start_given), optarg, 0, "64", ARG_INT,
                check_ambiguity, override, 0, 0,
                "quality-start", '-',
                additional_error))
              goto failure;
          
          }
          /* Minimum quality. A base with lesser quality becomes an N.  */
          else if (strcmp (long_options[option_index].name, "min-quality") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->min_quality_arg), 
                 &(args_info->min_quality_orig), &(args_info->min_quality_given),
                &(local_args_info.min_quality_given), optarg, 0, "0", ARG_INT,
                check_ambiguity, override, 0, 0,
                "min-quality", '-',
                additional_error))
              goto failure;
          
          }
          /* Hash function binary matrix.  */
          else if (strcmp (long_options[option_index].name, "matrix") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->matrix_arg), 
                 &(args_info->matrix_orig), &(args_info->matrix_given),
                &(local_args_info.matrix_given), optarg, 0, 0, ARG_STRING,
                check_ambiguity, override, 0, 0,
                "matrix", '-',
                additional_error))
              goto failure;
          
          }
          /* Print timing information.  */
          else if (strcmp (long_options[option_index].name, "timing") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->timing_arg), 
                 &(args_info->timing_orig), &(args_info->timing_given),
                &(local_args_info.timing_given), optarg, 0, 0, ARG_STRING,
                check_ambiguity, override, 0, 0,
                "timing", '-',
                additional_error))
              goto failure;
          
          }
          /* Number of buffers per thread.  */
          else if (strcmp (long_options[option_index].name, "buffers") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->buffers_arg), 
                 &(args_info->buffers_orig), &(args_info->buffers_given),
                &(local_args_info.buffers_given), optarg, 0, 0, ARG_LONG,
                check_ambiguity, override, 0, 0,
                "buffers", '-',
                additional_error))
              goto failure;
          
          }
          /* Size of buffers.  */
          else if (strcmp (long_options[option_index].name, "buffer-size") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->buffer_size_arg), 
                 &(args_info->buffer_size_orig), &(args_info->buffer_size_given),
                &(local_args_info.buffer_size_given), optarg, 0, "8192", ARG_LONG,
                check_ambiguity, override, 0, 0,
                "buffer-size", '-',
                additional_error))
              goto failure;
          
          }
          /* Size of output buffer per thread.  */
          else if (strcmp (long_options[option_index].name, "out-buffer-size") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->out_buffer_size_arg), 
                 &(args_info->out_buffer_size_orig), &(args_info->out_buffer_size_given),
                &(local_args_info.out_buffer_size_given), optarg, 0, "20000000", ARG_LONG,
                check_ambiguity, override, 0, 0,
                "out-buffer-size", '-',
                additional_error))
              goto failure;
          
          }
          
          break;
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", MER_COUNTER_CMDLINE_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */



  if (check_required)
    {
      error += mer_counter_cmdline_required2 (args_info, argv[0], additional_error);
    }

  mer_counter_cmdline_release (&local_args_info);

  if ( error )
    return (EXIT_FAILURE);

  if (optind < argc)
    {
      int i = 0 ;
      int found_prog_name = 0;
      /* whether program name, i.e., argv[0], is in the remaining args
         (this may happen with some implementations of getopt,
          but surely not with the one included by gengetopt) */

      i = optind;
      while (i < argc)
        if (argv[i++] == argv[0]) {
          found_prog_name = 1;
          break;
        }
      i = 0;

      args_info->inputs_num = argc - optind - found_prog_name;
      args_info->inputs =
        (char **)(malloc ((args_info->inputs_num)*sizeof(char *))) ;
      while (optind < argc)
        if (argv[optind++] != argv[0])
          args_info->inputs[ i++ ] = gengetopt_strdup (argv[optind-1]) ;
    }

  return 0;

failure:
  
  mer_counter_cmdline_release (&local_args_info);
  return (EXIT_FAILURE);
}
