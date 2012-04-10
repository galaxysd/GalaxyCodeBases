/*
  File autogenerated by gengetopt version 2.22.4
  generated with the following command:
  gengetopt --show-required --default-option -c cc -H hpp -F hash_fastq_merge_cmdline -f hash_fastq_merge_cmdline -a hash_fastq_merge_args --unamed-opts=database.jf

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

#include "hash_fastq_merge_cmdline.hpp"

const char *hash_fastq_merge_args_purpose = "Merge quake databases";

const char *hash_fastq_merge_args_usage = "Usage: jellyfish merge [OPTIONS]... [database.jf]...";

const char *hash_fastq_merge_args_description = "";

const char *hash_fastq_merge_args_full_help[] = {
  "  -h, --help                  Print help and exit",
  "      --full-help             Print help, including hidden options, and exit",
  "  -V, --version               Print version and exit",
  "  -s, --size=LONG             Merged hash table size (mandatory)",
  "  -m, --mer-len=INT           Mer length (mandatory)",
  "  -o, --output=STRING         Output file  (default=`merged.jf')",
  "  -p, --reprobes=INT          Maximum number of reprobes  (default=`62')",
  "      --out-buffer-size=LONG  Size of output buffer per thread  \n                                (default=`20000000')",
    0
};

static void
init_help_array(void)
{
  hash_fastq_merge_args_help[0] = hash_fastq_merge_args_full_help[0];
  hash_fastq_merge_args_help[1] = hash_fastq_merge_args_full_help[1];
  hash_fastq_merge_args_help[2] = hash_fastq_merge_args_full_help[2];
  hash_fastq_merge_args_help[3] = hash_fastq_merge_args_full_help[3];
  hash_fastq_merge_args_help[4] = hash_fastq_merge_args_full_help[4];
  hash_fastq_merge_args_help[5] = hash_fastq_merge_args_full_help[5];
  hash_fastq_merge_args_help[6] = hash_fastq_merge_args_full_help[6];
  hash_fastq_merge_args_help[7] = 0; 
  
}

const char *hash_fastq_merge_args_help[8];

typedef enum {ARG_NO
  , ARG_STRING
  , ARG_INT
  , ARG_LONG
} hash_fastq_merge_cmdline_arg_type;

static
void clear_given (struct hash_fastq_merge_args *args_info);
static
void clear_args (struct hash_fastq_merge_args *args_info);

static int
hash_fastq_merge_cmdline_internal (int argc, char **argv, struct hash_fastq_merge_args *args_info,
                        struct hash_fastq_merge_cmdline_params *params, const char *additional_error);

static int
hash_fastq_merge_cmdline_required2 (struct hash_fastq_merge_args *args_info, const char *prog_name, const char *additional_error);

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct hash_fastq_merge_args *args_info)
{
  args_info->help_given = 0 ;
  args_info->full_help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->size_given = 0 ;
  args_info->mer_len_given = 0 ;
  args_info->output_given = 0 ;
  args_info->reprobes_given = 0 ;
  args_info->out_buffer_size_given = 0 ;
}

static
void clear_args (struct hash_fastq_merge_args *args_info)
{
  FIX_UNUSED (args_info);
  args_info->size_orig = NULL;
  args_info->mer_len_orig = NULL;
  args_info->output_arg = gengetopt_strdup ("merged.jf");
  args_info->output_orig = NULL;
  args_info->reprobes_arg = 62;
  args_info->reprobes_orig = NULL;
  args_info->out_buffer_size_arg = 20000000;
  args_info->out_buffer_size_orig = NULL;
  
}

static
void init_args_info(struct hash_fastq_merge_args *args_info)
{

  init_help_array(); 
  args_info->help_help = hash_fastq_merge_args_full_help[0] ;
  args_info->full_help_help = hash_fastq_merge_args_full_help[1] ;
  args_info->version_help = hash_fastq_merge_args_full_help[2] ;
  args_info->size_help = hash_fastq_merge_args_full_help[3] ;
  args_info->mer_len_help = hash_fastq_merge_args_full_help[4] ;
  args_info->output_help = hash_fastq_merge_args_full_help[5] ;
  args_info->reprobes_help = hash_fastq_merge_args_full_help[6] ;
  args_info->out_buffer_size_help = hash_fastq_merge_args_full_help[7] ;
  
}

void
hash_fastq_merge_cmdline_print_version (void)
{
  printf ("%s %s\n",
     (strlen(HASH_FASTQ_MERGE_CMDLINE_PACKAGE_NAME) ? HASH_FASTQ_MERGE_CMDLINE_PACKAGE_NAME : HASH_FASTQ_MERGE_CMDLINE_PACKAGE),
     HASH_FASTQ_MERGE_CMDLINE_VERSION);
}

static void print_help_common(void) {
  hash_fastq_merge_cmdline_print_version ();

  if (strlen(hash_fastq_merge_args_purpose) > 0)
    printf("\n%s\n", hash_fastq_merge_args_purpose);

  if (strlen(hash_fastq_merge_args_usage) > 0)
    printf("\n%s\n", hash_fastq_merge_args_usage);

  printf("\n");

  if (strlen(hash_fastq_merge_args_description) > 0)
    printf("%s\n\n", hash_fastq_merge_args_description);
}

void
hash_fastq_merge_cmdline_print_help (void)
{
  int i = 0;
  print_help_common();
  while (hash_fastq_merge_args_help[i])
    printf("%s\n", hash_fastq_merge_args_help[i++]);
}

void
hash_fastq_merge_cmdline_print_full_help (void)
{
  int i = 0;
  print_help_common();
  while (hash_fastq_merge_args_full_help[i])
    printf("%s\n", hash_fastq_merge_args_full_help[i++]);
}

void
hash_fastq_merge_cmdline_init (struct hash_fastq_merge_args *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);

  args_info->inputs = 0;
  args_info->inputs_num = 0;
}

void
hash_fastq_merge_cmdline_params_init(struct hash_fastq_merge_cmdline_params *params)
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

struct hash_fastq_merge_cmdline_params *
hash_fastq_merge_cmdline_params_create(void)
{
  struct hash_fastq_merge_cmdline_params *params = 
    (struct hash_fastq_merge_cmdline_params *)malloc(sizeof(struct hash_fastq_merge_cmdline_params));
  hash_fastq_merge_cmdline_params_init(params);  
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
hash_fastq_merge_cmdline_release (struct hash_fastq_merge_args *args_info)
{
  unsigned int i;
  free_string_field (&(args_info->size_orig));
  free_string_field (&(args_info->mer_len_orig));
  free_string_field (&(args_info->output_arg));
  free_string_field (&(args_info->output_orig));
  free_string_field (&(args_info->reprobes_orig));
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
hash_fastq_merge_cmdline_dump(FILE *outfile, struct hash_fastq_merge_args *args_info)
{
  int i = 0;

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot dump options to stream\n", HASH_FASTQ_MERGE_CMDLINE_PACKAGE);
      return EXIT_FAILURE;
    }

  if (args_info->help_given)
    write_into_file(outfile, "help", 0, 0 );
  if (args_info->full_help_given)
    write_into_file(outfile, "full-help", 0, 0 );
  if (args_info->version_given)
    write_into_file(outfile, "version", 0, 0 );
  if (args_info->size_given)
    write_into_file(outfile, "size", args_info->size_orig, 0);
  if (args_info->mer_len_given)
    write_into_file(outfile, "mer-len", args_info->mer_len_orig, 0);
  if (args_info->output_given)
    write_into_file(outfile, "output", args_info->output_orig, 0);
  if (args_info->reprobes_given)
    write_into_file(outfile, "reprobes", args_info->reprobes_orig, 0);
  if (args_info->out_buffer_size_given)
    write_into_file(outfile, "out-buffer-size", args_info->out_buffer_size_orig, 0);
  

  i = EXIT_SUCCESS;
  return i;
}

int
hash_fastq_merge_cmdline_file_save(const char *filename, struct hash_fastq_merge_args *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", HASH_FASTQ_MERGE_CMDLINE_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  i = hash_fastq_merge_cmdline_dump(outfile, args_info);
  fclose (outfile);

  return i;
}

void
hash_fastq_merge_cmdline_free (struct hash_fastq_merge_args *args_info)
{
  hash_fastq_merge_cmdline_release (args_info);
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
hash_fastq_merge_cmdline (int argc, char **argv, struct hash_fastq_merge_args *args_info)
{
  return hash_fastq_merge_cmdline2 (argc, argv, args_info, 0, 1, 1);
}

int
hash_fastq_merge_cmdline_ext (int argc, char **argv, struct hash_fastq_merge_args *args_info,
                   struct hash_fastq_merge_cmdline_params *params)
{
  int result;
  result = hash_fastq_merge_cmdline_internal (argc, argv, args_info, params, 0);

  if (result == EXIT_FAILURE)
    {
      hash_fastq_merge_cmdline_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
hash_fastq_merge_cmdline2 (int argc, char **argv, struct hash_fastq_merge_args *args_info, int override, int initialize, int check_required)
{
  int result;
  struct hash_fastq_merge_cmdline_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = hash_fastq_merge_cmdline_internal (argc, argv, args_info, &params, 0);

  if (result == EXIT_FAILURE)
    {
      hash_fastq_merge_cmdline_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
hash_fastq_merge_cmdline_required (struct hash_fastq_merge_args *args_info, const char *prog_name)
{
  int result = EXIT_SUCCESS;

  if (hash_fastq_merge_cmdline_required2(args_info, prog_name, 0) > 0)
    result = EXIT_FAILURE;

  if (result == EXIT_FAILURE)
    {
      hash_fastq_merge_cmdline_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
hash_fastq_merge_cmdline_required2 (struct hash_fastq_merge_args *args_info, const char *prog_name, const char *additional_error)
{
  int error = 0;
  FIX_UNUSED (additional_error);

  /* checks for required options */
  if (! args_info->size_given)
    {
      fprintf (stderr, "%s: '--size' ('-s') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->mer_len_given)
    {
      fprintf (stderr, "%s: '--mer-len' ('-m') option required%s\n", prog_name, (additional_error ? additional_error : ""));
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
 * @param check_ambiguity @see hash_fastq_merge_cmdline_params.check_ambiguity
 * @param override @see hash_fastq_merge_cmdline_params.override
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
               hash_fastq_merge_cmdline_arg_type arg_type,
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
hash_fastq_merge_cmdline_internal (
  int argc, char **argv, struct hash_fastq_merge_args *args_info,
                        struct hash_fastq_merge_cmdline_params *params, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error = 0;
  struct hash_fastq_merge_args local_args_info;
  
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
    hash_fastq_merge_cmdline_init (args_info);

  hash_fastq_merge_cmdline_init (&local_args_info);

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
        { "size",	1, NULL, 's' },
        { "mer-len",	1, NULL, 'm' },
        { "output",	1, NULL, 'o' },
        { "reprobes",	1, NULL, 'p' },
        { "out-buffer-size",	1, NULL, 0 },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVs:m:o:p:", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          hash_fastq_merge_cmdline_print_help ();
          hash_fastq_merge_cmdline_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
          hash_fastq_merge_cmdline_print_version ();
          hash_fastq_merge_cmdline_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 's':	/* Merged hash table size.  */
        
        
          if (update_arg( (void *)&(args_info->size_arg), 
               &(args_info->size_orig), &(args_info->size_given),
              &(local_args_info.size_given), optarg, 0, 0, ARG_LONG,
              check_ambiguity, override, 0, 0,
              "size", 's',
              additional_error))
            goto failure;
        
          break;
        case 'm':	/* Mer length.  */
        
        
          if (update_arg( (void *)&(args_info->mer_len_arg), 
               &(args_info->mer_len_orig), &(args_info->mer_len_given),
              &(local_args_info.mer_len_given), optarg, 0, 0, ARG_INT,
              check_ambiguity, override, 0, 0,
              "mer-len", 'm',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* Output file.  */
        
        
          if (update_arg( (void *)&(args_info->output_arg), 
               &(args_info->output_orig), &(args_info->output_given),
              &(local_args_info.output_given), optarg, 0, "merged.jf", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output", 'o',
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

        case 0:	/* Long option with no short option */
          if (strcmp (long_options[option_index].name, "full-help") == 0) {
            hash_fastq_merge_cmdline_print_full_help ();
            hash_fastq_merge_cmdline_free (&local_args_info);
            exit (EXIT_SUCCESS);
          }

          /* Size of output buffer per thread.  */
          if (strcmp (long_options[option_index].name, "out-buffer-size") == 0)
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
          fprintf (stderr, "%s: option unknown: %c%s\n", HASH_FASTQ_MERGE_CMDLINE_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */



  if (check_required)
    {
      error += hash_fastq_merge_cmdline_required2 (args_info, argv[0], additional_error);
    }

  hash_fastq_merge_cmdline_release (&local_args_info);

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
  
  hash_fastq_merge_cmdline_release (&local_args_info);
  return (EXIT_FAILURE);
}
