/** @file cite_cmdline.hpp
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.4
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef CITE_CMDLINE_H
#define CITE_CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CITE_CMDLINE_PACKAGE
/** @brief the program name (used for printing errors) */
#define CITE_CMDLINE_PACKAGE "jellyfish cite"
#endif

#ifndef CITE_CMDLINE_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CITE_CMDLINE_PACKAGE_NAME "jellyfish cite"
#endif

#ifndef CITE_CMDLINE_VERSION
/** @brief the program version */
#define CITE_CMDLINE_VERSION VERSION
#endif

/** @brief Where the command line options are stored */
struct cite_args
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  int bibtex_flag;	/**< @brief Bibtex format (default=off).  */
  const char *bibtex_help; /**< @brief Bibtex format help description.  */
  char * output_arg;	/**< @brief Output file (default='/dev/fd/1').  */
  char * output_orig;	/**< @brief Output file original value given at command line.  */
  const char *output_help; /**< @brief Output file help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int bibtex_given ;	/**< @brief Whether bibtex was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */

} ;

/** @brief The additional parameters to pass to parser functions */
struct cite_cmdline_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure cite_args (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure cite_args (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *cite_args_purpose;
/** @brief the usage string of the program */
extern const char *cite_args_usage;
/** @brief all the lines making the help output */
extern const char *cite_args_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cite_cmdline (int argc, char **argv,
  struct cite_args *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cite_cmdline_ext() instead
 */
int cite_cmdline2 (int argc, char **argv,
  struct cite_args *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cite_cmdline_ext (int argc, char **argv,
  struct cite_args *args_info,
  struct cite_cmdline_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cite_cmdline_dump(FILE *outfile,
  struct cite_args *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cite_cmdline_file_save(const char *filename,
  struct cite_args *args_info);

/**
 * Print the help
 */
void cite_cmdline_print_help(void);
/**
 * Print the version
 */
void cite_cmdline_print_version(void);

/**
 * Initializes all the fields a cite_cmdline_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cite_cmdline_params_init(struct cite_cmdline_params *params);

/**
 * Allocates dynamically a cite_cmdline_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cite_cmdline_params structure
 */
struct cite_cmdline_params *cite_cmdline_params_create(void);

/**
 * Initializes the passed cite_args structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cite_cmdline_init (struct cite_args *args_info);
/**
 * Deallocates the string fields of the cite_args structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cite_cmdline_free (struct cite_args *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cite_cmdline_required (struct cite_args *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CITE_CMDLINE_H */
