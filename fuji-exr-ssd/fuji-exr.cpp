#include <argp.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "subcommands.h"

const char *argp_program_version = VERSION " " __DATE__;

//-------------------------------------------------------------------
// Command-line option parsing.
// Example at https://gist.github.com/sam-github/57f49711cd9073b35d22
//

// ## Top-level parser
static char doc_toplevel[1000];

static error_t parse_toplevel (int key, char *arg, struct argp_state *state) {
  switch (key) {
    case ARGP_KEY_ARG:
      assert( arg );
      if(strcmp(arg, "ssd") == 0) {
        run_ssd(state);
      }
      else if(strcmp(arg, "linear") == 0) {
        run_linear(state);
      }
      else if(strcmp(arg, "rotate") == 0) {
        run_rotate(state);
      }
      else {
        argp_error(state, "%s is not a valid command", arg);
      }
      break;

    case ARGP_KEY_END:
      if (state->arg_num < 2)
        argp_usage (state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
static struct argp argp_toplevel = {
  0, // no options
  parse_toplevel,
  "<COMMAND> <ARGS>",
  doc_toplevel
};
#pragma GCC diagnostic pop


int main(int argc, char **argv) {
  sprintf(
    doc_toplevel,
    "\n"
    "Utilities for processing Fuji EXR sensor data\n"
    "Version: %s\n"
    "\n"
    "Command: linear  interpolate channels without debayering\n"
    "         ssd     self-similarity-driven debayering\n"
    "         db      Duran-Buades debayering\n"
    "\n",
    argp_program_version
  );

  argp_parse (&argp_toplevel, argc, argv, ARGP_IN_ORDER, NULL, NULL);

  return(EXIT_SUCCESS);
}

