#include <argp.h>

// ----------------------------------------------------------------------------------
// ## SDD command parser
//
struct arg_sdd {
  bool merged_cfa;
  char* geometry;
  char* input_file_0;
  char* input_file_1;
  char* input_file_2;
  char* output_file;
};

static char args_doc_sdd[] = "[-m WxH r.tiff g.tiff b.tiff | bayer_0.tiff bayer_1.tiff] output.tiff";

static char doc_sdd[] =
"\n"
"Self-similarity-driven debayering\n"
"\n"
"Input:\n"
"  Two raw Bayer frames extracted with dcraw from\n"
"  an HR (high-resolution) EXR image:\n"
"\n"
"    dcraw -v -w -d -s all -4 -T <source.RAF>\n"
"\n"
"  Or, if the -m option is given, image geometry followed\n"
"  by the three color planes of a merged HR Bayer array.\n"
"\n"
"  Use the -m option to operate on preprocessed inputs.\n"
"\n"
"Output:\n"
"  Interpolated and filtered TIFF image"
"\n"
"\v"
"The algorithm proceeds as follows:\n"
"\n"
"  1. The two input frames are rotated 45° CCW and merged\n"
"     (interleaved) to reconstruct the high-resoluttion\n"
"     EXR matrix.\n"
"\n"
"  2. An algorithm analogous to Adams-Hamilton but with\n"
"     EXR-specific stencils is used to do directional\n"
"     interpolation of the green channel. Then bilinear\n"
"     interpolation is applied to B-G and R-G differences.\n"
"\n"
"  3. A non-local means filter is applied to each channel,\n"
"     using the weighted average of the channel's raw values.\n"
"\n"
"  4. Chromatic noise is suppressed by a median filter.\n"
"\n"
"  5. The interpolated image is rotated to restore its\n"
"     photographic orientation.\n"
"\n"
"Author: Gene Selkov\n"
"\n"
"Idea and portions of code from:\n"
"\n"
"  Antoni Buades, Bartomeu Coll,\n"
"  Jean-Michel Morel, and Catalina Sbert,\n"
"  Self-similarity Driven Demosaicking,\n"
"  Image Processing On Line, 1 (2011).\n"
"  http://dx.doi.org/10.5201/ipol.2011.bcms-ssdd\n"
"\n"
;

static error_t parse_sdd_command(int key, char* arg, struct argp_state* state) {
  struct arg_sdd* arguments = (struct arg_sdd*)state->input;
  char **nonopt;

  assert( arguments );

  switch(key) {
    case 'm':
      arguments->merged_cfa = true;
      break;

    case ARGP_KEY_NO_ARGS:
      argp_usage (state);

    case ARGP_KEY_ARG: // non-option argument
      // Here we know that state->arg_num == 0, since we force option parsing
      // to end before any non-option arguments can be seen
      nonopt = &state->argv[state->next];
      state->next = state->argc; // we're done

      if (arguments->merged_cfa) {
        arguments->geometry = arg;
        arguments->input_file_0 = nonopt[0];
        arguments->input_file_1 = nonopt[1];
        arguments->input_file_2 = nonopt[2];
        arguments->output_file = nonopt[3];
      }
      else {
        arguments->input_file_0 = arg;
        arguments->input_file_1 = nonopt[0];
        arguments->output_file = nonopt[1];
      }
      break;

    case ARGP_KEY_END:
      if (arguments->merged_cfa) {
        if (state->arg_num < 5) {
          argp_error(state, "Not enough arguments");
        }
        if (state->arg_num > 5) {
          argp_error(state, "Extra arguments");
        }
      }
      else {
        if (state->arg_num < 3) {
          argp_error(state, "Not enough arguments");
        }
        if (state->arg_num > 3) {
          argp_error(state, "Extra arguments");
        }
      }
      break;

    default:
      return ARGP_ERR_UNKNOWN;
  }

  return 0;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
static struct argp_option options_sdd[] = {
  {"merged", 'm', 0, 0, "Input is a merged HR Bayer array" },
  { 0 }
};

static struct argp argp_sdd = {
  options_sdd,
  parse_sdd_command,
  args_doc_sdd,
  doc_sdd
};
#pragma GCC diagnostic pop

#define PARSE_ARGS_SDD \
  struct arg_sdd args; \
  int    argc = state->argc - state->next + 1; \
  char** argv = &state->argv[state->next - 1]; \
  char*  argv0 =  argv[0]; \
  argv[0] = (char *)malloc(strlen((char *)(state->name)) + strlen(" sdd") + 1); \
  if (!argv[0]) argp_failure(state, 1, ENOMEM, 0); \
  sprintf(argv[0], "%s sdd", state->name); \
  args.merged_cfa = false; \
  argp_parse(&argp_sdd, argc, argv, ARGP_IN_ORDER, &argc, &args); \
  free(argv[0]); \
  argv[0] = argv0; \
  state->next += argc - 1;

