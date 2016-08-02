#include <argp.h>

// ----------------------------------------------------------------------------------
// ## SSDD command parser
//
struct arg_ssdd {
  bool interlaced_cfa;
  char* geometry;
  char* input_file_0;
  char* input_file_1;
  char* input_file_2;
  char* output_file;
};

static char args_doc_ssdd[] = "[-x WxH r.tiff g.tiff b.tiff | bayer_0.tiff bayer_1.tiff] output.tiff";

static char doc_ssdd[] =
"\n"
"Self-similarity-driven debayering\n"
"\n"
"Input:\n"
"  Two raw Bayer frames extracted with dcraw from\n"
"  an HR (high-resolution) EXR image:\n"
"\n"
"    dcraw -v -w -d -s all -4 -T <source.RAF>\n"
"\n"
"  Or, if the -x option is given, the next three arguments\n"
"  must be the file names of the three color planes (R, G, B)\n"
"  of an interlaced high-resolution EXR array.\n"
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

static error_t parse_ssdd_command(int key, char* arg, struct argp_state* state) {
  struct arg_ssdd* arguments = (struct arg_ssdd*)state->input;
  char **nonopt;

  assert( arguments );

  switch(key) {
    case 'x':
      arguments->interlaced_cfa = true;
      arguments->geometry = arg;
      break;

    case ARGP_KEY_NO_ARGS:
      argp_usage (state);

    case ARGP_KEY_ARG: // non-option argument
      // Here we know that state->arg_num == 0, since we force option parsing
      // to end before any non-option arguments can be seen
      nonopt = &state->argv[state->next];
      state->next = state->argc; // we're done

      arguments->input_file_0 = arg;
      if (arguments->interlaced_cfa) {
        arguments->input_file_1 = nonopt[0];
        arguments->input_file_2 = nonopt[1];
        arguments->output_file = nonopt[2];
      }
      else {
        arguments->input_file_1 = nonopt[0];
        arguments->output_file = nonopt[1];
      }
      break;

    case ARGP_KEY_END:
      if (arguments->interlaced_cfa) {
        if (state->arg_num < 4) {
          argp_error(state, "Not enough arguments");
        }
        if (state->arg_num > 4) {
          argp_error(state, "Extra arguments");
        }
      }
      else {
        if (state->arg_num < 2) {
          argp_error(state, "Not enough arguments");
        }
        if (state->arg_num > 2) {
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
static struct argp_option options_ssdd[] = {
  {"highres-exr", 'x', "WxH", 0, "Input is an interlaced high-resolution EXR array with the CFA geometry of WxH" },
  { 0 }
};

static struct argp argp_ssdd = {
  options_ssdd,
  parse_ssdd_command,
  args_doc_ssdd,
  doc_ssdd
};
#pragma GCC diagnostic pop

#define PARSE_ARGS_SSDD \
  struct arg_ssdd args; \
  int    argc = state->argc - state->next + 1; \
  char** argv = &state->argv[state->next - 1]; \
  char*  argv0 =  argv[0]; \
  argv[0] = (char *)malloc(strlen((char *)(state->name)) + strlen(" ssdd") + 1); \
  if (!argv[0]) argp_failure(state, 1, ENOMEM, 0); \
  sprintf(argv[0], "%s ssdd", state->name); \
  args.interlaced_cfa = false; \
  argp_parse(&argp_ssdd, argc, argv, ARGP_IN_ORDER, &argc, &args); \
  free(argv[0]); \
  argv[0] = argv0; \
  state->next += argc - 1;

