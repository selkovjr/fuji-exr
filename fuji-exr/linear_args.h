#include <argp.h>

// ------------------------
// ## linear command parser
//
struct arg_linear {
  bool hr;
  char* input_file_0;
  char* input_file_1;
  char* output_file;
};

static char args_doc_linear[] = "[bayer.tiff | -h bayer_0.tiff bayer_1.tiff] output.tiff";

static char doc_linear[] =
"\n"
"Merge and interpolate the two subframes of the EXR Bayer array,\n"
"or interpolate a single frame."
"\n"
"Input:\n"
"  One or two raw Bayer frames extracted with dcraw from\n"
"  an EXR image:\n"
"\n"
"    dcraw -v -w -d -s all -4 -T <source.RAF>\n"
"\n"
"  The number of input frames depends on the EXR mode and rendering intent.\n"
"  The HR (high-resolution) images can be assembled from frames shot in any\n"
"  mode, and in that case, the frames must be interpolated together."
"\n"
"Output:\n"
"  Interpolated TIFF image"
"\n"
"\v"
"The algorithm proceeds as follows:\n"
"\n"
"  1. The two input frames are rotated 45° CCW and merged\n"
"     (interleaved) to reconstruct the high-resoluttion\n"
"     EXR matrix. This step is not done in the case of a\n"
"     single input frame."
"\n"
"  2. A simple linear interpolation with symmetrical unbiased\n"
"     stencils is used to fill the missing pixels in each color\n"
"     plane.\n"
"\n"
"  When assembling an HR image, if any transformations need to be\n"
"  applied to the input frames (such as exposure adjustment to\n"
"  frames shot in the SN mode), such transformations must precede\n"
"  interpolation.\n"
"\n"
;

static error_t parse_linear_command(int key, char* arg, struct argp_state* state) {
  struct arg_linear* arguments = (struct arg_linear*)state->input;
  char **nonopt;

  assert( arguments );

  switch(key) {
    case 'h':
      arguments->hr = true;
      break;

    case ARGP_KEY_NO_ARGS:
      argp_usage (state);

    case ARGP_KEY_ARG: // non-option argument
      // Here we know that state->arg_num == 0, since we force option parsing
      // to end before any non-option arguments can be seen
      nonopt = &state->argv[state->next];
      state->next = state->argc; // we're done

      if (arguments->hr) {
        arguments->input_file_0 = arg;
        arguments->input_file_1 = nonopt[0];
        arguments->output_file = nonopt[1];
      }
      else {
        arguments->input_file_0 = arg;
        arguments->output_file = nonopt[0];
      }
      break;

    case ARGP_KEY_END:
      if (arguments->hr) {
        if (state->arg_num < 3) {
          argp_error(state, "Not enough arguments");
        }
        if (state->arg_num > 3) {
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
static struct argp_option options_linear[] = {
  {"high-res", 'h', 0, 0, "merge input frames into a tilted HR Bayer array" },
  { 0 }
};

static struct argp argp_linear = {
  options_linear,
  parse_linear_command,
  args_doc_linear,
  doc_linear
};
#pragma GCC diagnostic pop

#define PARSE_ARGS_LINEAR \
  struct arg_linear args; \
  int    argc = state->argc - state->next + 1; \
  char** argv = &state->argv[state->next - 1]; \
  char*  argv0 =  argv[0]; \
  argv[0] = (char *)malloc(strlen((char *)(state->name)) + strlen(" linear") + 1); \
  if (!argv[0]) argp_failure(state, 1, ENOMEM, 0); \
  sprintf(argv[0], "%s linear", state->name); \
  args.hr = false; \
  argp_parse(&argp_linear, argc, argv, ARGP_IN_ORDER, &argc, &args); \
  free(argv[0]); \
  argv[0] = argv0; \
  state->next += argc - 1; \

