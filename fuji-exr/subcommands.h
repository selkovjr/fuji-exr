#include <argp.h>

void run_linear (struct argp_state* state);
void run_rotate (struct argp_state* state);
void run_sdd (struct argp_state* state);

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

