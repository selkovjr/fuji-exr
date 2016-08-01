#include <argp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <error.h>
#include <omp.h>
#include <ctime>

#include "cfa_mask.h"
#include "io_tiff.h"
#include "tiffio.h"
#include "libAuxiliary.h"

#define DIAG 1.4142136
#define DIAG12 2.236 // sqrt(5)

#define DUMP_STAGES
#undef DEBUG_GREEN

// --------------------
// ## rotate command parser

static char doc_rotate[] =
"\n"
"Merge and interpolate the two halves of the EXR-HR Bayer array \n"
"\n"
"Input:\n"
"  Two raw Bayer frames extracted with dcraw from\n"
"  an HR (high-resolution) EXR image:\n"
"\n"
"    dcraw -v -w -d -s all -4 -T <source.RAF>\n"
"\n"
"Output:\n"
"  Interpolated TIFF image"
"\n"
"\v"
"The algorithm proceeds as follows:\n"
"\n"
"  1. The two input frames are rotated 45° CCW and merged\n"
"     (interleaved) to reconstruct the high-resoluttion\n"
"     EXR matrix.\n"
"\n"
"  2. A simple linear interpolation with symmetrical stencils\n"
"     is used to fill the missing pixels in each color plane.\n"
"\n"
;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
static error_t parse_rotate_command(int key, char* arg, struct argp_state* state) {
  switch(key) {
    case ARGP_KEY_END:
      if (state->arg_num > 0) {
        argp_error(state, "Extra arguments");
      }
      break;

    default:
      return ARGP_ERR_UNKNOWN;
  }

  return 0;
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
static struct argp argp_rotate = {
  0, // no options
  parse_rotate_command,
  0, // no arguments
  doc_rotate
};
#pragma GCC diagnostic pop

// --------------------

void run_rotate (struct argp_state* state) {
  // command-line stuff
  int    argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  char*  argv0 =  argv[0];
  argv[0] = (char *)malloc(strlen((char *)(state->name)) + strlen(" rotate") + 1);
  if (!argv[0]) argp_failure(state, 1, ENOMEM, 0);
  sprintf(argv[0], "%s rotate", state->name);
  argp_parse(&argp_rotate, argc, argv, ARGP_IN_ORDER, &argc, NULL);
  free(argv[0]);
  argv[0] = argv0;
  state->next += argc - 1;

  // Processing starts here
  long cfaWidth = 20;
  long cfaHeight = 15;
  long width = cfaWidth + cfaHeight;
  float *data_in, *data_out;

  if (NULL == (data_in = (float *) malloc(sizeof(float) * width * width * 3))) {
    fprintf(stderr, "allocation error. not enough memory?\n");
    exit(EXIT_FAILURE);
  }

  for (long i = 0; i < width * width * 3; i++) {
    data_in[i] = 0;
  }

  unsigned char *mask = exr_cfa_mask(width, width, cfaWidth, cfaHeight);

  for (long y = 0; y < width; y++) {
    fprintf(stderr, "%03ld: ", y);
    for (long x = 0; x < width; x++) {
      long p = y * width + x;
      switch (mask[p]) {
       case REDPOSITION:
         fprintf(stderr, "r ");
         data_in[p] = 65535;
         break;
       case GREENPOSITION:
         fprintf(stderr, "g ");
         data_in[p + width * width] = 65535;
         break;
       case BLUEPOSITION:
         fprintf(stderr, "b ");
         data_in[p + 2 * width * width] = 65535;
         break;
       default:
         fprintf(stderr, ". ");
      }
    }
    fprintf(stderr, "\n");
  }

  write_tiff_rgb_f32("mask.tiff", data_in, width, width);

  delete[] mask;

  mask = bggr_cfa_mask(cfaWidth, cfaHeight);

  for (long y = 0; y < cfaHeight; y++) {
    fprintf(stderr, "%03ld: ", y);
    for (long x = 0; x < cfaWidth; x++) {
      long p = y * cfaWidth + x;
      switch (mask[p]) {
       case REDPOSITION:
         fprintf(stderr, "r ");
         data_in[p] = 65535;
         break;
       case GREENPOSITION:
         fprintf(stderr, "g ");
         data_in[p + cfaWidth * cfaHeight] = 65535;
         break;
       case BLUEPOSITION:
         fprintf(stderr, "b ");
         data_in[p + 2 * cfaWidth * cfaHeight] = 65535;
         break;
       default:
         fprintf(stderr, ". ");
      }
    }
    fprintf(stderr, "\n");
  }

  write_tiff_rgb_f32("mask.tiff", data_in, cfaWidth, cfaHeight);
  exit(0);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
{

  int height = width;

  int i, row, col;
  double step = sqrt(0.5); // Horizontal or vertical CFA step projected onto
                           // source-plane axes
  float r, c;              // Y- and X-coords in the source plane
  unsigned ur, uc;         // Y- and X-coords of the nearest source pixel
  float fr, fc;            // Y- and X-distance from (r, c) to nearest pixel
  ushort rotWidth, rotHeight;

  fprintf (stderr, "Rotating image 45 degrees...\n");

  // Inflated (√2) target image co-ordinates
  rotWidth = cfaWidth / step;
  rotHeight = (height - cfaWidth) / step;

  // img = (ushort (*)[4]) calloc (rotHeight, rotWidth * sizeof *img);
  if (NULL == (data_out = (float *) malloc(sizeof(float) * rotWidth * rotHeight * 3))) {
    fprintf(stderr, "allocation error. not enough memory?\n");
    exit(EXIT_FAILURE);
  }

  // Row and col are co-ordinates in the inflated target image.
  for (row = 0; row < rotHeight; row++) {
    printf("%d: ", row);
    for (col = 0; col < rotWidth; col++) {
      // Reverse mapping: find co-ordinates (r, c) in the rotated
      // CFA plane whose ushort casts (ur, uc) point to the source
      // CFA pixel.
      ur = r = cfaWidth + (row - col) * step;
      uc = c = (row + col) * step;

      // leave margins in the source image for the stencil
      if (ur > (unsigned)(height - 2) || uc > (unsigned)(width - 2)) continue;

      fr = r - ur;
      fc = c - uc;

      for (i = 0; i < 3; i++) { // for each color plane

        // David Coffin's original stencil (on an array of pixels)
        //
        //   pix = img + ur * iwidth + uc;
        //   img[row * wide + col][i] =
        //     (/* + */ pix[    0][i]*(1 - fc) + /* E  */ pix[        1][i] * fc) * (1 - fr) +
        //     (/* S */ pix[width][i]*(1 - fc) + /* SE */ pix[width + 1][i] * fc) * fr;
        //
        // Same stencil reformulated for stacked color planes
        //
        data_out[row * rotWidth + col + i * rotWidth * rotHeight] =
          (1 - fr) * (
            (1 - fc) * data_in[ur * width + uc + i * width * height]              // +
            +
                  fc * data_in[ur * width + uc + i * width * height + 1]          // E
          )
          +
          fr * (
            (1 - fc) * data_in[ur * width + uc + i * width * height + width]      // S
            +
                  fc * data_in[ur * width + uc + i * width * height + width + 1]  // SE
          )
          ;
      } // each color plane
    }
    printf("\n");
  }
  write_tiff_rgb_f32("out.tiff", data_out, rotWidth, rotHeight);
}
#pragma GCC diagnostic pop

  delete[] mask;

  exit(EXIT_SUCCESS);

} // run_rotate()


