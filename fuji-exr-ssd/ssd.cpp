/*
* Copyright (c) 2015, Gene Selkov <selkovjr@gmail.com>
*/

/*
* Portions copyright (c) 2009-2011, A. Buades <toni.buades@uib.es>
* All rights reserved.
*/

/*
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

/**
 * @mainpage Self Similarity Driven Demosaicking
 *
 * README.txt:
 * @verbinclude README.txt
 */

#include <argp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <error.h>
#include <algorithm>

#include "libdemosaic.h"
#include "io_tiff.h"
#include "tiffio.h"
#include "libAuxiliary.h"

#define BLANK 0
#define REDPOSITION 1
#define GREENPOSITION 2
#define BLUEPOSITION 3

#define DIAG 1.4142136
#define DIAG12 2.236 // sqrt(5)

#define DUMP_STAGES
#undef DEBUG_GREEN

// --------------------
// ## SSD command parser

struct arg_ssd {
  bool merged_cfa;
  char* geometry;
  char* input_file_0;
  char* input_file_1;
  char* input_file_2;
  char* output_file;
};

static char args_doc_ssd[] = "[-m WxH r.tiff g.tiff b.tiff | bayer_0.tiff bayer_1.tiff] output.tiff";

static char doc_ssd[] =
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

static error_t parse_ssd(int key, char* arg, struct argp_state* state) {
  struct arg_ssd* arguments = (struct arg_ssd*)state->input;
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
static struct argp_option options_ssd[] = {
  {"merged", 'm', 0, 0, "Input is a merged HR Bayer array" },
  { 0 }
};

static struct argp argp_ssd = {
  options_ssd,
  parse_ssd,
  args_doc_ssd,
  doc_ssd
};
#pragma GCC diagnostic pop


void run_ssd (struct argp_state* state) {
  struct arg_ssd args;
  int    argc = state->argc - state->next + 1;
  char** argv = &state->argv[state->next - 1];
  char*  argv0 =  argv[0];

  size_t nx0 = 0, ny0 = 0;
  size_t nx1 = 0, ny1 = 0;
  size_t nx2 = 0, ny2 = 0;
  char *description;
  unsigned long origWidth;
  unsigned long origHeight;
  unsigned long width, height;
  float *frame0, *frame1, *frame2; // Bayer EXR frames or TCA-corrected R, G, B
  float *data_in, *data_out;
  float *out_ptr, *end_ptr;
  bool landscape = false;


  argv[0] = (char *)malloc(strlen((char *)(state->name)) + strlen(" ssd") + 1);

  if (!argv[0]) {
    argp_failure(state, 1, ENOMEM, 0);
  }

  sprintf(argv[0], "%s ssd", state->name);

  // Default value
  args.merged_cfa = false;
  argp_parse(&argp_ssd, argc, argv, ARGP_IN_ORDER, &argc, &args);

  free(argv[0]);

  argv[0] = argv0;

  state->next += argc - 1;

  if (args.merged_cfa) {
    printf("geometry: %s\n", args.geometry);
    printf("red input file: %s\n", args.input_file_0);
    printf("green input file: %s\n", args.input_file_1);
    printf("blue input file: %s\n", args.input_file_2);

    if (sscanf(args.geometry, "%lux%lu", &origWidth, &origHeight) < 0) {
      fprintf(stderr, "error parsing image geometry '%s'\n", args.geometry);
      exit(EXIT_FAILURE);
    }
    width = height = origWidth + origHeight;

    /* TIFF 16-bit grayscale -> float input */
    if (NULL == (frame0 = read_tiff_gray16_f32(args.input_file_0, &nx0, &ny0, &description))) {
      fprintf(stderr, "error while reading from %s\n", args.input_file_0);
      exit(EXIT_FAILURE);
    }
    if (NULL == (frame1 = read_tiff_gray16_f32(args.input_file_1, &nx1, &ny1, &description))) {
      fprintf(stderr, "error while reading from %s\n", args.input_file_1);
      exit(EXIT_FAILURE);
    }
    if (NULL == (frame2 = read_tiff_gray16_f32(args.input_file_2, &nx2, &ny2, &description))) {
      fprintf(stderr, "error while reading from %s\n", args.input_file_2);
      exit(EXIT_FAILURE);
    }
    if (nx0 != nx1 or nx0 != nx2 or nx1 != nx2 or ny0 != ny1 or ny0 != ny2 or ny1 != ny2) {
      fprintf(stderr, "Input frames must have identical size. Got %ldx%ld, %ldx%ld, %ldx%ld\n", nx0, ny0, nx1, ny1, nx2, ny2);
      exit(EXIT_FAILURE);
    }
    if (nx0 != width) {
      fprintf(stderr, "Stated image geometry (%ldx%ld) does not fit input frames (%ldx%ld)\n", origWidth, origHeight, nx0, ny0);
      exit(EXIT_FAILURE);
    }
    printf("read three %ldx%ld input frames.\n", origWidth, origHeight);

    if (NULL == (data_in = (float *) malloc(sizeof(float) * width * width * 3))) {
      fprintf(stderr, "allocation error. not enough memory?\n");
      exit(EXIT_FAILURE);
    }

    // CFA Mask indicating which color each sensor pixel has
    unsigned char* mask = (unsigned char *) malloc(width * height * sizeof(unsigned char));

    for (unsigned y = 0; y < height; y++) {
      for (unsigned x = 0; x < width; x++) {
        unsigned p = y * width + x;
        if (
          x + y >= origWidth - 1 &&                   // NW boundary
          y > x - origWidth - 1 &&                    // NE boundary
          x + y < origWidth + 2 * origHeight - 1 &&   // SE boundary
          x > y - origWidth                           // SW boundary
        ) {
          if (y % 2 == 0) {
            mask[p] = GREENPOSITION;
          }
          else {
            if ((x + y - 1) % 4 == 0 || (x + y - 1) % 4 == 1) {
              mask[p] = REDPOSITION;
            }
            else {
              mask[p] = BLUEPOSITION;
            }
          }
        }
        else {
          mask[p] = BLANK;
        }
      }
    }

    for (unsigned long i = 0; i < width * width * 3; i++) {
      data_in[i] = 0;
    }
    for (unsigned long i = 0; i < (unsigned long)origWidth * origHeight; i++) {
      if (origWidth > origHeight) {
        // Landscape
        //
        // B........G
        // ..........
        // ..........
        // G........R
        //
        landscape = true;
        unsigned long x = i % origWidth + (unsigned long)(i / origWidth); // points to the first pixel in a pair
        unsigned long y = (origWidth - i % origWidth - 1) + (i / origWidth);
        unsigned long rix = y * width + x;
        unsigned long gix = y * width + x + width * width;
        unsigned long bix = y * width + x + width * width * 2;
        // data_in[y * width + x] = frame0[i];
        if (y % 2 == 0) {
          data_in[gix] = frame1[gix];
          data_in[gix + 1] = frame1[gix + 1];
        }
        else {
          if ((x + y - 1) % 4 == 0 || (x + y - 1) % 4 == 1) {
            data_in[rix] = frame0[rix];
            data_in[rix + 1] = frame0[rix + 1];
          }
          else {
            data_in[bix] = frame2[bix];
            data_in[bix + 1] = frame2[bix + 1];
          }
        }
      }
      else {
        // Portrait 270° CW
        //
        //  G.....R
        //  .......
        //  .......
        //  .......
        //  B.....G
        //
        unsigned long x0 = origHeight - 1 + i % origWidth - (unsigned long)(i / origWidth);
        unsigned long x1 = x0 + 1; // the second frame (fn == 1) is shifted 1px to the right
        unsigned long y = i % origWidth + (unsigned long)(i / origWidth);
        data_in[y * width + x0] = frame0[i];
        data_in[y * width + x1] = frame1[i];
        data_in[y * width + x0 + width * width] = frame0[i];
        data_in[y * width + x1 + width * width] = frame1[i];
        data_in[y * width + x0 + width * width * 2] = frame0[i];
        data_in[y * width + x1 + width * width * 2] = frame1[i];
      }
    }
  } // merged CFA on input

  else { // Raw EXR Bayer frames
    printf("input file 1: %s\n", args.input_file_0);
    printf("input file 2: %s\n", args.input_file_1);

    /* TIFF 16-bit grayscale -> float input */
    if (NULL == (frame0 = read_tiff_gray16_f32(args.input_file_0, &nx0, &ny0, &description))) {
      fprintf(stderr, "error while reading from %s\n", args.input_file_0);
      exit(EXIT_FAILURE);
    }
    if (NULL == (frame1 = read_tiff_gray16_f32(args.input_file_1, &nx1, &ny1, &description))) {
      fprintf(stderr, "error while reading from %s\n", args.input_file_1);
      exit(EXIT_FAILURE);
    }
    if (nx0 != nx1 or ny0 != ny1) {
      fprintf(stderr, "Input frames must have identical size. Got %ldx%ld vs. %ldx%ld\n", nx0, ny0, nx1, ny1);
      exit(EXIT_FAILURE);
    }
    origWidth = nx0;
    origHeight = ny0;
    width = origWidth + origHeight;

    if (NULL == (data_in = (float *) malloc(sizeof(float) * width * width * 3))) {
      fprintf(stderr, "allocation error. not enough memory?\n");
      exit(EXIT_FAILURE);
    }

    for (unsigned long i = 0; i < width * width * 3; i++) {
      data_in[i] = 0;
    }
    for (unsigned long i = 0; i < (unsigned long)origWidth * origHeight; i++) {
      if (origWidth > origHeight) {
        // Landscape
        //
        // B........G
        // ..........
        // ..........
        // G........R
        //
        landscape = true;
        unsigned long x0 = i % origWidth + (unsigned long)(i / origWidth);
        unsigned long x1 = x0 + 1; // the second frame (fn == 1) is shifted 1px to the right
        unsigned long y = (origWidth - i % origWidth - 1) + (i / origWidth);
        data_in[y * width + x0] = frame0[i];
        data_in[y * width + x1] = frame1[i];
        data_in[y * width + x0 + width * width] = frame0[i];
        data_in[y * width + x1 + width * width] = frame1[i];
        data_in[y * width + x0 + width * width * 2] = frame0[i];
        data_in[y * width + x1 + width * width * 2] = frame1[i];
      }
      else {
        // Portrait 270° CW
        //
        //  G.....R
        //  .......
        //  .......
        //  .......
        //  B.....G
        //
        unsigned long x0 = origHeight - 1 + i % origWidth - (unsigned long)(i / origWidth);
        unsigned long x1 = x0 + 1; // the second frame (fn == 1) is shifted 1px to the right
        unsigned long y = i % origWidth + (unsigned long)(i / origWidth);
        data_in[y * width + x0] = frame0[i];
        data_in[y * width + x1] = frame1[i];
        data_in[y * width + x0 + width * width] = frame0[i];
        data_in[y * width + x1 + width * width] = frame1[i];
        data_in[y * width + x0 + width * width * 2] = frame0[i];
        data_in[y * width + x1 + width * width * 2] = frame1[i];
      }
    }
  } // Raw EXR Bayer frames

  fprintf(stderr, "merged input frames\n");

  write_tiff_rgb_f32("input-merged.tif", data_in, width, width);

  if (NULL == (data_out = (float *) malloc(sizeof(float) * width * width * 3))) {
    fprintf(stderr, "allocation error. not enough memory?\n");
    exit(EXIT_FAILURE);
  }

  /* process */
  ssd_demosaic_chain(
    data_in,
    data_in + width * width,
    data_in + 2 * width * width,
    data_out,
    data_out + width * width,
    data_out + 2 * width * width,
    (int) width,
    (int) width,
    landscape ? origWidth : origHeight,
    landscape ? origHeight : origWidth
  );

  /* limit to 0-65535 */
  out_ptr = data_out;
  end_ptr = out_ptr + 3 * width * width;
  while (out_ptr < end_ptr) {
    if ( 0 > *out_ptr)
      *out_ptr = 0;
    if ( 65535 < *out_ptr)
      *out_ptr = 65535;
    out_ptr++;
  }

  fprintf(stderr, "writing output to %s\n", args.output_file);
  write_tiff_rgb_f32(args.output_file, data_out, width, width);

  free(data_in);
  free(data_out);

  exit(EXIT_SUCCESS);

} // run_ssd()

