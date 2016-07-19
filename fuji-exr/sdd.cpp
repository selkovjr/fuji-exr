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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <error.h>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <iomanip>

#include "sdd_args.h"
#include "termcolor.h"
#include "cfa_mask.h"
#include "libdemosaic.h"
#include "io_tiff.h"
#include "write_tiff.h"
#include "libAuxiliary.h"

using namespace std;
using namespace termcolor;

#define DIAG 1.4142136
#define DIAG12 2.236 // sqrt(5)

void run_sdd (struct argp_state* state) {
  PARSE_ARGS_SDD;

  clock_t start_time, end_time;
  double elapsed;

  size_t nx0 = 0, ny0 = 0;
  size_t nx1 = 0, ny1 = 0;
  size_t nx2 = 0, ny2 = 0;
  char *description;
  unsigned long cfaWidth;
  unsigned long cfaHeight;
  unsigned long width, height;
  float *frame0, *frame1, *frame2; // Bayer EXR frames or TCA-corrected R, G, B
  float *data_in, *data_out, *data_rot;
  float *out_ptr, *end_ptr;
  ushort *u_data_out;
  bool landscape = false;

  if (args.merged_cfa) {
    printf("geometry: %s\n", args.geometry);
    printf("red input file: %s\n", args.input_file_0);
    printf("green input file: %s\n", args.input_file_1);
    printf("blue input file: %s\n", args.input_file_2);

    if (sscanf(args.geometry, "%lux%lu", &cfaWidth, &cfaHeight) < 0) {
      fprintf(stderr, "error parsing image geometry '%s'\n", args.geometry);
      exit(EXIT_FAILURE);
    }
    width = height = cfaWidth + cfaHeight;

    start_time = clock();
    {
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
      printf("read three %ldx%ld input color planes (rotated %ldx%ld).\n", width, height, cfaWidth, cfaHeight);
    }
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    fprintf(stderr, "%6.3f seconds to read input\n", elapsed);

    start_time = clock();
    {
      if (nx0 != nx1 or nx0 != nx2 or nx1 != nx2 or ny0 != ny1 or ny0 != ny2 or ny1 != ny2) {
        fprintf(stderr, "Input color planes must have identical size. Got %ldx%ld, %ldx%ld, %ldx%ld\n", nx0, ny0, nx1, ny1, nx2, ny2);
        exit(EXIT_FAILURE);
      }
      if (nx0 != width) {
        fprintf(stderr, "Stated image geometry (%ldx%ld) does not fit input color planes (%ldx%ld)\n", cfaWidth, cfaHeight, nx0, ny0);
        exit(EXIT_FAILURE);
      }

      if (NULL == (data_in = (float *) malloc(sizeof(float) * width * width * 3))) {
        fprintf(stderr, "allocation error. not enough memory?\n");
        exit(EXIT_FAILURE);
      }
      if (NULL == (data_out = (float *) malloc(sizeof(float) * width * width * 3))) {
        fprintf(stderr, "allocation error. not enough memory?\n");
        exit(EXIT_FAILURE);
      }

      for (unsigned long i = 0; i < width * width * 3; i++) {
        data_in[i] = 0;
      }
    }
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    fprintf(stderr, "%6.3f seconds to allocate and zero-set memory\n", elapsed);

    start_time = clock();
    for (unsigned long i = 0; i < (unsigned long)cfaWidth * cfaHeight; i++) {
      if (cfaWidth > cfaHeight) {
        // Landscape
        //
        // B........G
        // ..........
        // ..........
        // G........R
        //
        landscape = true;
        unsigned long x = i % cfaWidth + (unsigned long)(i / cfaWidth); // points to the first pixel in a pair
        unsigned long y = (cfaWidth - i % cfaWidth - 1) + (i / cfaWidth);
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
        unsigned long x0 = cfaHeight - 1 + i % cfaWidth - (unsigned long)(i / cfaWidth);
        unsigned long x1 = x0 + 1; // the second frame (fn == 1) is shifted 1px to the right
        unsigned long y = i % cfaWidth + (unsigned long)(i / cfaWidth);
        data_in[y * width + x0] = frame0[i];
        data_in[y * width + x1] = frame1[i];
        data_in[y * width + x0 + width * width] = frame0[i];
        data_in[y * width + x1 + width * width] = frame1[i];
        data_in[y * width + x0 + width * width * 2] = frame0[i];
        data_in[y * width + x1 + width * width * 2] = frame1[i];
      }
    }
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    fprintf(stderr, "%6.3f seconds to merge input color planes\n", elapsed);
    // write_tiff_rgb_f32("input-merged.tif", data_in, width, width);
  } // merged CFA on input

  else { // Raw EXR Bayer frames
    start_time = clock();
    {
      /* TIFF 16-bit grayscale -> float input */
      fprintf(stderr, "input file 0: %s\n", args.input_file_0);
      if (NULL == (frame0 = read_tiff_gray16_f32(args.input_file_0, &nx0, &ny0, &description))) {
        fprintf(stderr, "error while reading from %s\n", args.input_file_0);
        exit(EXIT_FAILURE);
      }

      fprintf(stderr, "input file 1: %s\n", args.input_file_1);
      if (NULL == (frame1 = read_tiff_gray16_f32(args.input_file_1, &nx1, &ny1, &description))) {
        fprintf(stderr, "error while reading from %s\n", args.input_file_1);
        exit(EXIT_FAILURE);
      }
    }
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    fprintf(stderr, "%6.3f seconds to read input\n", elapsed);

    start_time = clock();
    {
      if (nx0 != nx1 or ny0 != ny1) {
        fprintf(stderr, "Input frames must have identical size. Got %ldx%ld vs. %ldx%ld\n", nx0, ny0, nx1, ny1);
        exit(EXIT_FAILURE);
      }
      cfaWidth = nx0;
      cfaHeight = ny0;
      width = height = cfaWidth + cfaHeight;
      landscape = cfaWidth > cfaHeight ? true : false;

      if (NULL == (data_in = (float *) malloc(sizeof(float) * width * width * 3))) {
        fprintf(stderr, "allocation error. not enough memory?\n");
        exit(EXIT_FAILURE);
      }

      for (unsigned long i = 0; i < width * width * 3; i++) {
        data_in[i] = 0;
      }
    }
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    fprintf(stderr, "%6.3f seconds to allocate and zero-set memory\n", elapsed);

    start_time = clock();
    for (unsigned long i = 0; i < (unsigned long)cfaWidth * cfaHeight; i++) {
      if (cfaWidth > cfaHeight) {
        // Landscape
        //
        // B........G
        // ..........
        // ..........
        // G........R
        //
        landscape = true;
        unsigned long x0 = i % cfaWidth + (unsigned long)(i / cfaWidth);
        unsigned long x1 = x0 + 1; // the second frame (fn == 1) is shifted 1px to the right
        unsigned long y = (cfaWidth - i % cfaWidth - 1) + (i / cfaWidth);
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
        unsigned long x0 = cfaHeight - 1 + i % cfaWidth - (unsigned long)(i / cfaWidth);
        unsigned long x1 = x0 + 1; // the second frame (fn == 1) is shifted 1px to the right
        unsigned long y = i % cfaWidth + (unsigned long)(i / cfaWidth);
        data_in[y * width + x0] = frame0[i];
        data_in[y * width + x1] = frame1[i];
        data_in[y * width + x0 + width * width] = frame0[i];
        data_in[y * width + x1 + width * width] = frame1[i];
        data_in[y * width + x0 + width * width * 2] = frame0[i];
        data_in[y * width + x1 + width * width * 2] = frame1[i];
      }
    }
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    fprintf(stderr, "%6.3f seconds to merge input frames\n", elapsed);

    //write_tiff_rgb_f32("input-merged.tif", data_in, width, width);
  } // Raw EXR Bayer frames


  if (NULL == (data_out = (float *) malloc(sizeof(float) * width * width * 3))) {
    fprintf(stderr, "allocation error. not enough memory?\n");
    exit(EXIT_FAILURE);
  }

  start_time = clock();
  unsigned char *mask = exr_cfa_mask(width, width, cfaWidth, cfaHeight);
  end_time = clock();
  elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
  fprintf(stderr, "%6.3f seconds to compute CFA mask\n", elapsed);

  /* process */
  start_time = clock();
  sdd_demosaic_chain (
    data_in,
    data_in + width * width,
    data_in + 2 * width * width,
    data_out,
    data_out + width * width,
    data_out + 2 * width * width,
    (int) width,
    (int) width,
    landscape ? cfaWidth : cfaHeight,
    landscape ? cfaHeight : cfaWidth,
    mask
  );
  end_time = clock();
  elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
  fprintf(stderr, "%6.3f seconds to complete debayering\n", elapsed);

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

  // write_tiff_rgb_f32("result.tiff", data_out, width, width);

  // ---------------------------------------------------------------------------
  // Rotate the interpolated result 45°
  //
  start_time = clock();

  int i, row, col;
  double step = sqrt(0.5); // Horizontal or vertical CFA step projected onto
                           // source-plane axes
  float r, c;              // Y- and X-coords in the source plane
  unsigned ur, uc;         // Y- and X-coords of the nearest source pixel
  float fr, fc;            // Y- and X-distance from (r, c) to nearest pixel
  ushort rotWidth, rotHeight;


  // Inflated (√2) target image co-ordinates
  rotWidth = cfaWidth / step;
  rotHeight = (height - cfaWidth) / step;

  if (NULL == (data_rot = (float *) malloc(sizeof(float) * rotWidth * rotHeight * 3))) {
    fprintf(stderr, "allocation error. not enough memory?\n");
    exit(EXIT_FAILURE);
  }

  // Row and col are co-ordinates in the inflated target image.
  for (row = 0; row < rotHeight; row++) {
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

        // David Coffin's original stencil (chunked configuration)
        //
        //   pix = img + ur * iwidth + uc;
        //   img[row * wide + col][i] =
        //     (/* + */ pix[    0][i]*(1 - fc) + /* E  */ pix[        1][i] * fc) * (1 - fr) +
        //     (/* S */ pix[width][i]*(1 - fc) + /* SE */ pix[width + 1][i] * fc) * fr;
        //
        // Same stencil reformulated for planar configuration
        //
        data_rot[row * rotWidth + col + i * rotWidth * rotHeight] =
          (1 - fr) * (
            (1 - fc) * data_out[ur * width + uc + i * width * height]              // +
            +
                  fc * data_out[ur * width + uc + i * width * height + 1]          // E
          )
          +
          fr * (
            (1 - fc) * data_out[ur * width + uc + i * width * height + width]      // S
            +
                  fc * data_out[ur * width + uc + i * width * height + width + 1]  // SE
          )
          ;
      } // each color plane
    }
  }
  end_time = clock();
  elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
  fprintf(stderr, "%6.3f seconds to rotate\n", elapsed);


  cerr << grey << "writing output to " << white << args.output_file << reset << endl;
  start_time = clock();
  {
    // Not using libtiff to write output because it creates invalid TIFF directories.

    // Convert from planar to chunked
    if (NULL == (u_data_out = (ushort *) malloc(sizeof(float) * rotWidth * rotHeight * 3))) {
      cerr << on_red << "allocation error: not enough memory" << reset << endl;
      exit(EXIT_FAILURE);
    }
    for (unsigned x = 0; x < rotWidth; x++) {
      for (unsigned y = 0; y < rotHeight; y++) {
        unsigned p = y * rotWidth + x;
        u_data_out[y * 3 * rotWidth + x * 3] = data_rot[p];
        u_data_out[y * 3 * rotWidth + x * 3 + 1] = data_rot[p + rotWidth * rotHeight];
        u_data_out[y * 3 * rotWidth + x * 3 + 2] = data_rot[p + 2 * rotWidth * rotHeight];
      }
    }
    write_tiff_img(args.output_file, (unsigned char *)u_data_out, rotWidth, rotHeight, 16, 3, 0 /* chunked */);
  }
  end_time = clock();
  elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
  cerr << yellow << setw(7) << setprecision(2) << elapsed << "s" << white << " writing" << reset << endl;

  delete[] mask;
  free(data_in);
  free(data_out);
  free(data_rot);
  free(u_data_out);

  exit(EXIT_SUCCESS);

} // run_sdd()

