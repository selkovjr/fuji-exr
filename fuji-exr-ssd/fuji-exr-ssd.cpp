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
#include <stdlib.h>
#include <string.h>

#include "libdemosaic.h"
#include "io_tiff.h"
#include "tiffio.h"

/**
 * @file   fuji-exr-ssd.cpp
 * @brief  Main executable file
 *
 *
 *
 * @author Antoni Buades <toni.buades@uib.es>
 */

#define BLANK 0
#define REDPOSITION 1
#define GREENPOSITION 2
#define BLUEPOSITION 3

int main(int argc, char **argv) {
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
  bool corrected = false; // TCA-corrected input

  /* version info */
  if (argc >= 2 && strcmp("-v", argv[1]) == 0) {
    fprintf(stdout, "%s version " __DATE__ "\n", argv[0]);
    return EXIT_SUCCESS;
  }

  if (argc >= 2 && strcmp("-c", argv[1]) == 0) {
    corrected = true;
  }

  /* sanity check */
  if ((corrected and argc != 8) or (!corrected and argc != 4)) {
    fprintf(stderr, "usage: %s input_frame_0.tiff input_frame_1.tiff output.tiff\n", argv[0]);
    fprintf(stderr, "       %s -c width height tca-corrected-r.tiff g.tiff tca-corrected-b.tiff output.tiff\n", argv[0]);
    return EXIT_FAILURE;
  }


  if (corrected) {
    origWidth = atoi(argv[2]);
    origHeight = atoi(argv[3]);
    width = height = origWidth + origHeight;

    /* TIFF 16-bit grayscale -> float input */
    if (NULL == (frame0 = read_tiff_gray16_f32(argv[4], &nx0, &ny0, &description))) {
      fprintf(stderr, "error while reading from %s\n", argv[4]);
      return EXIT_FAILURE;
    }
    if (NULL == (frame1 = read_tiff_gray16_f32(argv[5], &nx1, &ny1, &description))) {
      fprintf(stderr, "error while reading from %s\n", argv[5]);
      return EXIT_FAILURE;
    }
    if (NULL == (frame2 = read_tiff_gray16_f32(argv[6], &nx2, &ny2, &description))) {
      fprintf(stderr, "error while reading from %s\n", argv[6]);
      return EXIT_FAILURE;
    }
    if (nx0 != nx1 or nx0 != nx2 or nx1 != nx2 or ny0 != ny1 or ny0 != ny2 or ny1 != ny2) {
      fprintf(stderr, "Input frames must have identical size. Got %ldx%ld, %ldx%ld, %ldx%ld\n", nx0, ny0, nx1, ny1, nx2, ny2);
      return EXIT_FAILURE;
    }
    if (nx0 != width) {
      fprintf(stderr, "Stated image geometry (%ldx%ld) does not fit input frames (%ldx%ld)\n", origWidth, origHeight, nx0, ny0);
      return EXIT_FAILURE;
    }
    printf("read three %ldx%ld input frames.\n", origWidth, origHeight);

    if (NULL == (data_in = (float *) malloc(sizeof(float) * width * width * 3))) {
      fprintf(stderr, "allocation error. not enough memory?\n");
      return EXIT_FAILURE;
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
  } // TCA-corrected

  else { // Raw EXR Bayer frames
    /* TIFF 16-bit grayscale -> float input */
    if (NULL == (frame0 = read_tiff_gray16_f32(argv[1], &nx0, &ny0, &description))) {
      fprintf(stderr, "error while reading from %s\n", argv[1]);
      return EXIT_FAILURE;
    }
    if (NULL == (frame1 = read_tiff_gray16_f32(argv[2], &nx1, &ny1, &description))) {
      fprintf(stderr, "error while reading from %s\n", argv[2]);
      return EXIT_FAILURE;
    }
    if (nx0 != nx1 or ny0 != ny1) {
      fprintf(stderr, "Input frames must have identical size. Got %ldx%ld vs. %ldx%ld\n", nx0, ny0, nx1, ny1);
      return EXIT_FAILURE;
    }
    origWidth = nx0;
    origHeight = ny0;
    width = origWidth + origHeight;

    if (NULL == (data_in = (float *) malloc(sizeof(float) * width * width * 3))) {
      fprintf(stderr, "allocation error. not enough memory?\n");
      return EXIT_FAILURE;
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
  }
  fprintf(stderr, "merged input frames\n");

  write_tiff_rgb_f32("input-merged.tif", data_in, width, width);

  if (NULL == (data_out = (float *) malloc(sizeof(float) * width * width * 3))) {
    fprintf(stderr, "allocation error. not enough memory?\n");
    return EXIT_FAILURE;
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

  fprintf(stderr, "writing output to %s\n", argv[3]);
  write_tiff_rgb_f32(argv[3], data_out, width, width);

  free(data_in);
  free(data_out);

  return EXIT_SUCCESS;
}

