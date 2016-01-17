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


int main(int argc, char **argv) {
  size_t nx0 = 0, ny0 = 0;
  size_t nx1 = 0, ny1 = 0;
  char *description;
  unsigned long origWidth;
  unsigned long origHeight;
  unsigned long width;
  float *frame0, *frame1;
  float *data_in, *data_out;
  float *out_ptr, *end_ptr;
  bool landscape = false;

  /* version info */
  if (2 <= argc && 0 == strcmp("-v", argv[1])) {
    fprintf(stdout, "%s version " __DATE__ "\n", argv[0]);
    return EXIT_SUCCESS;
  }

  /* sanity check */
  if (4 != argc) {
    fprintf(stderr, "usage : %s input_0.tiff input_1.tiff output.tiff\n", argv[0]);
    return EXIT_FAILURE;
  }

  /* TIFF 16-bit grayscale -> float input */
  if (NULL == (frame0 = read_tiff_gray16_f32(argv[1], &nx0, &ny0, &description))) {
    fprintf(stderr, "error while reading from %s\n", argv[1]);
    return EXIT_FAILURE;
  }
  if (NULL == (frame1 = read_tiff_gray16_f32(argv[1], &nx1, &ny1, &description))) {
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
  if (NULL == (data_out = (float *) malloc(sizeof(float) * width * width * 3))) {
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
  fprintf(stderr, "merged input frames\n");

  write_tiff_rgb_f32("input-merged.tif", data_in, width, width);

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

