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
  size_t nx = 0, ny = 0;
  char *description;
  int origWidth;
  int origHeight;
  float *data_in, *data_out;
  float *out_ptr, *end_ptr;

  /* version info */
  if (2 <= argc && 0 == strcmp("-v", argv[1])) {
    fprintf(stdout, "%s version " __DATE__ "\n", argv[0]);
    return EXIT_SUCCESS;
  }

  /* sanity check */
  if (3 != argc) {
    fprintf(stderr, "usage : %s input.tiff output.tiff\n", argv[0]);
    return EXIT_FAILURE;
  }

  /* TIFF 16-bit grayscale -> float input */
  if (NULL == (data_in = read_tiff_gray16_f32(argv[1], &nx, &ny, &description))) {
    fprintf(stderr, "error while reading from %s\n", argv[1]);
    return EXIT_FAILURE;
  }

  if (sscanf(description, "width = %d, height = %d", &origWidth, &origHeight) != 2) {
    fprintf(stderr, "image description (%s) does not contain a size expression (width = nnnn, height = nnnn)\n", description);
    return EXIT_FAILURE;
  }

  if (NULL == (data_out = (float *) malloc(sizeof(float) * nx * ny * 4))) {
    fprintf(stderr, "allocation error. not enough memory?\n");
    free(data_in);
    return EXIT_FAILURE;
  }


  /* process */
  ssd_demosaic_chain(
    data_in,
    data_in + nx * ny,
    data_in + 2 * nx * ny,
    data_out,
    data_out + nx * ny,
    data_out + 2 * nx * ny,
    (int) nx,
    (int) ny,
    origWidth,
    origHeight
  );

  /* limit to 0-65535 */
  out_ptr = data_out;
  end_ptr = out_ptr + 3 * nx * ny;
  while (out_ptr < end_ptr) {
    if ( 0 > *out_ptr)
      *out_ptr = 0;
    if ( 65535 < *out_ptr)
      *out_ptr = 65535;
    out_ptr++;
  }

  write_tiff_rgb_f32(argv[2], data_out, nx, ny);

  free(data_in);
  free(data_out);

  return EXIT_SUCCESS;
}

