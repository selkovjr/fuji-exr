/*
* Copyright (c) 2009-2011, A. Buades <toni.buades@uib.es>
* All rights reserved.
*
*
* Patent warning:
*
* This file implements algorithms possibly linked to the patents
*
* # J. Hamilton Jr and J. Adams Jr, “Adaptive color plan interpolation
* in single sensor color electronic camera,” 1997, US Patent 5,629,734.
*
* # D. Cok, “Signal processing method and apparatus for producing
* interpolated chrominance values in a sampled color image signal”,
* 1987, US Patent 4,642,678.
*
* # A. Buades, T. Coll and J.M. Morel, Image data processing method by
* reducing image noise, and camera integrating means for implementing
* said method, EP Patent 1,749,278 (Feb. 7, 2007).
*
* This file is made available for the exclusive aim of serving as
* scientific tool to verify the soundness and completeness of the
* algorithm description. Compilation, execution and redistribution
* of this file may violate patents rights in certain countries.
* The situation being different for every country and changing
* over time, it is your responsibility to determine which patent
* rights restrictions apply to you before you compile, use,
* modify, or redistribute this file. A patent lawyer is qualified
* to make this determination.
* If and only if they don't conflict with any patent terms, you
* can benefit from the following license terms attached to this
* file.
*
* License:
*
* This program is provided for scientific and educational only:
* you can use and/or modify it for these purposes, but you are
* not allowed to redistribute this work or derivative works in
* source or executable form. A license must be obtained from the
* patent right holders for any other use.
*
*
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
    fprintf(stderr, "usage : %s input.tiff output.tiff\n",
        argv[0]);
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

  if (NULL == (data_out = (float *) malloc(sizeof(float) * nx * ny * 3))) {
    fprintf(stderr, "allocation error. not enough memory?\n");
    free(data_in);
    return EXIT_FAILURE;
  }


  /* process */
  ssd_demosaic_chain(
    data_in,
    data_out,
    data_out + nx * ny,
    data_out + 2 * nx * ny,
    (int) nx,
    (int) ny,
    origWidth,
    origHeight
  );

  printf("demosaic_chain() -> %f\n", data_out[8990]);

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

  /* copy alpha channel */
  memcpy(data_out + 3 * nx * ny, data_in + 3 * nx * ny, nx * ny * sizeof(float));

  /* TIFF RGB float->8bit output */
  write_tiff_rgb_f32(argv[2], data_out, nx, ny);

  free(data_in);
  free(data_out);

  return EXIT_SUCCESS;
}

