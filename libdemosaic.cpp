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

#include <algorithm>
#include "libdemosaic.h"

#include "io_tiff.h"
#include "tiffio.h"


#define BLANK 0
#define GREENPOSITION 1
#define REDPOSITION 2
#define BLUEPOSITION 3

#define DIAG 1.4142136
#define DIAG12 2.236 // sqrt(5)

#define DUMP_STAGES
#undef DEBUG_GREEN

/**
 * @file   libdemosaic.cpp
 * @brief  Demosaicking functions: Hamilton-Adams algorithm, NLmeans-based demosaicking, Chromatic components filtering
 *
 *
 * @author Antoni Buades <toni.buades@uib.es>
 */


/**
 * \brief  Classical Adams-Hamilton demosaicking algorithm
 *
 *  The green channel is interpolated directionally depending on the green first and red and blue second directional derivatives.
 *  The red and blue differences with the green channel are interpolated bilinearly.
 *
 * @param[in]  ired, igreen, iblue  original cfa image
 * @param[out] ored, ogreen, oblue  demosaicked output
 * @param[in]  threshold value to consider horizontal and vertical variations equivalent and average both estimates
 * @param[in]  width, height size of the image
 *
 */

void adams_hamilton(
  float threshold,
  float *input,
  float *ored,
  float *ogreen,
  float *oblue,
  int width,
  int height,
  int origWidth,
  int origHeight
) {
#ifdef DUMP_STAGES
  float *data_out;
  if (NULL == (data_out = (float *) malloc(sizeof(float) * width * height * 3))) {
    fprintf(stderr, "allocation error. not enough memory?\n");
    exit(0);
  }
#endif

  wxCopy(input, ored, width * height);
  wxCopy(input, ogreen, width * height);
  wxCopy(input, oblue, width * height);

  // CFA Mask indicating which color each sensor pixel has
  unsigned char* mask = (unsigned char *) malloc(width * height * sizeof(unsigned char));

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int p = y * width + x;
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

  // Interpolate the green channel in the 4-pixel-wide boundary by inverse
  // distance weighting.
  for (int x = 0; x < width; x++) {
    for (int y = 0; y < height; y++) {
      int p = y * width + x;
      if (
        (
         mask[p] == BLUEPOSITION ||
         mask[p] == REDPOSITION
        )
        &&
        (
         x + y < origWidth + 3 ||                    // NW boundary
         x >= y + origWidth - 3 ||                   // NE boundary
         x + y >= origWidth + 2 * origHeight - 5 ||  // SE boundary
         y >= x + origWidth - 4                      // SW boundary
        )
      ) {
        float avg = 0;
        float weight = 0;
        int
          n = (y - 1) * width + x,
          s = (y + 1) * width + x,
          ne = (y - 1) * width + x + 1,
          se = (y + 1) * width + x + 1,
          sw = (y + 1) * width + x - 1,
          nw = (y - 1) * width + x - 1;

        // East and west corners are special cases because
        // in the east corner the [ne] and [se] pixels are undefined,
        // and in the west corner the [nw] and [sw] pixels are undefined
        // and cannot be tested.
        if (x == 0) { // West corner red
          ogreen[p] = (ogreen[ne] + ogreen[se]) / 2;
        }
        else if (x == width - 1) { // East corner blue
          ogreen[p] = (ogreen[nw] + ogreen[sw]) / 2;
        }
        else {
          // Each non-green pixel is surrounded by 6 green pixels,
          // except on the boundaries.
          //
          // G G G
          //   *
          // G G G
          //
          if (mask[nw] != BLANK) {
            avg += ogreen[nw] / DIAG;
            weight += 1 / DIAG;
          }
          if (mask[n] != BLANK) {
            avg += ogreen[n];
            weight += 1;
          }
          if (mask[ne] != BLANK) {
            avg += ogreen[ne] / DIAG;
            weight += 1 / DIAG;
          }

          if (mask[sw] != BLANK) {
            avg += ogreen[sw] / DIAG;
            weight += 1 / DIAG;
          }
          if (mask[s] != BLANK) {
            avg += ogreen[s];
            weight += 1;
          }
          if (mask[se] != BLANK) {
            avg += ogreen[se] / DIAG;
            weight += 1 / DIAG;
          }
          ogreen[p] = avg / weight;
        }
      }
    }
  }

  // Interpolate the green by Adams-Hamilton algorithm inside the image.


  // **************************** G R E E N **************************************
  // First interpolate the green channel directionally, using adaptive color plane
  // interpolation.

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int p = y * width + x;
      if (
        mask[p] != GREENPOSITION &&
        x + y >= origWidth - 1 + 4 &&                   // NW boundary
        y > x - origWidth - 1 + 4 &&                    // NE boundary
        x + y < origWidth + 2 * origHeight - 1 - 4 &&   // SE boundary
        x > y - origWidth + 4                           // SW boundary
      ) {
        int
          n = (y - 1) * width + x,
          n2 = (y - 2) * width + x,
          s = (y + 1) * width + x,
          s2 = (y + 2) * width + x,
          ne = (y - 1) * width + x + 1,
          se = (y + 1) * width + x + 1,
          sw = (y + 1) * width + x - 1,
          nw = (y - 1) * width + x - 1,
          ne2 = (y - 2) * width + x + 2,
          se2 = (y + 2) * width + x + 2,
          sw2 = (y + 2) * width + x - 2,
          nw2 = (y - 2) * width + x - 2;

        // Compute gradients in the green channel.
        //
        // N - S
        float gn = fabsf(ogreen[n] - ogreen[s]); // diagonal distance is longer

        // NW - SE
        float gnw = fabsf(ogreen[nw] - ogreen[se]) / DIAG; // diagonal distance is longer

        // NE - SW
        float gne = fabsf(ogreen[ne] - ogreen[sw]) / DIAG;

        float d2nw, d2n, d2ne;

        // Compute diagonal second derivatives in the same channel as current
        // pixel.

        if (mask[p] == REDPOSITION) {
          //
          //     NW          N           N          NE
          //               even         odd
          // R . . . .   . R . . .   . . . R .   . . . . R
          // . . . . .   . . . . .   . . . . .   . . . . .
          // . . * . .   . . * . .   . . * . .   . . * . .
          // . . . . .   . . . . .   . . . . .   . . . . .
          // . . . . R   . R . . .   . . . R .   R . . . .
          //
          d2nw = (2.0 * ored[p] - ored[nw2] - ored[se2]) / 8.0; // (2 * DIAG) squared
          d2ne = (2.0 * ored[p] - ored[ne2] - ored[sw2]) / 8.0;
          // almost vertical
          if (x % 2) {
            d2n  = (2.0 * ored[p] - ored[n2 - 1] - ored[s2 - 1]) / 5.0; // (DIAG12 squared
          }
          else {
            d2n  = (2.0 * ored[p] - ored[n2 - 1] - ored[s2 - 1]) / 5.0; // (DIAG12 squared
          }
        }
        else { // BLUEPOSITION
          //
          //     NW          N          NE
          //
          // B . . . .   . . R . .   . . . . B
          // . . . . .   . . . . .   . . . . .
          // . . * . .   . . * . .   . . * . .
          // . . . . .   . . . . .   . . . . .
          // . . . . B   . . R . .   B . . . .
          //
          d2nw = (2.0 * oblue[p] - oblue[nw2] - oblue[se2]) / 8.0;
          d2ne = (2.0 * oblue[p] - oblue[ne2] - oblue[sw2]) / 8.0;
          if (x % 2) {
            d2n  = (2.0 * oblue[p] - oblue[n2 - 1] - oblue[s2 - 1]) / 5.0; // (DIAG12 squared
          }
          else {
            d2n  = (2.0 * oblue[p] - oblue[n2 - 1] - oblue[s2 - 1]) / 5.0; // (DIAG12 squared
          }
        }

        // Add second differences to gradients
        gnw += fabsf(d2nw);
        gn  += fabsf(d2n);
        gne += fabsf(d2ne);

#ifdef DEBUG_GREEN
        int
          e = p + 1,
          w = p - 1;

        if (mask[p] == BLUEPOSITION) {
          printf("\x1b[34m--------------------------------------------------------------\x1b[0m\n");
          printf("\x1b[34mblue\x1b[0m, %s (%d, %d)\n", x % 2 ? "odd" : "even", x, y);
          if (x % 2 == 0) {
            printf("\x1b[34m%6.0f%6.0f\x1b[31m%6.0f%6.0f\x1b[34m%6.0f\x1b[0m\n",   input[nw2], input[nw2 + 1], input[n2], input[ne2 - 1], input[ne2]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    input[nw],      input[n],  input[ne]);
            printf( "      \x1b[31m%6.0f\x1b[1m\x1b[34m%6.0f\x1b[0m\x1b[34m%6.0f\x1b[0m\n",                                    input[w],       input[p],  input[e]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    input[sw],      input[s],  input[se]);
            printf("\x1b[34m%6.0f%6.0f\x1b[31m%6.0f%6.0f\x1b[34m%6.0f\x1b[0m\n",   input[sw2], input[sw2 + 1], input[s2], input[se2 - 1], input[se2]);
          }
          else {
            printf("\x1b[34m%6.0f\x1b[31m%6.0f%6.0f\x1b[34m%6.0f%6.0f\x1b[0m\n",   input[nw2], input[nw2 + 1], input[n2], input[ne2 - 1], input[ne2]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    input[nw],      input[n],  input[ne]);
            printf( "      \x1b[34m%6.0f\x1b[1m%6.0f\x1b[0m\x1b[31m%6.0f\x1b[0m\n",                                    input[w],       input[p],  input[e]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    input[sw],      input[s],  input[se]);
            printf("\x1b[34m%6.0f\x1b[31m%6.0f%6.0f\x1b[34m%6.0f%6.0f\x1b[0m\n",   input[sw2], input[sw2 + 1], input[s2], input[se2 - 1], input[se2]);
          }
        }
        else { // red
          printf("\x1b[31m--------------------------------------------------------------\x1b[0m\n");
          printf("\x1b[31mred\x1b[0m, %s (%d, %d)\n", x % 2 ? "odd" : "even", x, y);
          if (x % 2 == 0) {
            printf("\x1b[31m%6.0f%6.0f\x1b[34m%6.0f%6.0f\x1b[31m%6.0f\x1b[0m\n",   input[nw2], input[nw2 + 1], input[n2], input[ne2 - 1], input[ne2]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    input[nw],      input[n],  input[ne]);
            printf( "      \x1b[34m%6.0f\x1b[1m\x1b[31m%6.0f\x1b[0m\x1b[31m%6.0f\x1b[0m\n",                                    input[w],       input[p],  input[e]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    input[sw],      input[s],  input[se]);
            printf("\x1b[31m%6.0f%6.0f\x1b[34m%6.0f%6.0f\x1b[31m%6.0f\x1b[0m\n",   input[sw2], input[sw2 + 1], input[s2], input[se2 - 1], input[se2]);
          }
          else {
            printf("\x1b[31m%6.0f\x1b[34m%6.0f%6.0f\x1b[31m%6.0f%6.0f\x1b[0m\n",   input[nw2], input[nw2 + 1], input[n2], input[ne2 - 1], input[ne2]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    input[nw],      input[n],  input[ne]);
            printf( "      \x1b[31m%6.0f\x1b[1m%6.0f\x1b[0m\x1b[34m%6.0f\x1b[0m\n",                                    input[w],       input[p],  input[e]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    input[sw],      input[s],  input[se]);
            printf("\x1b[31m%6.0f\x1b[34m%6.0f%6.0f\x1b[31m%6.0f%6.0f\x1b[0m\n",   input[sw2], input[sw2 + 1], input[s2], input[se2 - 1], input[se2]);
          }
        }
#endif

        d2n = 0;
        d2nw = 0;
        d2ne = 0;

        float gmin = fminf(gnw, fminf(gn, gne));

        // If all differences are similar, compute an isotropic average.
        if ( 1 || (
          fabsf(gnw - gmin) < threshold &&
          fabsf(gn - gmin) < threshold &&
          fabsf(gne - gmin) < threshold
        )) {
#ifdef DEBUG_GREEN
          printf(  "grad # (\\ %5.0f, | %5.0f, / %5.0f)\n", gnw, gn, gne);
#endif
          ogreen[p] =
            (
             ogreen[nw] / DIAG + ogreen[n] + ogreen[ne] / DIAG +
             ogreen[se] / DIAG + ogreen[s] + ogreen[sw] / DIAG
            ) / (2 + 4 / DIAG) +
            (d2nw + d2n + d2ne) / 3;
        }
        else if (gmin == gnw) {
          // NW average
#ifdef DEBUG_GREEN
          printf(  "grad \\ (\\ %5.0f, | %5.0f, / %5.0f)\n", gnw, gn, gne);
#endif
          ogreen[p] = (ogreen[nw] + ogreen[se]) / 2.0 + d2nw;
        }
        else if (gmin == gn) {
          // N
#ifdef DEBUG_GREEN
          printf(  "grad | (\\ %5.0f, | %5.0f, / %5.0f)\n", gnw, gn, gne);
#endif
          ogreen[p] = (ogreen[n] + ogreen[s]) / 2.0 + d2n;
        }
        else if (gmin == gne) {
          // NE
#ifdef DEBUG_GREEN
          printf(  "grad / (\\ %5.0f, | %5.0f, / %5.0f)\n", gnw, gn, gne);
#endif
          ogreen[p] = (ogreen[ne] + ogreen[sw]) / 2.0 + d2ne;
        }
#ifdef DEBUG_GREEN
        printf("    \x1b[32mavg: %5.0f\x1b[0m\n", ogreen[p]);
#endif
      }
    }
  }

#ifdef DUMP_STAGES
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int p = y * width + x;
      data_out[p] = ored[p];
      data_out[p + width * height] = 0;
      data_out[p + 2 * width * height] = 0;
    }
  }
  /* TIFF RGB float->8bit output */
  write_tiff_rgb_f32("red-init.tiff", data_out, width, height);

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int p = y * width + x;
      data_out[p] = 0;
      data_out[p + width * height] = ogreen[p];
      data_out[p + 2 * width * height] = 0;
    }
  }
  /* TIFF RGB float->8bit output */
  write_tiff_rgb_f32("green-interpolated.tiff", data_out, width, height);

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int p = y * width + x;
      data_out[p] = 0;
      data_out[p + width * height] = 0;
      data_out[p + 2 * width * height] = oblue[p];
    }
  }
  /* TIFF RGB float->8bit output */
  write_tiff_rgb_f32("blue-init.tiff", data_out, width, height);
#endif

  // compute the bilinear on the differences of the red and blue with the already interpolated green
  bilinear_red_blue(ored, ogreen, oblue, width, height, origWidth, origHeight);

  free(mask);
}




/**
 * \brief  Classical bilinear interpolation of red and blue differences with the green channel
 *
 *
 * @param[in]  ored, ogreen, oblue  original cfa image with green interpolated
 * @param[out] ored, ogreen, oblue  demosaicked output
 * @param[in]  width, height size of the merged image
 * @param[in]  origWidth, origHeight original image size
 *
 */

void bilinear_red_blue(
  float *ored,
  float *ogreen,
  float *oblue,
  int width,
  int height,
  int origWidth,
  int origHeight
) {
#ifdef DUMP_STAGES
  float *data_out;
  if (NULL == (data_out = (float *) malloc(sizeof(float) * width * height * 3))) {
    fprintf(stderr, "allocation error. not enough memory?\n");
    exit(0);
  }
#endif

  // Mask of color per pixel
  unsigned char* mask = (unsigned char *) malloc(width*height*sizeof(unsigned char));

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      if (
        x + y >= origWidth - 1 &&                   // NW boundary
        y > x - origWidth - 1 &&                    // NE boundary
        x + y < origWidth + 2 * origHeight - 1 &&   // SE boundary
        x > y - origWidth                           // SW boundary
      ) {
        int p = y * width + x;
        // The raw image has zeroes in all channels outside the sensor area
        if (y % 2 == 0) {
          mask[p] = GREENPOSITION;
          // ored[p] = 500;
          // ogreen[p] = 65535;
          // oblue[p] = 500;
        }
        else {
          if ((x + y - 1) % 4 == 0 || (x + y - 1) % 4 == 1) {
            mask[p] = REDPOSITION;
            // ored[p] = 65535;
            // ogreen[p] = 500;
            // oblue[p] = 500;
          }
          else {
            mask[p] = BLUEPOSITION;
            // ored[p] = 500;
            // ogreen[p] = 500;
            // oblue[p] = 65535;
          }
        }
      }
      else {
        mask[y * width + x] = BLANK;
      }
    }
  }

  // Compute the differences
  for (int i = 0; i < width * height; i++) {
    ored[i] -= ogreen[i];
    oblue[i] -= ogreen[i];
  }
#ifdef DUMP_STAGES
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int p = y * width + x;
      data_out[p] = ored[p];
      data_out[p + width * height] = 0;
      data_out[p + 2 * width * height] = 0;
    }
  }
  /* TIFF RGB float->8bit output */
  write_tiff_rgb_f32("red-diff.tiff", data_out, width, height);

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int p = y * width + x;
      data_out[p] = 0;
      data_out[p + width * height] = 0;
      data_out[p + 2 * width * height] = oblue[p];
    }
  }
  /* TIFF RGB float->8bit output */
  write_tiff_rgb_f32("blue-diff.tiff", data_out, width, height);
#endif

  // Interpolate blue differences making the average of possible values depending on the CFA location
  for (int x = 0; x < width; x++) {
    for (int y = 0; y < height; y++) {
      int p = y * width + x;
      if (mask[p] != BLUEPOSITION && mask[p] != BLANK) {
        int
          n = (y - 1) * width + x,
          s = (y + 1) * width + x,
          e = p + 1,
          w = p - 1,
          e2 = p + 2,
          w2 = p - 2,
          ne = (y - 1) * width + x + 1,
          se = (y + 1) * width + x + 1,
          sw = (y + 1) * width + x - 1,
          nw = (y - 1) * width + x - 1,
          n2 = (y - 2) * width + x,
          s2 = (y + 2) * width + x;

        if (x + y == origWidth - 1) {  // NW boundary
          if (mask[p] == GREENPOSITION) {
            // IDW-average the closest blue pixels.
            //
            //     B B . .
            //   * e . . .
            // B B . . B b
            // . . . . . .
            // . . b b . .
            //
            oblue[p] = (oblue[sw]/DIAG   + oblue[s]      + oblue[ne]/DIAG + oblue[ne + 1]/DIAG12) / (1/DIAG + 1 + 1/DIAG + 1/DIAG12);
            oblue[e] = (oblue[sw]/DIAG12 + oblue[s]/DIAG + oblue[ne]      + oblue[ne + 1]/DIAG + oblue[se + 2]/DIAG12) / (1/DIAG12 + 1/DIAG + 1 + 1/DIAG + 1/DIAG12);
          }
        }
        else if (y == x - origWidth) {  // NE boundary
          if (mask[p] == GREENPOSITION) {
            // IDW-average the underlying red block to the inner green pixel
            // (w) and copy it to the outer green pixel (*). All other blue
            // pixels are too far to be useful.
            //
            // b b . .
            // . . . w *
            // . . B B . .
            // . . . . . .
            // b b . . b b
            //
            oblue[w] = (oblue[sw - 1]/DIAG + oblue[sw]) / (1 + 1/DIAG);
            oblue[p] = oblue[w];
          }
          else { // REDPOSITION
            // Each red block can be interpolated from the closest two blue blocks.
            //
            //  Inner       Outer
            //
            // . . .       . . .
            // . B w *     . B - *
            // . . | \ .   . . / | .
            // . . B B .   . . B B .
            // . . . . .   . . . . .
            //
            oblue[p] = (oblue[s2]/2 + oblue[s2 - 1]/DIAG12 + oblue[w2]/2) / (1 + 1 / DIAG12);
            oblue[w] = (oblue[s2 - 1]/2 + oblue[s2]/DIAG12 + oblue[w2]) / (1.5 + 1 / DIAG12);
          }
        }
        else if (y == x + origWidth - 1) {  // SW boundary
          if (mask[p] == GREENPOSITION) {
            // IDW-average the closest blue pixels.
            //
            // . . B B . .
            // . . . . . .
            // B B . . B B
            //   * e . . .
            //     B B . .
            //
            oblue[p] = (oblue[nw]/DIAG   + oblue[n]      + oblue[se]/DIAG + oblue[se + 1]/DIAG12) / (1/DIAG + 1 + 1/DIAG + 1/DIAG12);
            oblue[e] = (oblue[nw]/DIAG12 + oblue[n]/DIAG + oblue[se]      + oblue[se + 1]/DIAG + oblue[ne + 2]/DIAG12) / (1/DIAG12 + 1/DIAG + 1 + 1/DIAG + 1/DIAG12);
          }
          // else { // REDPOSITION
          //   // Each red block can be interpolated from the closest two blue blocks.
          //   //
          //   //   Outer        Inner
          //   //
          //   // . . . . .    . . . . .
          //   // . B B . .    . B B . .
          //   // . | / . .    . \ | . .
          //   //   * - B .      * R B .
          //   //     . . .        . . .
          //   //
          //   oblue[p] = (oblue[n2] / 2 + oblue[n2 + 1] / DIAG12 + oblue[e2] / 2) / (1 + 1 / DIAG12);
          //   oblue[e] = (oblue[n2 - 1] / DIAG12 + oblue[n2] / 2 + oblue[e]) / (1.5 + 1 / DIAG12);
          // }
        }
        else if (x + y == origWidth + 2 * origHeight - 2) {  // SE boundary
          if (mask[p] == GREENPOSITION) {
            // IDW-average the overlying blue block to the inner green pixel
            // (w) and copy it to the outer green pixel (*). All other blue
            // pixels are too far to be useful.
            //
            // b b . . b b
            // . . . . . .
            // . . B B . .
            // . . . w *
            // b b . .
            //
            oblue[w] = (oblue[nw] + oblue[nw - 1] / DIAG) / (1 + 1 / DIAG);
            oblue[p] = oblue[w];
          }
          else { // REDPOSITION
            // Each red block can be interpolated from the closest two blue blocks.
            //
            //   Outer        Inner
            //
            // . . . . .    . . . . .
            // . . B B .    . . B B .
            // . . \ | .    . . | / .
            // . B - *      . B w *
            // . . .        . . .
            //
            oblue[p] = (oblue[n2] / 2 + oblue[n2 - 1] / DIAG12 + oblue[w2] / 2) / (1 + 1 / DIAG12);
            oblue[w] = oblue[p];
          }
        }
        else if ( // Green interior and east boundaries
          mask[p] == GREENPOSITION &&
          x + y >= origWidth - 1 + 2 &&  // exclude NW boundary (even though it should have filled right)
          x > y - origWidth + 2 &&       // exclude SW boundary (even though it should have filled right)
          y > x - origWidth - 1 + 2 &&                // exclude NE boundary
          x + y < origWidth + 2 * origHeight - 1 - 2  // exclude SE boundary
        ) {
          // The interpolation pattern for green pixels depends on their horizontal position.
          if ((x + y + 3) % 4 == 0) {
            //
            //   0 1 2 3
            // . . . . . . .
            // B B . . . . .
            // . * . . . . .
            // . . B . . . .
            // . . . . . . .
            //
            oblue[p] = (oblue[nw] / DIAG + oblue[se] / DIAG + oblue[n]) / (1 + 2 / DIAG);
          }
          if ((x + y + 3) % 4 == 1) {
            //
            //   0 1 2 3
            // . . . . . . .
            // . B . . . . .
            // . . * . . . .
            // . . B B . . .
            // . . . . . . .
            //
            oblue[p] = (oblue[nw] / DIAG + oblue[se] / DIAG + oblue[s]) / (1 + 2 / DIAG);
          }
          if ((x + y + 3) % 4 == 2) {
            //
            //   0 1 2 3
            // . . . . . . .
            // . . . . B . .
            // . . . * . . .
            // . . B B . . .
            // . . . . . . .
            //
            oblue[p] = (oblue[ne] / DIAG + oblue[sw] / DIAG + oblue[s]) / (1 + 2 / DIAG);
          }
          if ((x + y + 3) % 4 == 3) {
            //
            //   0 1 2 3
            // . . . . . . .
            // . . . . B B .
            // . . . . * . .
            // . . . B . . .
            // . . . . . . .
            //
            oblue[p] = (oblue[ne] / DIAG + oblue[sw] / DIAG + oblue[n]) / (1 + 2 / DIAG);
          }
        }

        else if ( // Red interior (boundaries have already been excluded above)
          mask[p] == REDPOSITION
        ) {
          if (x % 2 == 0) {
            oblue[p] = (oblue[n2] / 2 + oblue[e2] / 2 + oblue[s2] / 2 + oblue[w]) / 2.5;
          }
          else {
            oblue[p] = (oblue[n2] / 2 + oblue[e] + oblue[s2] / 2 + oblue[w2] / 2) / 2.5;
          }
        }
      }
    }
  }
#ifdef DUMP_STAGES
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int p = y * width + x;
      data_out[p] = 0;
      data_out[p + width * height] = 0;
      data_out[p + 2 * width * height] = oblue[p];
    }
  }
  /* TIFF RGB float->8bit output */
  write_tiff_rgb_f32("blue-diff-interpolated.tiff", data_out, width, height);
#endif


  // Interpolate red differences making the average of possible values depending on the CFA location
  for (int x = 0; x < width; x++) {
    for (int y = 0; y < height; y++) {
      int p = y * width + x;
      if (mask[p] != REDPOSITION && mask[p] != BLANK) {
        int
          n = (y - 1) * width + x,
          s = (y + 1) * width + x,
          e = p + 1,
          w = p - 1,
          e2 = p + 2,
          w2 = p - 2,
          ne = (y - 1) * width + x + 1,
          se = (y + 1) * width + x + 1,
          sw = (y + 1) * width + x - 1,
          nw = (y - 1) * width + x - 1,
          n2 = (y - 2) * width + x,
          s2 = (y + 2) * width + x;


        if (x + y == origWidth - 1) {  // NW boundary
          if (mask[p] == GREENPOSITION) {
            // IDW-average the underlying red block to the inner green pixel
            // (G) and copy it to the outer green pixel (*). All other red
            // pixels are too far to be useful.
            //
            //     . . r r
            //   * e . . .
            // . . R R . .
            // . . . . . .
            // r r . . r r
            //
            ored[e] = (ored[se] + ored[se + 1] / DIAG) / (1 + 1 / DIAG);
            ored[p] = ored[e];
          }
          else { // BLUEPOSITION
            // Each blue block can be interpolated from the closest two red blocks.
            //
            //  Inner       Outer
            //
            //     . . .       . . .
            //   * e R .     * - R .
            // . / | . .   . | \ . .
            // . R R . .   . R R . .
            // . . . . .   . . . . .
            //
            ored[e] = (ored[s2] / DIAG12 + ored[s2 + 1] / 2 + ored[e2]) / (1.5 + 1 / DIAG12);
            ored[p] = (ored[s2] / 2 + ored[s2 + 1] / DIAG12 + ored[e2] / 2) / (1 + 1 / DIAG12);
          }
        }
        else if (x == y - origWidth + 1) { // SW boundary
          if (mask[p] == GREENPOSITION) {
            // IDW-average the overlying red block to the inner green pixel
            // (e) and copy it to the outer green pixel (*). All other red
            // pixels are too far to be useful.
            //
            // r r . . r r
            // . . . . . .
            // . . R R . .
            //   * e . . .
            //     . . r r
            //
            ored[e] = (ored[ne] + ored[ne + 1] / DIAG) / (1 + 1 / DIAG);
            ored[p] = ored[e];
          }
          else { // BLUEPOSITION
            // Each blue block can be interpolated from the closest two red blocks.
            //
            //   Outer        Inner
            //
            // . . . . .    . . . . .
            // . R R . .    . R R . .
            // . | / . .    . \ | . .
            //   * - R .      * e R .
            //     . . .        . . .
            //
            ored[p] = (ored[n2] / 2 + ored[n2 + 1] / DIAG12 + ored[e2] / 2) / (1 + 1 / DIAG12);
            ored[e] = (ored[n2] / DIAG12 + ored[n2 + 1] / 2 + ored[e2]) / (1.5 + 1 / DIAG12);
          }
        }
        else if ( // Green interior and east boundaries
          mask[p] == GREENPOSITION &&
          x + y >= origWidth + 1 &&  // exclude NW boundary (it has been filled)
          x > y - origWidth + 2      // exclude SW boundary (it has been filled)
        ) {
          // The interpolation pattern for green pixels depends on their horizontal position.
          if ((x + y + 1) % 4 == 0) {
            //
            //   0 1 2 3
            // . . . . . . .
            // R R . . . . .
            // . * . . . . .
            // . . R . . . .
            // . . . . . . .
            //
            ored[p] = (ored[nw] / DIAG + ored[se] / DIAG + ored[n]) / (1 + 2 / DIAG);
          }
          if ((x + y + 1) % 4 == 1) {
            //
            //   0 1 2 3
            // . . . . . . .
            // . R . . . . .
            // . . * . . . .
            // . . R R . . .
            // . . . . . . .
            //
            ored[p] = (ored[nw] / DIAG + ored[se] / DIAG + ored[s]) / (1 + 2 / DIAG);
          }
          if ((x + y + 1) % 4 == 2) {
            //
            //   0 1 2 3
            // . . . . . . .
            // . . . . R . .
            // . . . * . . .
            // . . R R . . .
            // . . . . . . .
            //
            ored[p] = (ored[ne] / DIAG + ored[sw] / DIAG + ored[s]) / (1 + 2 / DIAG);
          }
          if ((x + y + 1) % 4 == 3) {
            //
            //   0 1 2 3
            // . . . . . . .
            // . . . . R R .
            // . . . . * . .
            // . . . R . . .
            // . . . . . . .
            //
            ored[p] = (ored[ne] / DIAG + ored[sw] / DIAG + ored[n]) / (1 + 2 / DIAG);
          }
        }

        else if ( // Blue interior (blue does not occur on the west boundaries)
          mask[p] == BLUEPOSITION &
          x + y >= origWidth + 2 &&  // exclude NW boundary (even though it should have filled right)
          x > y - origWidth + 2      // exclude SW boundary (even though it should have filled right)
        ) {
          if (x % 2 == 0) {
            ored[p] = (ored[n2] / 2 + ored[e2] / 2 + ored[s2] / 2 + ored[w]) / 2.5;
          }
          else {
            ored[p] = (ored[n2] / 2 + ored[e] + ored[s2] / 2 + ored[w2] / 2) / 2.5;
          }
        }
      }
    }
  }
#ifdef DUMP_STAGES
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int p = y * width + x;
      data_out[p] = ored[p];
      data_out[p + width * height] = 0;
      data_out[p + 2 * width * height] = 0;
    }
  }
  /* TIFF RGB float->8bit output */
  write_tiff_rgb_f32("red-diff-interpolated.tiff", data_out, width, height);
#endif


  // Make back the differences
  for (int i=0; i < width * height; i++){
    ored[i] += ogreen[i];
    ored[i] *= 1.565476;
    oblue[i] += ogreen[i];
    oblue[i] *= 1.845238;
  }

  free(mask);
}



/**
 * \brief  NLmeans based demosaicking
 *
 * For each value to be filled, a weigthed average of original CFA values of the same channel is performed.
 * The weight depends on the difference of a 3x3 color patch
 *
 * @param[in]  ired, igreen, iblue  initial demosaicked image
 * @param[out] ored, ogreen, oblue  demosaicked output
 * @param[in]  (redx, redy)  coordinates of the red pixel: (0,0), (0,1), (1,0), (1,1)
 * @param[in]  bloc  research block of size (2+bloc+1) x (2*bloc+1)
 * @param[in]  h kernel bandwidth
 * @param[in]  width, height size of the image
 *
 */


void demosaic_nlmeans(int bloc, float h,int redx,int redy,float *ired,float *igreen,float *iblue,float *ored,float *ogreen,float *oblue,int width,int height)
{


  // Initializations
  int bluex = 1 - redx;
  int bluey = 1 - redy;



  // CFA Mask of color per pixel
  unsigned char *cfamask = new unsigned char[width*height];


  for(int x=0;x<width;x++)
    for(int y=0;y<height;y++){

      if (x%2 == redx && y%2 == redy) cfamask[y*width+x] = REDPOSITION;
      else  if (x%2 == bluex && y%2 == bluey) cfamask[y*width+x] = BLUEPOSITION;
      else  cfamask[y*width+x] = GREENPOSITION;

    }


  wxCopy(ired,ored,width*height);
  wxCopy(igreen,ogreen,width*height);
  wxCopy(iblue,oblue,width*height);



  // Tabulate the function Exp(-x) for x>0.
  int luttaille = (int) (LUTMAX*LUTPRECISION);
  float *lut = new float[luttaille];

  sFillLut(lut, luttaille);




  // for each pixel
  for(int y=2; y <height-2; y++)
    for(int x=2; x<width-2; x++)
    {
      // index of current pixel
      int l=y*width+x;


      // Learning zone depending on the window size
      int imin=MAX(x-bloc,1);
      int jmin=MAX(y-bloc,1);

      int imax=MIN(x+bloc,width-2);
      int jmax=MIN(y+bloc,height-2);


      // auxiliary variables for computing average
      float red=0.0;
      float green=0.0;
      float blue=0.0;

      float rweight=0.0;
      float gweight=0.0;
      float bweight=0.0;


      // for each pixel in the neighborhood
      for(int j=jmin;j<=jmax;j++)
        for(int i=imin;i<=imax;i++) {

          // index of neighborhood pixel
          int l0=j*width+i;

          // We only interpolate channels differents of the current pixel channel
          if (cfamask[l]!=cfamask[l0]) {


            // Distances computed on color
            float some = 0.0;

            some = l2_distance_r1(ired,  x, y, i,
                        j, width);
            some += l2_distance_r1(igreen,  x, y, i,
                         j, width);
            some += l2_distance_r1(iblue,  x, y, i,
                         j, width);



            // Compute weight
            some= some / (27.0 * h);

            float weight = sLUT(some,lut);

            // Add pixel to corresponding channel average

            if (cfamask[l0] == GREENPOSITION)  {

              green += weight*igreen[l0];
              gweight+= weight;

            } else if (cfamask[l0] == REDPOSITION) {

              red += weight*ired[l0];
              rweight+= weight;

            } else {

              blue += weight*iblue[l0];
              bweight+= weight;
            }

          }

        }


      // Set value to current pixel
      if (cfamask[l] != GREENPOSITION && gweight > fTiny)  ogreen[l]  =   green / gweight;
      else  ogreen[l] = igreen[l];

      if ( cfamask[l] != REDPOSITION && rweight > fTiny)  ored[l]  =  red / rweight ;
      else    ored[l] = ired[l];

      if  (cfamask[l] != BLUEPOSITION && bweight > fTiny)   oblue[l] =  blue / bweight;
      else  oblue[l] = iblue[l];


    }

  delete[] cfamask;
  delete[] lut;

}




/**
 * \brief  Iterate median filter on chromatic components of the image
 *
 *
 * @param[in]  ired, igreen, iblue  initial  image
 * @param[in]  iter  number of iteracions
 * @param[out] ored, ogreen, oblue  filtered output
 * @param[in]  (redx, redy)  coordinates of the red pixel: (0,0), (0,1), (1,0), (1,1)
 * @param[in]  side  median in a (2*side+1) x (2*side+1) window
 * @param[in]  projflag if not zero, values of the original CFA are kept
 * @param[in]  width, height size of the image
 *
 */


void chromatic_median(int iter,int redx,int redy,int projflag,float side,float *ired,float *igreen, float *iblue,float *ored,float *ogreen,float *oblue,int width,int height)
{


  int size=height*width;

  // Auxiliary variables for computing chromatic components
  float *y=new float[size];
  float *u=new float[size];
  float *v=new float[size];
  float *u0=new float[size];
  float *v0=new float[size];

  int bluex=1-redx;
  int bluey=1-redy;


  // For each iteration
  for(int i=1;i<=iter;i++){


    // Transform to YUV
    wxRgb2Yuv(ired,igreen,iblue,y,u,v,width,height);

    // Perform a Median on YUV component
    wxMedian(u,u0,side,1 ,width,height);
    wxMedian(v,v0,side,1 ,width,height);


    // Transform back to RGB
    wxYuv2Rgb(ored,ogreen,oblue,y,u0,v0,width,height);

    // If projection flag activated put back original CFA values
    if (projflag)
    {

      for(int x=0;x<width;x++)
        for(int y=0;y<height;y++){

          int l=y*width+x;

          if (x%2==redx && y%2==redy) ored[l]=ired[l];

          else if (x%2==bluex && y%2==bluey) oblue[l]=iblue[l];

          else ogreen[l]=igreen[l];

        }

    }


    wxCopy(ored,ired,size);
    wxCopy(ogreen,igreen,size);
    wxCopy(oblue,iblue,size);

  }


  // delete auxiliary memory
  delete[] y;
  delete[] u;
  delete[] u0;
  delete[] v;
  delete[] v0;

}


/**
 * \brief Demosaicking chain
 *
 *
 *
 * Compute initialization by Adams-Hamilton algorithm (u0)
 *
 * for h in {16,4,1} do
 * {
 *
 *    u <- NL_h(u0);      Apply NLmeans demosaicking
 *
 *    u <- CR(u);       Apply chromatic regularization
 *
 *    u0 <- u;
 *
 * }
 *
 * Output <- u;
 *
 *
 * @param[in]  ired, igreen, iblue  initial  image
 * @param[out] ored, ogreen, oblue  filtered output
 * @param[in]  (redx, redy)  coordinates of the red pixel: (0,0), (0,1), (1,0), (1,1)
 * @param[in]  width, height size of the image
 *
 */



void ssd_demosaic_chain(
  float *input,
  float *ored,
  float *ogreen,
  float *oblue,
  int width,
  int height,
  int origWidth,
  int origHeight
) {

  ////////////////////////////////////////////// Process

  float h;
  int dbloc = 7;
  float side = 1.5;
  int iter = 1;
  int projflag = 1;
  float threshold = 2.0;


  adams_hamilton(threshold, input, ored, ogreen, oblue, width, height, origWidth, origHeight);



//  h = 16.0;
//  demosaic_nlmeans(dbloc,h,redx,redy,ored,ogreen,oblue,ired,igreen,iblue,width,height);
//   chromatic_median(iter,redx,redy,projflag,side,ired,igreen,iblue,ored,ogreen,oblue,width,height);
//
//
//
//   h = 4.0;
//   demosaic_nlmeans(dbloc,h,redx,redy,ored,ogreen,oblue,ired,igreen,iblue,width,height);
//   chromatic_median(iter,redx,redy,projflag,side,ired,igreen,iblue,ored,ogreen,oblue,width,height);
//
//
//
//   h = 1.0;
//   demosaic_nlmeans(dbloc,h,redx,redy,ored,ogreen,oblue,ired,igreen,iblue,width,height);
//   chromatic_median(iter,redx,redy,projflag,side,ired,igreen,iblue,ored,ogreen,oblue,width,height);
}


