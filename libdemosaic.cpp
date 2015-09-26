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


#include <algorithm>
#include "libdemosaic.h"

#include "io_tiff.h"
#include "tiffio.h"


#define BLANK 0
#define REDPOSITION 1
#define GREENPOSITION 2
#define BLUEPOSITION 3

#define DIAG 1.4142136
#define DIAG12 2.236 // sqrt(5)

#define DUMP_STAGES
#undef DEBUG_GREEN

/**
 * @file   libdemosaic.cpp
 * @brief  Demosaicking functions: Hamilton-Adams algorithm, NLmeans-based demosaicking, Chromatic components filtering
 *
 * @author Antoni Buades <toni.buades@uib.es>
 * @author Gene Selkov <selkovjr@gmail.com> (Fuji EXR)
 */


/**
 * \brief  Classical Adams-Hamilton demosaicking algorithm (adapted for Fuji EXR)
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

void g_directional(
  float threshold,
  float *ired,
  float *igreen,
  float *iblue,
  float *ored,
  float *ogreen,
  float *oblue,
  int width,
  int height,
  int origWidth,
  int origHeight
) {
  wxCopy(ired, ored, width * height);
  wxCopy(igreen, ogreen, width * height);
  wxCopy(iblue, oblue, width * height);

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
            printf("\x1b[34m%6.0f%6.0f\x1b[31m%6.0f%6.0f\x1b[34m%6.0f\x1b[0m\n",   ired[nw2], ired[nw2 + 1], ired[n2], ired[ne2 - 1], ired[ne2]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    ired[nw],      ired[n],  ired[ne]);
            printf( "      \x1b[31m%6.0f\x1b[1m\x1b[34m%6.0f\x1b[0m\x1b[34m%6.0f\x1b[0m\n",                                    ired[w],       ired[p],  ired[e]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    ired[sw],      ired[s],  ired[se]);
            printf("\x1b[34m%6.0f%6.0f\x1b[31m%6.0f%6.0f\x1b[34m%6.0f\x1b[0m\n",   ired[sw2], ired[sw2 + 1], ired[s2], ired[se2 - 1], ired[se2]);
          }
          else {
            printf("\x1b[34m%6.0f\x1b[31m%6.0f%6.0f\x1b[34m%6.0f%6.0f\x1b[0m\n",   ired[nw2], ired[nw2 + 1], ired[n2], ired[ne2 - 1], ired[ne2]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    ired[nw],      ired[n],  ired[ne]);
            printf( "      \x1b[34m%6.0f\x1b[1m%6.0f\x1b[0m\x1b[31m%6.0f\x1b[0m\n",                                    ired[w],       ired[p],  ired[e]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    ired[sw],      ired[s],  ired[se]);
            printf("\x1b[34m%6.0f\x1b[31m%6.0f%6.0f\x1b[34m%6.0f%6.0f\x1b[0m\n",   ired[sw2], ired[sw2 + 1], ired[s2], ired[se2 - 1], ired[se2]);
          }
        }
        else { // red
          printf("\x1b[31m--------------------------------------------------------------\x1b[0m\n");
          printf("\x1b[31mred\x1b[0m, %s (%d, %d)\n", x % 2 ? "odd" : "even", x, y);
          if (x % 2 == 0) {
            printf("\x1b[31m%6.0f%6.0f\x1b[34m%6.0f%6.0f\x1b[31m%6.0f\x1b[0m\n",   ired[nw2], ired[nw2 + 1], ired[n2], ired[ne2 - 1], ired[ne2]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    ired[nw],      ired[n],  ired[ne]);
            printf( "      \x1b[34m%6.0f\x1b[1m\x1b[31m%6.0f\x1b[0m\x1b[31m%6.0f\x1b[0m\n",                                    ired[w],       ired[p],  ired[e]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    ired[sw],      ired[s],  ired[se]);
            printf("\x1b[31m%6.0f%6.0f\x1b[34m%6.0f%6.0f\x1b[31m%6.0f\x1b[0m\n",   ired[sw2], ired[sw2 + 1], ired[s2], ired[se2 - 1], ired[se2]);
          }
          else {
            printf("\x1b[31m%6.0f\x1b[34m%6.0f%6.0f\x1b[31m%6.0f%6.0f\x1b[0m\n",   ired[nw2], ired[nw2 + 1], ired[n2], ired[ne2 - 1], ired[ne2]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    ired[nw],      ired[n],  ired[ne]);
            printf( "      \x1b[31m%6.0f\x1b[1m%6.0f\x1b[0m\x1b[34m%6.0f\x1b[0m\n",                                    ired[w],       ired[p],  ired[e]);
            printf( "      \x1b[32m%6.0f%6.0f%6.0f\x1b[0m\n",                                    ired[sw],      ired[s],  ired[se]);
            printf("\x1b[31m%6.0f\x1b[34m%6.0f%6.0f\x1b[31m%6.0f%6.0f\x1b[0m\n",   ired[sw2], ired[sw2 + 1], ired[s2], ired[se2 - 1], ired[se2]);
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
  fprintf(stderr, "green channel interpolated\n");

  // compute the bilinear on the differences of the red and blue with the already interpolated green
  bilinear_red_blue(ored, ogreen, oblue, width, height, origWidth, origHeight);

  free(mask);
  fprintf(stderr, "demosaic complete\n");
} // g_directional();




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
        mask[y * width + x] = BLANK;
      }
    }
  }

  // Compute the differences
  for (int i = 0; i < width * height; i++) {
    ored[i] -= ogreen[i];
    oblue[i] -= ogreen[i];
  }

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
  fprintf(stderr, "R-G interpolated\n");


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
          mask[p] == BLUEPOSITION &&
          x + y >= origWidth + 2 and  // exclude NW boundary (even though it should have filled right)
          x > y - origWidth + 2       // exclude SW boundary (even though it should have filled right)
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
  fprintf(stderr, "B-G interpolated\n");

  // Make back the differences
  for (int i = 0; i < width * height; i++){
    ored[i] += ogreen[i];
    // ored[i] *= 1.565476;
    oblue[i] += ogreen[i];
    // oblue[i] *= 1.845238;
  }

  free(mask);
} // bilinear_red_blue()



/**
 * \brief  NLmeans-based demosaicking
 *
 * For each value to be filled, a weigthed average of original CFA values of the same channel is performed.
 * The weight depends on the difference of a 3x3 color patch
 *
 *  NLM(u)(x) = Σ|y ∈ Ω| w(x, y) * u(y), where:
 *
 *  u is a color channel
 *  Ω is the domain of u (the entire color plane)
 *  w(x, y) = 1/wx exp(-L2((Patch(x) - Patch(y)) / 2 * h^2))
 *  wx = Σ|y ∈ Ω| exp(-L2((Patch(x) - Patch(y)) / 2 * h^2))
 *
 *  High similarty between patches Patch(x) and Patch(y) is reflected in
 *  a higher weight w(x, y).
 *
 * @param[in]  ired, igreen, iblue  initial demosaicked image
 * @param[out] ored, ogreen, oblue  demosaicked output
 * @param[in]  bloc search block of size (2+bloc+1) x (2*bloc+1)
 * @param[in]  h kernel bandwidth
 * @param[in]  width, height size of the image
 *
 */


void demosaic_nlmeans(
  int bloc,
  float h,
  float *ired,
  float *igreen,
  float *iblue,
  float *ored,
  float *ogreen,
  float *oblue,
  int width,
  int height,
  int origWidth,
  int origHeight
) {
  fprintf(stderr, "running NLM interpolation...\n");

  // CFA Mask indicating which color each sensor pixel has
  unsigned char *cfamask = new unsigned char[width*height];

  wxCopy(ired, ored, width * height);
  wxCopy(igreen, ogreen, width * height);
  wxCopy(iblue, oblue, width * height);


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
          cfamask[p] = GREENPOSITION;
        }
        else {
          if ((x + y - 1) % 4 == 0 || (x + y - 1) % 4 == 1) {
            cfamask[p] = REDPOSITION;
          }
          else {
            cfamask[p] = BLUEPOSITION;
          }
        }
      }
      else {
        cfamask[p] = BLANK;
      }
    }
  }


  // Tabulate the function Exp(-x) for x > 0.
  int luttaille = (int) (LUTMAX * LUTPRECISION);
  float *lut = new float[luttaille];
  sFillLut(lut, luttaille);


  // for each pixel in the interior
  for (int y = 2; y < height - 2; y++) {
    for (int x = 2; x < width - 2; x++) {
      int p = y * width + x;
      if (
        cfamask[p] != BLANK and
        (
         x + y >= origWidth + 3 + bloc - 1 and                    // NW boundary
         x < y + origWidth - 3 - bloc + 1 and                   // NE boundary
         x + y < origWidth + 2 * origHeight - 5 - bloc + 1 and  // SE boundary
         y < x + origWidth - 4 - bloc + 1                      // SW boundary
        )
      ) {
        // Learning zone depending on window size
        int imin = MAX(x - bloc, 1);
        int jmin = MAX(y - bloc, 1);

        int imax = MIN(x + bloc, width - 2);
        int jmax = MIN(y + bloc, height - 2);

        // auxiliary variables for computing average
        float red = 0.0;
        float green = 0.0;
        float blue = 0.0;

        float rweight = 0.0;
        float gweight = 0.0;
        float bweight = 0.0;


        // for each pixel in the neighborhood
        for (int j = jmin; j <= jmax; j++) {
          for (int i = imin; i <= imax; i++) {

            // index of neighborhood pixel
            int n = j * width + i;

            // We only interpolate channels other than the current pixel channel
            if (cfamask[p] != cfamask[n]) {

              // Distances computed on color
              float sum = 0.0;

              sum = l2_distance_r1(ired,  x, y, i, j, width);
              sum += l2_distance_r1(igreen,  x, y, i, j, width);
              sum += l2_distance_r1(iblue,  x, y, i, j, width);

              // Compute weight
              sum /= (65536 * 27.0 * h); // The original was probably tuned to 8-bit images (so the sum is 256^2 larger)
              // sum /= (8192 * 27.0 * h); // this seems to produce a more agreeable denoising on red

              // weight = exp(-sum)
              float weight = sLUT(sum, lut);

              // Add pixel to corresponding channel average
              if (cfamask[n] == GREENPOSITION)  {
                green += weight * igreen[n];
                gweight += weight;
              }
              else if (cfamask[n] == REDPOSITION) {
                red += weight * ired[n];
                rweight += weight;
              }
              else {
                blue += weight * iblue[n];
                bweight += weight;
              }

            }

          }
        }


        // Set value to current pixel
        if (cfamask[p] != GREENPOSITION && gweight > fTiny) ogreen[p] = green / gweight;
        else ogreen[p] = igreen[p];

        if ( cfamask[p] != REDPOSITION && rweight > fTiny) ored[p] = red / rweight;
        else ored[p] = ired[p];

        if (cfamask[p] != BLUEPOSITION && bweight > fTiny) oblue[p] = blue / bweight;
        else  oblue[p] = iblue[p];
      }
    }
  }

  delete[] cfamask;
  delete[] lut;
  fprintf(stderr, "NLM denoising done for h = %f\n", h);
}


/**
 * \brief  Iterate median filter on chromatic components of the image
 *
 *
 * @param[in]  ired, igreen, iblue  initial  image
 * @param[in]  iter  number of iteracions
 * @param[out] ored, ogreen, oblue  filtered output
 * @param[in]  side  median in a (2*side+1) x (2*side+1) window
 * @param[in]  projflag if not zero, values of the original CFA are kept
 * @param[in]  width, height size of the image
 *
 */


void chromatic_median(
  int iter,
  int projflag,
  float side,
  float *ired,
  float *igreen,
  float *iblue,
  float *ored,
  float *ogreen,
  float *oblue,
  int width,
  int height,
  int origWidth,
  int origHeight
) {
  fprintf(stderr, "%d iterations of chromatic median\n", iter);

  int size = height * width;

  // Auxiliary variables for computing chromatic components
  float *Y = new float[size];
  float *U = new float[size];
  float *V = new float[size];
  float *U0 = new float[size];
  float *V0 = new float[size];

  // For each iteration
  for (int i = 1; i <= iter; i++) {
    // Transform to YUV
    wxRgb2Yuv(ired, igreen, iblue, Y, U, V, width, height, origWidth, origHeight);

    // Perform a Median on YUV component.
    // The filtered image U0 is copied back to U inside. So is V0 -> V.
    wxMedian(U, U0, side, 1, width, height, origWidth, origHeight);
    wxMedian(V, V0, side, 1, width, height, origWidth, origHeight);

    // Transform back to RGB
    wxYuv2Rgb(ored, ogreen, oblue, Y, U0, V0, width, height);

    // If projection flag is set, put back original CFA values
    if (projflag) {
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
              ogreen[p] = igreen[p];
            }
            else {
              if ((x + y - 1) % 4 == 0 || (x + y - 1) % 4 == 1) {
                // cfamask[p] = REDPOSITION;
                ored[p] = ired[p];
              }
              else {
                oblue[p] = iblue[p];
              }
            }
          }
        }
      }
    }

    wxCopy(ored, ired, size);
    wxCopy(ogreen, igreen, size);
    wxCopy(oblue, iblue, size);
  }

  // delete auxiliary memory
  delete[] Y;
  delete[] U;
  delete[] U0;
  delete[] V;
  delete[] V0;
  fprintf(stderr, "chromatic median done\n");
}


/** \brief Demosaicking chain
 *
 *
 *
 * Compute initialization by Adams-Hamilton algorithm (u0)
 *
 * for h in {16,4,1} do {
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
 * Sketch of the SSD algorithm:
 *
 *  * For u the red, green or blue channel, let Ωu be the domain of u.
 *
 *  * For x ∉ Ωu, SID-NLM(u)(x) = Σ|y ∈ Ωu| w(x, y) * u(y), where w is computed
 *  on an initial color estimate obtained by a standard demosaicking algorithm.
 *
 *  * To reduce artefacts caused by the initial estimate, a coarse-to-fine
 *  strategy is used by iteratively applying SSD-NLM with a decreasing
 *  sequence of h together with a color regularization step (median filtering
 *  on chromaticity components).
 *
 *
 * @param[in]  ired, igreen, iblue  initial  image
 * @param[out] ored, ogreen, oblue  filtered output
 *
 */


void ssd_demosaic_chain(
  float *ired,
  float *igreen,
  float *iblue,
  float *ored,
  float *ogreen,
  float *oblue,
  int width,
  int height,
  int origWidth,
  int origHeight
) {

  ////////////////////////////////////////////// Process

  int dbloc = 7;
  float side = 1.5;
  int iter = 2;
  int projflag = 1;
  float threshold = 2.0;

  // g_directional(threshold,     ired, igreen, iblue,  ored, ogreen, oblue,  width, height, origWidth, origHeight);
  // write_image("demosaicked.tiff",                    ored, ogreen, oblue,  width, height);
  // //                                          ______/      /      /
  // //                                        /      _______/      /
  // //                                       /      /      _______/
  // //                                      /      /      /
  // chromatic_median(iter, projflag, side,  ored, ogreen, oblue,  ired, igreen, iblue,  width, height, origWidth, origHeight);
  // //                                                              \    |    /
  // //                                                               \   |   /
  // //                                                                 output
  // write_image("median.tiff",                                    ired, igreen, iblue,  width, height);

  g_directional(threshold,     ired, igreen, iblue,  ored, ogreen, oblue,  width, height, origWidth, origHeight);
  write_image("demosaicked.tiff",                    ored, ogreen, oblue,  width, height);
  //                                  ________________/      /      /
  //                                /      _________________/      /
  //                               /      /      _________________/
  //                              /      /      /
  demosaic_nlmeans(dbloc, 16,  ored, ogreen, oblue,  ired, igreen, iblue,  width, height, origWidth, origHeight);
  write_image("nlmeans-16.tiff",                     ired, igreen, iblue,  width, height);
  //                                            ______/      /      /
  //                                           /      ______/      /
  //                                          /      /      ______/
  //                                         /      /      /
  chromatic_median(iter, projflag, side,  ired, igreen, iblue,  ored, ogreen, oblue,  width, height, origWidth, origHeight);
  write_image("median-16.tiff",                                 ored, ogreen, oblue,  width, height);

  // demosaic_nlmeans(dbloc, 4,  ored, ogreen, oblue,  ired, igreen, iblue,  width, height, origWidth, origHeight);
  // write_image("nlmeans-4.tiff",                     ired, igreen, iblue,  width, height);
  // chromatic_median(iter, projflag, side,  ired, igreen, iblue,  ored, ogreen, oblue,  width, height, origWidth, origHeight);
  // write_image("median-4.tiff",                                  ored, ogreen, oblue,  width, height);

  // demosaic_nlmeans(dbloc, 1,  ored, ogreen, oblue,  ired, igreen, iblue,  width, height, origWidth, origHeight);
  // write_image("nlmeans-1.tiff",                     ired, igreen, iblue,  width, height);
  // chromatic_median(iter, projflag, side,  ired, igreen, iblue,  ored, ogreen, oblue,  width, height, origWidth, origHeight);
  // write_image("median-1.tiff",                                  ored, ogreen, oblue,  width, height);
}


