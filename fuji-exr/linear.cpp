#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <error.h>
#include <ctime>
#include <iostream>
#include <iomanip>

#include "linear_args.h"
#include "termcolor.h"
#include "cfa_mask.h"
#include "io_tiff.h"
#include "libAuxiliary.h" // wxCopy()

#define DIAG 1.4142136
#define DIAG12 2.236 // sqrt(5)

using namespace std;
using namespace termcolor;

void interpolate_hires_linear (
  float *ired,
  float *igreen,
  float *iblue,
  float *ored,
  float *ogreen,
  float *oblue,
  int width,
  int height,
  int cfaWidth,
  int cfaHeight,
  unsigned char *mask
);

void interpolate_subframe_linear (
  float *ired,
  float *igreen,
  float *iblue,
  float *ored,
  float *ogreen,
  float *oblue,
  int cfaWidth,
  int cfaHeight,
  unsigned char *mask
);

// --------------------------------------------------------------------
void run_linear (struct argp_state* state) {
  PARSE_ARGS_LINEAR;

  clock_t start_time, end_time;
  double elapsed;

  size_t nx0 = 0, ny0 = 0;
  size_t nx1 = 0, ny1 = 0;
  char *description;
  unsigned long cfaWidth;
  unsigned long cfaHeight;
  unsigned long width, height;
  float *frame0, *frame1;
  float *data_in, *data_out;
  float *out_ptr, *end_ptr;
  bool landscape;
  unsigned long i, x0, x1, y;

  cerr.setf(ios::fixed, ios::floatfield);

  /* TIFF 16-bit grayscale -> float input */
  start_time = clock();
  {
    if (args.hr) {
      cerr << grey << "input file 0: " << white << args.input_file_0 << reset << endl;
      if (NULL == (frame0 = read_tiff_gray16_f32(args.input_file_0, &nx0, &ny0, &description))) {
        cerr << on_red << "error while reading from " << args.input_file_0 << reset << endl;
        exit(EXIT_FAILURE);
      }

      cerr << grey << "input file 1: " << white << args.input_file_1 << reset << endl;
      if (NULL == (frame1 = read_tiff_gray16_f32(args.input_file_1, &nx1, &ny1, &description))) {
        cerr << on_red << "error while reading from " << args.input_file_1 << reset << endl;
        exit(EXIT_FAILURE);
      }
      if (nx0 != nx1 or ny0 != ny1) {
        cerr << on_red << "Input frames must have identical geometry. Got "
          << bold << nx0 << reset << on_red << "×" << bold << ny0 << reset
          << on_red << " and "
          << bold << nx1 << reset << on_red << "×" << bold << ny1 << reset
          << endl;
        exit(EXIT_FAILURE);
      }

      cfaWidth = nx0;
      cfaHeight = ny0;
      width = height = cfaWidth + cfaHeight;
    }
    else {
      cerr << grey << "input file 0: " << white << args.input_file_0 << reset << endl;
      if (NULL == (data_in = read_tiff_gray16_f32(args.input_file_0, &nx0, &ny0, &description))) {
        cerr << on_red << "error while reading from " << args.input_file_0 << reset << endl;
        exit(EXIT_FAILURE);
      }
      if (NULL == (data_out = (float *) malloc(sizeof(float) * nx0 * ny0 * 3))) {
        cerr << on_red << "allocation error: not enough memory" << reset << endl;
        exit(EXIT_FAILURE);
      }

      frame0 = frame1 = data_in; // to avoid -Wmaybe-uninitialized

      width = cfaWidth = nx0;
      height = cfaHeight = ny0;
    }

    landscape = cfaWidth > cfaHeight ? true : false;
  }
  end_time = clock();
  elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
  cerr << yellow << setw(7) << setprecision(2) << elapsed << "s" << white << " reading input" << reset << endl;

  if (args.hr) {
    start_time = clock();
    {
      if (NULL == (data_in = (float *) malloc(sizeof(float) * width * height * 3))) {
        cerr << on_red << "allocation error: not enough memory" << reset << endl;
        exit(EXIT_FAILURE);
      }
      if (NULL == (data_out = (float *) malloc(sizeof(float) * width * height * 3))) {
        cerr << on_red << "allocation error: not enough memory" << reset << endl;
        exit(EXIT_FAILURE);
      }

      for (i = 0; i < width * height * 3; i++) {
        data_in[i] = 0;
      }
    }
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    cerr << yellow << setw(7) << setprecision(2) << elapsed << "s" << white << " allocating and zero-setting memory" << reset << endl;

    start_time = clock();
    for (i = 0; i < (unsigned long)cfaWidth * cfaHeight; i++) {
      if (landscape) {
        //
        // B........G
        // ..........
        // ..........
        // G........R
        //
        x0 = i % cfaWidth + (unsigned long)(i / cfaWidth);
        x1 = x0 + 1; // the second frame (fn == 1) is shifted 1px to the right
        y = (cfaWidth - i % cfaWidth - 1) + (i / cfaWidth);
        data_in[y * width + x0] = frame0[i];
        data_in[y * width + x1] = frame1[i];
        data_in[y * width + x0 + width * height] = frame0[i];
        data_in[y * width + x1 + width * height] = frame1[i];
        data_in[y * width + x0 + width * height * 2] = frame0[i];
        data_in[y * width + x1 + width * height * 2] = frame1[i];
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
        x0 = cfaHeight - 1 + i % cfaWidth - (unsigned long)(i / cfaWidth);
        x1 = x0 + 1; // the second frame (fn == 1) is shifted 1px to the right
        y = i % cfaWidth + (unsigned long)(i / cfaWidth);
        data_in[y * width + x0] = frame0[i];
        data_in[y * width + x1] = frame1[i];
        data_in[y * width + x0 + width * height] = frame0[i];
        data_in[y * width + x1 + width * height] = frame1[i];
        data_in[y * width + x0 + width * height * 2] = frame0[i];
        data_in[y * width + x1 + width * height * 2] = frame1[i];
      }
    }
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    cerr << yellow << setw(7) << setprecision(2) << elapsed << "s" << white << " merging input frames" << reset << endl;

    // write_tiff_rgb_f32("input-merged.tif", data_in, width, height);

    start_time = clock();
    unsigned char *mask = exr_cfa_mask(width, height, cfaWidth, cfaHeight);
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    cerr << yellow << setw(7) << setprecision(2) << elapsed << "s" << white << " computing the CFA mask" << reset << endl;

    start_time = clock();
    interpolate_hires_linear (
      data_in,
      data_in + width * height,
      data_in + 2 * width * height,
      data_out,
      data_out + width * height,
      data_out + 2 * width * height,
      (int) width,
      (int) height,
      landscape ? cfaWidth : cfaHeight,
      landscape ? cfaHeight : cfaWidth,
      mask
    );

    delete[] mask;
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    cerr << yellow << setw(7) << setprecision(2) << elapsed << "s" << white << " interpolating" << reset << endl;
  } // HR

  else { // interpolate one BGGR subframe
    start_time = clock();
    unsigned char *mask = bggr_cfa_mask(width, height);
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    cerr << yellow << setw(7) << setprecision(2) << elapsed << "s" << white << " computing the CFA mask" << reset << endl;

    // It's bitching about data_in and data_out
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
    start_time = clock();
    interpolate_subframe_linear (
      data_in,
      data_in + width * height,
      data_in + 2 * width * height,
      data_out,
      data_out + width * height,
      data_out + 2 * width * height,
      (int)width,
      (int)height,
      mask
    );
#pragma GCC diagnostic pop

    delete[] mask;
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    cerr << yellow << setw(7) << setprecision(2) << elapsed << "s" << white << " interpolating" << reset << endl;
  }

  /* limit to 0-65535 */
  start_time = clock();
  {
    out_ptr = data_out;
    end_ptr = out_ptr + 3 * width * height;
    while (out_ptr < end_ptr) {
      if ( 0 > *out_ptr)
        *out_ptr = 0;
      if ( 65535 < *out_ptr)
        *out_ptr = 65535;
      out_ptr++;
    }
  }
  end_time = clock();
  elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
  cerr << yellow << setw(7) << setprecision(2) << elapsed << "s" << white << " clipping values" << reset << endl;

  cerr << grey << "writing output to " << white << args.output_file << reset << endl;
  start_time = clock();
  {
    write_tiff_rgb_f32(args.output_file, data_out, width, height);
  }
  end_time = clock();
  elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
  cerr << yellow << setw(7) << setprecision(2) << elapsed << "s" << white << " writing" << reset << endl;

  free(data_in);
  free(data_out);

  exit(EXIT_SUCCESS);

} // run_linear()


void interpolate_hires_linear (
  float *ired,
  float *igreen,
  float *iblue,
  float *ored,
  float *ogreen,
  float *oblue,
  int width,
  int height,
  int cfaWidth,
  int cfaHeight,
  unsigned char *mask
) {
  long x, y, p;

  wxCopy(ired, ored, width * height);
  wxCopy(igreen, ogreen, width * height);
  wxCopy(iblue, oblue, width * height);

  // Interpolate the green channel in the 4-pixel-wide boundary by inverse
  // distance weighting.
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      p = y * width + x;
      if (
        (
         mask[p] == BLUEPOSITION ||
         mask[p] == REDPOSITION
        )
        &&
        (
         x + y < cfaWidth + 3 ||                    // NW boundary
         x >= y + cfaWidth - 3 ||                   // NE boundary
         x + y >= cfaWidth + 2 * cfaHeight - 5 ||  // SE boundary
         y >= x + cfaWidth - 4                      // SW boundary
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

  // -------------------------------------------------------
  // Do simple linear interpolation for green inside the image
  // -------------------------------------------------------
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x++) {
      p = y * width + x;
      if (
          mask[p] != GREENPOSITION &&
          x + y >= cfaWidth - 1 + 4 &&                   // NW boundary
          y > x - cfaWidth - 1 + 4 &&                    // NE boundary
          x + y < cfaWidth + 2 * cfaHeight - 1 - 4 &&   // SE boundary
          x > y - cfaWidth + 4                           // SW boundary
         ) {
        int
          n = (y - 1) * width + x,
          s = (y + 1) * width + x;

        ogreen[p] = (ogreen[n] + ogreen[s]) / 2.0;
      }
    }
  }
  cerr << green << "green " << grey << "channel interpolated" << endl;

  // Interpolate blue making the average of possible values in
  // location-dependent patterns (different for G and R locations)
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      p = y * width + x;
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

        // First interpolate boundary values
        if (x + y == cfaWidth - 1) {  // NW boundary
          if (mask[p] == GREENPOSITION) {
            // IDW-average of the closest blue pixels.
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
        else if (y == x - cfaWidth) {  // NE boundary
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
        else if (y == x + cfaWidth - 1) {  // SW boundary
          if (mask[p] == GREENPOSITION) {
            // IDW-average of the closest blue pixels.
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
        }
        else if (x + y == cfaWidth + 2 * cfaHeight - 2) {  // SE boundary
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

        // Now interpolate inside the image
        else if ( // Green interior and east boundaries
          mask[p] == GREENPOSITION &&
          x + y >= cfaWidth - 1 + 2 &&  // exclude NW boundary (even though it should have filled right)
          x > y - cfaWidth + 2 &&       // exclude SW boundary (even though it should have filled right)
          y > x - cfaWidth - 1 + 2 &&                // exclude NE boundary
          x + y < cfaWidth + 2 * cfaHeight - 1 - 2  // exclude SE boundary
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
  cerr << blue << " blue " << grey << "channel interpolated" << endl;


  // Interpolate red making the average of possible values in
  // location-dependent patterns (different for G and B locations)
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      p = y * width + x;
      if (mask[p] != REDPOSITION && mask[p] != BLANK) {
        int
          // n = (y - 1) * width + x,
          // s = (y + 1) * width + x,
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


        // Do the boundaries
        if (x + y == cfaWidth - 1) {  // NW boundary
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
        else if (x == y - cfaWidth + 1) { // SW boundary
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

        // Do the interior
        else if ( // Green interior and east boundaries
          mask[p] == GREENPOSITION &&
          x + y >= cfaWidth + 1 &&  // exclude NW boundary (it has been filled)
          x > y - cfaWidth + 2      // exclude SW boundary (it has been filled)
        ) {
          // The interpolation pattern for green pixels depends on their horizontal position.
          if ((x + y + 1) % 4 == 0) {
            //
            //   0 1 2 3
            // . . . . . . .
            // R . . . . . .
            // . * . . . . .
            // . . R . . . .
            // . . . . . . .
            //
            ored[p] = (ored[nw] + ored[se]) / 2;
          }
          if ((x + y + 1) % 4 == 1) {
            //
            //   0 1 2 3
            // . . . . . . .
            // . R . . . . .
            // . . * . . . .
            // . . . R . . .
            // . . . . . . .
            //
            ored[p] = (ored[nw] + ored[se]) / 2;
          }
          if ((x + y + 1) % 4 == 2) {
            //
            //   0 1 2 3
            // . . . . . . .
            // . . . . R . .
            // . . . * . . .
            // . . R . . . .
            // . . . . . . .
            //
            ored[p] = (ored[ne] + ored[sw]) / 2;
          }
          if ((x + y + 1) % 4 == 3) {
            //
            //   0 1 2 3
            // . . . . . . .
            // . . . . . R .
            // . . . . * . .
            // . . . R . . .
            // . . . . . . .
            //
            ored[p] = (ored[ne] + ored[sw]) / 2;
          }
        }

        else if ( // Blue interior (blue does not occur on the west boundaries)
          mask[p] == BLUEPOSITION &&
          x + y >= cfaWidth + 2 and  // exclude NW boundary (even though it should have filled right)
          x > y - cfaWidth + 2       // exclude SW boundary (even though it should have filled right)
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
  cerr << red << "  red " << grey << "channel interpolated" << endl;

} // interpolate_hires_linear();


void interpolate_subframe_linear (
  float *ired,
  float *igreen,
  float *iblue,
  float *ored,
  float *ogreen,
  float *oblue,
  int width,
  int height,
  unsigned char *mask
) {
  long x, y, p;

  wxCopy(ired, ored, width * height);
  wxCopy(igreen, ogreen, width * height);
  wxCopy(iblue, oblue, width * height);

  // -------------------------------------------------
  // Do simple linear interpolation in the green plane
  // -------------------------------------------------
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x++) {
      p = y * width + x;
      if (mask[p] != GREENPOSITION) {
        int
          n = (y - 1) * width + x,
          s = (y + 1) * width + x,
          e = p + 1,
          w = p - 1;

        if (y == 0) { // N edge
          if (x == 0) { // NW corner
            ogreen[p] = (ogreen[e] + ogreen[s]) / 2.0;
          }
          else {
            ogreen[p] = (ogreen[e] + ogreen[w]) / 2.0;
          }
        }
        else if (x == 0 and y > 0) { // W edge
          ogreen[p] = (ogreen[n] + ogreen[s]) / 2.0;
        }
        else { // interior
          ogreen[p] = (ogreen[n] + ogreen[s] + ogreen[w] + ogreen[e]) / 4.0;
        }
      }
    }
  }
  cerr << green << "green " << grey << "channel interpolated" << endl;

  // Interpolate blue making the average of possible values in
  // location-dependent patterns (different for G and R locations)
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      p = y * width + x;
      if (mask[p] == GREENPOSITION) {
        int
          n = (y - 1) * width + x,
          s = (y + 1) * width + x,
          e = p + 1,
          w = p - 1;

        if (x == width - 1) {  // E edge
          oblue[p] = oblue[w];
        }
        else if (y == height - 1) { // S edge
          oblue[p] = oblue[n];
        }
        else if (x % 2) { // odd interior column
          oblue[p] = (oblue[e] + oblue[w]) / 2.0;
        }
        else { // even interior column
          oblue[p] = (oblue[n] + oblue[s]) / 2.0;
        }
      }

      if (mask[p] == REDPOSITION) {
        int
          ne = (y - 1) * width + x + 1,
          se = (y + 1) * width + x + 1,
          sw = (y + 1) * width + x - 1,
          nw = (y - 1) * width + x - 1;

        if (x == width - 1 or y == height - 1) {  // E and S edges
          oblue[p] = oblue[nw];
        }
        else { // interior
          oblue[p] = (oblue[nw] + oblue[ne] + oblue[se] + oblue[sw]) / 4.0;
        }
      }
    }
  }
  cerr << blue << " blue " << grey << "channel interpolated" << endl;

  // Interpolate red making the average of possible values in
  // location-dependent patterns (different for G and B locations)
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      p = y * width + x;
      if (mask[p] == GREENPOSITION) {
        int
          n = (y - 1) * width + x,
          s = (y + 1) * width + x,
          e = p + 1,
          w = p - 1;

        if (x == 0) {  // W edge
          ored[p] = ored[e];
        }
        else if (y == 0) { // N edge
          ored[p] = ored[s];
        }
        else if (x % 2) { // odd interior column
          ored[p] = (ored[e] + ored[w]) / 2.0;
        }
        else { // even interior column
          ored[p] = (ored[n] + ored[s]) / 2.0;
        }
      }

      if (mask[p] == BLUEPOSITION) {
        int
          ne = (y - 1) * width + x + 1,
          se = (y + 1) * width + x + 1,
          sw = (y + 1) * width + x - 1,
          nw = (y - 1) * width + x - 1;

        if (x == 0 or y == 0) {  // W and N edges
          ored[p] = ored[se];
        }
        else { // interior
          ored[p] = (ored[nw] + ored[ne] + ored[se] + ored[sw]) / 4.0;
        }
      }
    }
  }
  cerr << red << "  red " << grey << "channel interpolated" << endl;

} // interpolate_subframe_linear();
