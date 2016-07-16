#include <argp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <error.h>
#include <algorithm>

#include "ssd.h"
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

const char *argp_program_version = VERSION " " __DATE__;

//-------------------------------------------------------------------
// Command-line option parsing.
// Example at https://gist.github.com/sam-github/57f49711cd9073b35d22
//

// ## Top-level parser
static char doc_toplevel[1000];

static error_t parse_toplevel (int key, char *arg, struct argp_state *state) {
  switch (key) {
    case ARGP_KEY_ARG:
      assert( arg );
      if(strcmp(arg, "ssd") == 0) {
        run_ssd(state);
      }
      else if(strcmp(arg, "linear") == 0) {
        // run_linear(state);
      }
      else {
        argp_error(state, "%s is not a valid command", arg);
      }
      break;

    case ARGP_KEY_END:
      if (state->arg_num < 2)
        argp_usage (state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
static struct argp argp_toplevel = {
  0, // no options
  parse_toplevel,
  "<COMMAND> <ARGS>",
  doc_toplevel
};
#pragma GCC diagnostic pop


void interpolate (
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

  // -------------------------------------------------------
  // Do simple linear interpolation for green inside the image
  // -------------------------------------------------------
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
          s = (y + 1) * width + x;

        ogreen[p] = (ogreen[n] + ogreen[s]) / 2.0;
      }
    }
  }
  fprintf(stderr, "green channel interpolated\n");


  // Interpolate blue making the average of possible values in
  // location-dependent patterns (different for G and R locations)
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

        // First interpolate boundary values
        if (x + y == origWidth - 1) {  // NW boundary
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

        // Now interpolate inside the image
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
  fprintf(stderr, "blue interpolated\n");


  // Interpolate red making the average of possible values in
  // location-dependent patterns (different for G and B locations)
  for (int x = 0; x < width; x++) {
    for (int y = 0; y < height; y++) {
      int p = y * width + x;
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

        // Do the interior
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
  fprintf(stderr, "red interpolated\n");

  free(mask);
} // interpolate();



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

  sprintf(
    doc_toplevel,
    "\n"
    "Utilities for processing Fuji EXR sensor data\n"
    "Version: %s\n"
    "\n"
    "Command: linear  interpolate channels without debayering\n"
    "         ssd     self-similarity-driven debayering\n"
    "         db      Duran-Buades debayering\n"
    "\n",
    argp_program_version
  );

  argp_parse (&argp_toplevel, argc, argv, ARGP_IN_ORDER, NULL, NULL);

  exit (0);



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
  //interpolate(ired, igreen, iblue,  ored, ogreen, oblue,  width, height, origWidth, origHeight);
  interpolate (
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
