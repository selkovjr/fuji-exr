#include "cfa_mask.h"


// CFA Mask indicating which color each sensor pixel has

long x, y, p;

unsigned char* cfa_mask(unsigned width, unsigned height, unsigned imageWidth, unsigned imageHeight) {
  unsigned char *mask = new unsigned char[width * height];

  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x++) {
      p = y * width + x;
      if (
          x + y >= imageWidth - 1 and                    // NW boundary
          y > x - imageWidth - 1 and                     // NE boundary
          x + y < imageWidth + 2 * imageHeight - 1 and   // SE boundary
          x > y - imageWidth                             // SW boundary
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

  return mask;
}

