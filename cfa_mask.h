#ifndef   CFA_MASK_H
#define   CFA_MASK_H

#define BLANK 0
#define REDPOSITION 1
#define GREENPOSITION 2
#define BLUEPOSITION 3

// CFA Mask indicating which color each sensor pixel has
unsigned char* exr_cfa_mask(unsigned width, unsigned height, unsigned imageWidth, unsigned imageHeight);
unsigned char* bggr_cfa_mask(unsigned width, unsigned height);

#endif
