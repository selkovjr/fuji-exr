/*
 * Copyright 2009, 2010 IPOL Image Processing On Line
 *  <http://www.ipol.im/>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <string.h>
#include "io_tiff.h"


/*
 * READ FUNCTIONS
 */

/**
 * @brief load the data from a TIFF image file as a float array
 *
 * The array is allocated by this function.
 *
 * @param fname the file name to read
 * @param nx, ny storage space for the image size
 *
 * @return the data array pointer, NULL if an error occured
 */
float *read_tiff_gray16_f32(const char *fname, size_t *nx, size_t *ny, char **description)
{
  TIFF *fp = NULL;
  uint32 width = 0,
         height = 0;
  char *c = NULL;
  uint32 config;
  uint32 row;
  uint32 i;
  tdata_t buf = NULL;
  uint16 s, nsamples;
  float *data = NULL;
  float *ptr_r, *ptr_g, *ptr_b;

  /* no warning messages */
  (void) TIFFSetWarningHandler(NULL);

  /* open the TIFF file and structure */
  if (NULL == (fp = TIFFOpen(fname, "r")))
    return NULL;

  /* read width and height and allocate the storage raster */
  if (
    1 != TIFFGetField(fp, TIFFTAG_IMAGEWIDTH, &width) ||
    1 != TIFFGetField(fp, TIFFTAG_IMAGELENGTH, &height) ||
    1 != TIFFGetField(fp, TIFFTAG_IMAGEDESCRIPTION, &c) ||
    1 != TIFFGetField(fp, TIFFTAG_PLANARCONFIG, &config) ||
    NULL == (buf = _TIFFmalloc(TIFFScanlineSize(fp))) ||
    NULL == (data = (float *) malloc(3 * width * height * sizeof(float)))
  ) {
    TIFFClose(fp);
    return NULL;
  }

  if (NULL != nx)
    *nx = (size_t) width;
  if (NULL != ny)
    *ny = (size_t) height;
  if (NULL != description)
    *description = strdup(c);

  TIFFGetField(fp, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
  printf("samples: %d, width: %d, height: %d\n", nsamples, width, height);

  /* setup the pointers */
  ptr_r = data;
  ptr_g = ptr_r + width * height;
  ptr_b = ptr_g + width * height;

  /*
   * Deinterlace the TIFF raster into three arrays (ptr_r, ptr_g, ptr_b)
   */

  for (s = 0; s < nsamples; s++) {
    for (row = 0; row < height; row++) {
      TIFFReadScanline(fp, buf, row, s);
      for (i = 0; i < width; i++) {
        *ptr_r++ = (float) ((uint16 *)buf)[i];
        *ptr_g++ = (float) ((uint16 *)buf)[i];
        *ptr_b++ = (float) ((uint16 *)buf)[i];
      }
    }
  }


  _TIFFfree(buf);
  TIFFClose(fp);

  return data;
}

/*
 * WRITE FUNCTIONS
 */

/**
 * @brief save an array into a LZW-compressed TIFF file
 *
 * @return 0 if OK, != 0 if an error occured
 */
static int write_tiff_rgb_raw(const char *fname, const uint16 *data_tiff, size_t nx, size_t ny) {
  TIFF *fp = NULL;

  printf("write_tiff_rgb_raw() -> %d\n", data_tiff[8990 * 3]);
  /*
   * ensure the data is allocated
   * and the width and height are within the limits
   * (tiff uses uint32, 2^32 - 1 = 4294967295)
   * and open the TIFF file and structure
   */
  if (NULL == data_tiff || 4294967295. < (double) nx || 4294967295. < (double) ny)
    return -1;

  /* open the TIFF file and structure */
  if (NULL == (fp = TIFFOpen(fname, "w")))
    return -1;

  /* insert tags into the TIFF structure */
  if (1 != TIFFSetField(fp, TIFFTAG_IMAGEWIDTH, nx)
    || 1 != TIFFSetField(fp, TIFFTAG_IMAGELENGTH, ny)
    || 1 != TIFFSetField(fp, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT)
    || 1 != TIFFSetField(fp, TIFFTAG_BITSPERSAMPLE, 16)
    || 1 != TIFFSetField(fp, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG)
    || 1 != TIFFSetField(fp, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(fp, nx * 3)) // strip size of the file <- size of one row of pixels
    || 1 != TIFFSetField(fp, TIFFTAG_SAMPLESPERPIXEL, 3)
    || 1 != TIFFSetField(fp, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB)
    || 1 != TIFFSetField(fp, TIFFTAG_COMPRESSION, COMPRESSION_LZW)
  ) {
    TIFFClose(fp);
    return -1;
  }

  // Write the image one strip at a time
  for (size_t row = 0; row < ny; row++) {
    if (TIFFWriteScanline(fp, (void *)&data_tiff[row * 3 * nx], row, 0) < 0)
      break;
  }

  TIFFClose(fp);

  return 0;
}


/**
 * @brief save four contiguous float arrays into a TIFF file
 *
 * @param fname TIFF file name
 * @param data input array of float values in [0 .. 65535]
 * @param nx ny array size
 *
 * @return 0 if OK, != 0 if an error occured
 */
int write_tiff_rgb_f32(const char *fname, const float *data, size_t nx, size_t ny) {
  const float *ptr_r, *ptr_g, *ptr_b;
  uint16 *data_tiff = NULL;
  uint16 *ptr_out;
  size_t i;
  int retval;

  printf("write_tiff_rgb() -> %f\n", data[8990]);

  /* check allocaton */
  if (NULL == data)
    return -1;

  /* create the tiff array */
  if (NULL == (data_tiff = (uint16*) malloc(3 * nx * ny * sizeof(uint16))))
    return -1;

  ptr_out = data_tiff;


  /* setup the pointers */
  ptr_r = data;
  ptr_g = ptr_r + nx * ny;
  ptr_b = ptr_g + nx * ny;

  /*
   * interlace three arrays (ptr_r, ptr_g, ptr_b)
   * into the TIFF raw array (ptr_out)
   */
  for (i = 0; i < nx * ny; i++) {
    *ptr_out++ = (uint16) (*ptr_r++ + .5);
    *ptr_out++ = (uint16) (*ptr_g++ + .5);
    *ptr_out++ = (uint16) (*ptr_b++ + .5);
  }

  /* write the file */
  retval = write_tiff_rgb_raw(fname, data_tiff, nx, ny);

  free(data_tiff);

  return retval;
}
