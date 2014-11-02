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


/**
 * @file io_tiff.cpp
 * @brief TIFF I/O routines
 *
 * @author Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
 *
 * @todo stdin/stdout handling
 * @todo TIFF float version
 *
 * These routines read raw sensor data as a 16-bit grayscale file
 * and write out the demosaicked file as a 16-bit RGB. The I/O data
 * are internally represented as concatenated float arrays.
 */

#include <stdlib.h>
#include <tiffio.h>


float *read_tiff_gray16_f32(const char *fname, size_t *nx, size_t *ny);
int write_tiff_rgb_f32(const char *fname, const float *data, size_t nx, size_t ny);




