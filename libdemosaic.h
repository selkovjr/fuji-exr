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


#ifndef _LIBDEMOSAIC_H_
#define _LIBDEMOSAIC_H_



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "libAuxiliary.h"

/**
 * @file   libdemosaic.cpp
 * @brief  Demosaicking functions: HAmilton-Adams algorithm, NLmeans based demosaicking, Chromatic components filtering
 *
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
  int origHeight,
  unsigned char* mask
);


/**
 * \brief  Classical bilinear interpolation of red and blue differences with the green channel
 *
 *
 * @param[in]  ored, ogreen, oblue  original cfa image with green interpolated
 * @param[out] ored, ogreen, oblue  demosaicked output
 * @param[in]  width, height size of the image
 *
 */
void bilinear_red_blue(
  float *ored,
  float *ogreen,
  float *oblue,
  int width,
  int height,
  int origWidth,
  int origHeight,
  unsigned char* mask
);



/**
 * \brief  NLmeans based demosaicking
 *
 * For each value to be filled, a weigthed average of original CFA values of the same channel is performed.
 * The weight depends on the difference of a 3x3 color patch
 *
 * @param[in]  ired, igreen, iblue  initial demosaicked image
 * @param[out] ored, ogreen, oblue  demosaicked output
 * @param[in]  bloc  research block of size (2+bloc+1) x (2*bloc+1)
 * @param[in]  h kernel bandwidth
 * @param[in]  width, height size of the image
 *
 */


void demosaic_nlmeans(
  int bloc,
  float h,
  float *ored,
  float *ogreen,
  float *oblue,
  float *ired,
  float *igreen,
  float *iblue,
  int width,
  int height,
  int origWidth,
  int origHeight,
  unsigned char* mask
);



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



void chromatic_median(int iter,int ,int projflag,float side,float *ired,float *igreen, float *iblue,float *ored,float *ogreen,float *oblue,int width,int height);




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
 * @param[in]  width, height size of the image
 *
 */


void ssdd_demosaic_chain(
  float *ired,
  float *igreen,
  float *iblue,
  float *ored,
  float *ogreen,
  float *oblue,
  int width,
  int height,
  int origWidth,
  int origHeight,
  unsigned char* mask
);

#endif
