/* s2e-cfsmooth.h
 *
 * Copyright 2006-2012 David G. Barnes, Paul Bourke, Christopher Fluke
 *
 * This file is part of S2PLOT.
 *
 * S2PLOT is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * S2PLOT is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with S2PLOT.  If not, see <http://www.gnu.org/licenses/>. 
 *
 * We would appreciate it if research outcomes using S2PLOT would
 * provide the following acknowledgement:
 *
 * "Three-dimensional visualisation was conducted with the S2PLOT
 * progamming library"
 *
 * and a reference to
 *
 * D.G.Barnes, C.J.Fluke, P.D.Bourke & O.T.Parry, 2006, Publications
 * of the Astronomical Society of Australia, 23(2), 82-93.
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "s2e-cfastro.h"

#if defined(__cplusplus)
extern "C" {
#endif


float cic2d(float x1, float x2, float y1, float y2, int M, XYZ *xyz, int N, float **grid);
float ngp2d(float x1, float x2, float y1, float y2, int M, XYZ *xyz, int N, float **grid);
float tsc2d(float x1, float x2, float y1, float y2, int M, XYZ *xyz, int N, float **grid);


float ngp3d(float x1, float x2, float y1, float y2, float z1, float z2,
                int M, XYZ *xyz, int N, float ***grid);
float cic3d(float x1, float x2, float y1, float y2, float z1, float z2,
                int M, XYZ *xyz, int N, float ***grid);
float tsc3d(float x1, float x2, float y1, float y2, float z1, float z2,
                int M, XYZ *xyz, int N, float ***grid);

XYZ *ngp3dRGB(float x1, float x2, float y1, float y2, float z1, float z2,
                int M, XYZ *xyz, int N, float ***grid);


#if defined(__cplusplus)
} // extern "C" {
#endif
