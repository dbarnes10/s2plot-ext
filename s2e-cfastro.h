/* s2e-cfastro.h
 * 
 * Copyright 2008- Christopher Fluke, David G Barnes
 *
 * This file is part of S2PLOT-EXT.
 *
 * S2PLOT-EXT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * S2PLOT-EXT is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with S2PLOT-EXT.  If not, see <http://www.gnu.org/licenses/>.
 *
 * We would appreciate it if research outcomes using S2PLOT-EXT would
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

#ifndef __S2ASTRO_H__
#define __S2ASTRO_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


#include "s2plot.h"

#if defined(__cplusplus)
extern "C" {
#endif


/* CONSTANTS */
#define D2R  0.017453292519943295       /* Conversion of degrees to radians */
#define RA2R 0.261799387799149408	/* Conversion of RA (hrs) to radians */
#define R2D 57.295779513082322865       /* Conversion of radians to degrees */

#define DEMIN			-90.0		/* Maximum Declination in degrees */
#define DEMAX                   +90.0   	/* Maximum Declination in degrees */
#define RAMIN                     0.0   	/* Minimum RA value in hours */
#define RAMAX                    24.0   	/* Maximum RA value in hours */


//#define S2EFITS
#if defined(S2EFITS)
#include "fitsio.h"

/* DATA STRUCTURES */

typedef struct {
   float **array;
   float dmin, dmax;
   float mean, stdev;
   int nx, ny;
} FITSImage;


/* FUNCTION PROTOTYPES */


/* Create a default FITS file called fname with dimensions nx * ny */
void createFITS(char *fname, int nx, int ny);

/* Read in a FITS file called fname */
FITSImage readFITS(char *fname);

/* Calculate image statistics for a FITS file */
void statisticsFITS(FITSImage *image);

/* Create an s2plot texture ID for a FITS file. Set anything below clip to 0 */
unsigned int us2tidFITS(char *fname, float *asp, char clip);

/* Draw a polygon with central position and aspect ratio asp, width is */
/* scale and transparency (trans) is one of 'o', 's', 't' */
void us2f4tFITSxy(int id, XYZ centre, float scale, float asp, char trans);
void us2f4tFITSxz(int id, XYZ centre, float scale, float asp, char trans);
void us2f4tFITSyz(int id, XYZ centre, float scale, float asp, char trans);
void us2f4tFITS(int id, XYZ *P, char trans);

#endif // #if defined(S2EFITS)

/* Coordinate functions */

/* Convert RA (hours) and Declination (degrees) into XYZ coordinates */
/*   for a sphere of given radius */
XYZ us2rade2xyz(float ra, float dec, float radius);

/* Plot N pairs of RA (hours) and Declination (degrees) points */
/*  Sphere has radius r, plot using symbol sym */
void us2radept(int N, float *ra, float *dec, float radius, int sym);

/* Plot a single RA (hours) and Declination (degrees) point */
/*  Sphere has radius r, plot using symbol sym */
void us2radept1(float ra, float dec, float radius, int sym);

/* Draw a line at constant RA (hours) between Declinations (degrees) */
/*   demin and demax. Sphere has radius r */
void us2raline(float ra, float demin, float demax, float r);

/* Draw a line at constant Declination (degrees) between RAs (hours) */
/*   ramin and ramax. Sphere has radius r */
void us2deline(float de, float ramin, float ramax, float r);

/* Draw a RA (hour) and Declination (degrees) grip with integer steps in */
/*  RA and Declination.  Sphere has radius r */
void  us2radegrid(int rastep, int destep, float radius);

/* Draw an arrow between radius r1 and radius r2 for constant RA (hour) */
/*   and Declination (degrees). Arrow has colour ci and size sz  */
void us2radearrow(double ra, double de, double r1, double r2, int ci, int sz);

XYZ unit(XYZ a);
float dotp(XYZ a, XYZ b);
float length(XYZ a);
float dota(XYZ a, XYZ b);
XYZ cross(XYZ a, XYZ b);
void us2radeline(float ra1, float de1, float ra2, float de2, float radius);

int countLines(FILE *fp);

float ***initVolume(int nx, int ny, int nz);
void equalizeVolume(float ***volume, int M, float vmin, float vmax);
void log10Volume(float ***volume, int M);
void normaliseVolume(float ***volume, int M, float vmin, float vmax,
	float scale);
void minmaxVolume(float ***volume, int M, float *rvmin, float *rvmax);
XYZ minmaxVolumeCell(float ***volume, int M, float *rvmin, float *rvmax);
void maskVolume(float ***volume, int M, float clip, float value);

void maskSlice(float **slice, int M, float clip, float value);
void minmaxSlice(float **slice, int M, float *rsmin, float *rsmax);
void normaliseSlice(float **slice, int M, float smin, float smax, float scale);
void log10Slice(float **slice, int M);
float **initSlice(int nx, int ny);
void equalizeSlice(float **slice, int M, float smin, float smax);
float hsl2rgb(float T, float S, float L);

#if defined(__cplusplus)
} // extern "C" {
#endif

#endif

