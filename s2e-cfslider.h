/* s2e-cfslider.h
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

#if !defined(S2ECFSLIDER_H)
#define S2ECFSLIDER_H

#if defined(__cplusplus)
extern "C" {
#endif

#define SLIDERX 0
#define SLIDERY 1

typedef struct {
   int id;                      /* Handle id */
   float dmin, dmax;            /* Data range */
   float smin, smax;            /* Screen coordinate range */
   float level;                 /* Current value */
   int Nhand;                   /* Number of handles */
   float step;                  /* Step size */
   float width;                 /* Width */
   int direct;                  /* Direction: X/Y */
   float base;                  /* Other screen coordinate */
   char *cmap;                  /* Name of colour map */
   int c1, c2;                  /* Colour index range */
   char *mylabel;
   float last;                  /* Last value of screen coordinates */
   int lower;
   int upper;
   int nstep;
   int toggle; 
   int visible;
} Slider;

Slider initSlider(int id, float dmin, float dmax, float smin, float smax,
		  float base, float width, int c1, int c2,
		  char *cmap, char *ilabel, int lower, int upper,
		  int direct, int nstep, int toggle, int visible);
void drawSlider(Slider sl, float scale, int Ns, int reinstallcmap);

#if defined(__cplusplus)
} // extern "C" {
#endif

#endif
