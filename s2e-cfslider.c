/* s2e-cfslider.cc
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

#include "s2e-cfslider.h"
#include "s2plot.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void drawSlider(Slider sl, float scale, int Ns, int reinstallcmap)
/* Draw an onscreen slider */
{
   float asp = ss2qar();			/* Aspect ratio of window */
   float x, y;
   float dd = (sl.smax - sl.smin)/(float)(Ns-1);
   float r1, g1, b1;
   float dc = (float)(sl.c2-sl.c1)/(float)(Ns-1);
   float ci;
   int i;
   XYZ P[4];
   COLOUR col = { 1.0, 1.0, 1.0 };

   s2scir(sl.c1, sl.c2);			/* Install colour map */
   if (reinstallcmap) {
     s2icm(sl.cmap, sl.c1, sl.c2);
   }

   if (sl.direct == SLIDERX) {			/* X-axis slider */
      P[0].y = P[3].y = sl.base-sl.width/2.0;
      P[1].y = P[2].y = sl.base+sl.width/2.0;
      P[0].z = P[1].z = P[2].z = P[3].z = 0.02;

      for (i=0;i<(Ns-1);i++) {
         x = sl.smin + i*dd;
         ci = sl.c1 + dc*i;
         s2qcr((int)ci, &r1, &g1, &b1);
         col.r = r1; col.g = g1; col.b = b1;
         P[0].x = P[1].x = x;
         P[2].x = P[3].x = x+dd;         
         ns2vf4(P, col);
      }
   } else if (sl.direct == SLIDERY) {		/* Y-axis slider */
      P[0].x = P[3].x = sl.base-sl.width/2.0;
      P[1].x = P[2].x = sl.base+sl.width/2.0;
      P[0].z = P[1].z = P[2].z = P[3].z = 0.02;

     for (i=0;i<(Ns-1);i++) {
         y = sl.smin + i*dd;
         ci = sl.c1 + dc*i;
         s2qcr((int)ci, &r1, &g1, &b1);
         col.r = r1; col.g = g1; col.b = b1;
         P[0].y = P[1].y = y;
         P[2].y = P[3].y = y+dd;
         ns2vf4(P, col);
      }
   }

   col.r = col.g = col.b = 1.0;

/* Write the label */
   dd /= 10.0;
   if (sl.direct == SLIDERX) {
      XYZ pos   = { sl.smax+dd, sl.base-sl.width/4.0, 0.1 };
      XYZ right = { scale, 0.0, 0.0 };
      XYZ norm  = { 0.0, scale*asp, 0.0 };
      ns2vtext(pos, right, norm, col, sl.mylabel);
   } else if (sl.direct == SLIDERY) {
      XYZ pos   = { sl.base+sl.width/4.0, sl.smax+dd, 0.1 };
      XYZ right = { 0.0, scale*asp, 0.0 };
      XYZ norm  = { -scale, 0.0, 0.0 };
      ns2vtext(pos, right, norm, col, sl.mylabel);
   }

}

Slider initSlider(int id, float dmin, float dmax, float smin, float smax,
                float base, float width, int c1, int c2,
                char *cmap, char *ilabel, int lower, int upper,
                int direct, int nstep, int toggle, int visible)
/* Initialise a slider data type */
{
   Slider slider;

   slider.id   = id;		/* ID for this slider */
   slider.dmin = dmin;		/* Data range minimum */
   slider.dmax = dmax;		/* Data range maximum */
   slider.smin = smin;		/* Screen position minimum */
   slider.smax = smax;		/* Screen position maximum */
   slider.base = base;		/* The other coordinate */
   slider.width = width;	/* Width of the slider */
   slider.c1 = c1;		/* Colour index of minimum */
   slider.c2 = c2;		/* Colour index of maximum */
   slider.level = slider.dmin;	/* Set the starting level */
   slider.last = ((slider.level-slider.dmin)/(slider.dmax-slider.dmin))*
                        (slider.smax-slider.smin) + slider.smin;
				/* Set the maximum slider level */

   slider.cmap  = (char *)calloc(255, sizeof(char));	
   sprintf(slider.cmap,"%s",cmap);	/* Colour map name */

   slider.mylabel = (char *)calloc(strlen(ilabel)+1, sizeof(char));
   sprintf(slider.mylabel,"%s",ilabel); 	/* Slider label */

   slider.toggle = toggle;	/* Whether slider has an activation toggle */
   slider.visible = visible;	/* Is the slider visible? */

   slider.lower = lower;	/* Only needed if more than one slider handle */
   slider.upper = upper;
   if ((lower < 0) && (upper < 0)) {
      slider.Nhand = 1;
   } else if (((lower <0) && (upper >= 0)) || ((lower >= 0) && (upper < 0))) {
      slider.Nhand = 2;
   } else {
      slider.Nhand = 3;
   }

   slider.direct = direct;	/* Axis direction: X or Y */
   slider.nstep  = nstep;	/* Number of discrete slider increments */
   slider.step   = (slider.smax-slider.smin)/(float)nstep;	
				/* Size of slider step */
   return slider;
}

