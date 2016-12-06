/* s2e-cfsmooth.c
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

#include "s2e-cfsmooth.h"


float ngp2d(float x1, float x2, float y1, float y2, int M, XYZ *xyz, int N, float **grid)
{
   int i, j, k;
   XYZ cen, diff;
   float dx, dy;
   float cw = 0;
   float wx, wy, w;
   int jmin, jmax;
   int kmin, kmax;
   int gx, gy;
   int match = 0;

   dx = (x2-x1)/M;
   dy = (y2-y1)/M;
   for (i=0;i<N;i++) {
      gx = (xyz[i].x-x1)/dx;            /* Guess at which grid cell is centre */
      gy = (xyz[i].y-y1)/dy;
      jmin = gx-2; if (jmin < 0) jmin = 0;      /* Only search near this */
      jmax = gx+2; if (jmax > M) jmax = M;
      kmin = gy-2; if (kmin < 0) kmin = 0;
      kmax = gy+2; if (kmax > M) kmax = M;
      for (j=jmin;j<jmax;j++) {
         cen.x = (j+0.5)*dx + x1;               /* Cell centre: x */
         diff.x = fabs(xyz[i].x - cen.x);       /* Distance to cell centre: x */
         wx = 0.0;                              /* Reset weight */
         if (diff.x/dx <= 0.5) wx = 1.0;        /* Weight condition */
         for (k=kmin;k<kmax;k++) {
            cen.y = (k+0.5)*dy + y1;            /* Cell centre: y */
            diff.y = fabs(xyz[i].y - cen.y);    /* Distance to cell centre: y */
            wy = 0;                             /* Reset weight */
            if (diff.y/dx <= 0.5) wy = 1.0;     /* Weight condition */
            w = wx*wy;
            if (w) {
               match++;                         /* Found something */
               cw += w;                         /* Cumulative weight */
               grid[j][k] += w;
            }
         }
      }
   }
   if (match != N) { fprintf(stderr,"NGP missed %d\n", N-match); }
   return cw;
}

float cic2d(float x1, float x2, float y1, float y2, int M, XYZ *xyz, int N, float **grid)
{
   int i, j, k;
   XYZ cen, diff;
   float dx, dy;
   float wx, wy, w;
   float ddx, ddy;
   int match = 0;
   int jmin, jmax;
   int kmin, kmax;
   int gx, gy;
   float cw = 0;

   dx = (x2-x1)/M;
   dy = (y2-y1)/M;
   for (i=0;i<N;i++) {
      gx = (xyz[i].x - x1)/dx;          /* Guess at which grid cell is centre */
      gy = (xyz[i].y - y1)/dy;          /* Guess at which grid cell is centre */
      jmin = gx-2; if (jmin < 0) jmin = 0;      /* Only search near this */
      jmax = gx+2; if (jmax > M) jmax = M;
      kmin = gy-2; if (kmin < 0) kmin = 0;
      kmax = gy+2; if (kmax > M) kmax = M;

      for (j=jmin;j<jmax;j++) {
         cen.x = (j+0.5)*dx + x1;               /* Centre of the current cell: x */
         diff.x = fabs(xyz[i].x - cen.x);       /* Distance to centre of current: x */
         wx = 0.0;                              /* Reset weight */
         ddx = diff.x/dx;
         if (ddx <= 1.0) { wx = 1.0 - ddx; };   /* Weight condition */
         for (k=kmin;k<kmax;k++) {
            cen.y = (k+0.5)*dy + y1;            /* Centre of the current cell: y */
            diff.y = fabs(xyz[i].y - cen.y);    /* Distance to centre of current: x */
            wy = 0;                             /* Reset weight */
            ddy = diff.y/dy;
            if (ddy <= 1.0) { wy = 1.0 - ddy; };        /* Weight condition */
            w = wx*wy;                          /* Weight for this cell */
            if (w) {                            /* Non-zero weight */
               if (((j-gx) == 0) && ((k-gy) == 0)) {
                   match++;                     /* Check that something was found */
               }
               grid[j][k] += w;                 /* Increase weight of this cell */
               cw += w;
            }
         }
      }
   }
   if (match != N) { fprintf(stderr,"CIC missed %d\n", N-match); }
   return cw;
}

float tsc2d(float x1, float x2, float y1, float y2, int M, XYZ *xyz, int N, float **grid)
{
   int i, j, k;
   XYZ cen, diff;
   float dx, dy;
   float wx, wy, w;
   float ddx, ddy;
   float cw = 0;
   int match;
   int jmin,jmax;
   int kmin,kmax;
   int gx,gy;

   dx = (x2-x1)/M;
   dy = (y2-y1)/M;
   match = 0;
   for (i=0;i<N;i++) {
      gx = (xyz[i].x-x1)/dx;            /* Guess at which grid cell is centre */
      gy = (xyz[i].y-y1)/dy;
      jmin = gx-2; if (jmin < 0) jmin = 0;      /* Only search near this */
      jmax = gx+2; if (jmax > M) jmax = M;
      kmin = gy-2; if (kmin < 0) kmin = 0;
      kmax = gy+2; if (kmax > M) kmax = M;
      for (j=jmin;j<jmax;j++) {
         cen.x = (j+0.5)*dx + x1;               /* Cell centre: x */
         diff.x = fabs(xyz[i].x - cen.x);       /* Distance to cell centre: x */
         wx = 0.0;                              /* Reset weight */
         ddx = diff.x/dx;
         if (ddx <= 0.5) { wx = 0.75 - ddx*ddx; }
         else if ((ddx > 0.5) && (ddx <= 1.5))  /* Weight conditions */
                { wx = 0.5*(1.5 - ddx)*(1.5-ddx); }
         for (k=kmin;k<kmax;k++) {
            cen.y = (k+0.5)*dy + y1;            /* Cell centre: y */
            diff.y = fabs(xyz[i].y - cen.y);    /* Distance to cell centre: y */
            wy = 0;                             /* Reset weight */
            ddy = diff.y/dy;
            if (ddy <= 0.5) { wy = 0.75 - ddy*ddy; }
            if ((ddy > 0.5) && (ddy <= 1.5))    /* Weight condition */
                { wy = 0.5*(1.5 - ddy)*(1.5-ddy); }
            w = wx*wy;
            if (w) {
               if (((j-gx) == 0) && ((k-gy) == 0)) {    /* Was something found */
                   match++;
               }
               grid[j][k] += w;
               cw += w;
            }
         }
      }
   }
   if (match != N) { fprintf(stderr,"TSC missed %d\n", N-match); }
   return cw;
}

float ngp3d(float x1, float x2, float y1, float y2, float z1, float z2,
                int M, XYZ *xyz, int N, float ***grid)
{
   int i, j, k, m;
   XYZ cen, diff;
   float dx, dy, dz;
   float cw = 0;
   float wx, wy, wz, w;
   int jmin, jmax;
   int kmin, kmax;
   int mmin, mmax;
   int gx, gy, gz;
   int match = 0;


   dx = (x2-x1)/M;
   dy = (y2-y1)/M;
   dz = (z2-z1)/M;
   fprintf(stderr,"%f %f %f %f %f %f\n",x2,x1,y2,y1,z2,z1);
   fprintf(stderr,"%f %f %f\n",dx,dy,dz);
   for (i=0;i<N;i++) {
      gx = (xyz[i].x-x1)/dx;            /* Guess at which grid cell is centre */
      gy = (xyz[i].y-y1)/dy;
      gz = (xyz[i].z-z1)/dz;
      jmin = gx-2; if (jmin < 0) jmin = 0;      /* Only search near this */
      jmax = gx+2; if (jmax > M) jmax = M;
      kmin = gy-2; if (kmin < 0) kmin = 0;
      kmax = gy+2; if (kmax > M) kmax = M;
      mmin = gz-2; if (mmin < 0) mmin = 0;
      mmax = gz+2; if (mmax > M) mmax = M;
      for (j=jmin;j<jmax;j++) {
         cen.x = (j+0.5)*dx + x1;               /* Cell centre: x */
         diff.x = fabs(xyz[i].x - cen.x);       /* Distance to cell centre: x */
         wx = 0.0;                              /* Reset weight */
         if (diff.x/dx <= 0.5) wx = 1.0;        /* Weight condition */
         for (k=kmin;k<kmax;k++) {
            cen.y = (k+0.5)*dy + y1;            /* Cell centre: y */
            diff.y = fabs(xyz[i].y - cen.y);    /* Distance to cell centre: y */
            wy = 0;                             /* Reset weight */
            if (diff.y/dx <= 0.5) wy = 1.0;     /* Weight condition */
            for (m=mmin;m<mmax;m++) {
               cen.z = (m+0.5)*dz + z1;         /* Cell centre: z */
               diff.z = fabs(xyz[i].z - cen.z); /* Distance to cell centre: y */
               wz = 0;                          /* Reset weight */
               if (diff.z/dz <= 0.5) wz = 1.0;  /* Weight condition */
               w = wx*wy*wz;
               if (w) {
                  match++;                              /* Found something */
                  cw += w;                              /* Cumulative weight */
                  grid[j][k][m] += w;
               }
            }
         }
      }
   }
   if (match != N) { fprintf(stderr,"NGP missed %d\n", N-match); }
   fprintf(stderr,"NGP done\n");
   return cw;
}


XYZ *ngp3dRGB(float x1, float x2, float y1, float y2, float z1, float z2,
                int M, XYZ *xyz, int N, float ***grid)
{
   XYZ *cell;
   int i, j, k, m;
   XYZ cen, diff;
   float dx, dy, dz;
   float cw = 0;
   float wx, wy, wz, w;
   int jmin, jmax;
   int kmin, kmax;
   int mmin, mmax;
   int gx, gy, gz;
   int match = 0;

   cell = (XYZ *)calloc(N, sizeof(XYZ));

   dx = (x2-x1)/M;
   dy = (y2-y1)/M;
   dz = (z2-z1)/M;
   for (i=0;i<N;i++) {
      gx = (xyz[i].x-x1)/dx;            /* Guess at which grid cell is centre */
      gy = (xyz[i].y-y1)/dy;
      gz = (xyz[i].z-z1)/dz;
      jmin = gx-2; if (jmin < 0) jmin = 0;      /* Only search near this */
      jmax = gx+2; if (jmax > M) jmax = M;
      kmin = gy-2; if (kmin < 0) kmin = 0;
      kmax = gy+2; if (kmax > M) kmax = M;
      mmin = gz-2; if (mmin < 0) mmin = 0;
      mmax = gz+2; if (mmax > M) mmax = M;
      for (j=jmin;j<jmax;j++) {
         cen.x = (j+0.5)*dx + x1;               /* Cell centre: x */
         diff.x = fabs(xyz[i].x - cen.x);       /* Distance to cell centre: x */
         wx = 0.0;                              /* Reset weight */
         if (diff.x/dx <= 0.5) wx = 1.0;        /* Weight condition */
         for (k=kmin;k<kmax;k++) {
            cen.y = (k+0.5)*dy + y1;            /* Cell centre: y */
            diff.y = fabs(xyz[i].y - cen.y);    /* Distance to cell centre: y */
            wy = 0;                             /* Reset weight */
            if (diff.y/dx <= 0.5) wy = 1.0;     /* Weight condition */
            for (m=mmin;m<mmax;m++) {
               cen.z = (m+0.5)*dz + z1;         /* Cell centre: z */
               diff.z = fabs(xyz[i].z - cen.z); /* Distance to cell centre: y */
               wz = 0;                          /* Reset weight */
               if (diff.z/dz <= 0.5) wz = 1.0;  /* Weight condition */
               w = wx*wy*wz;
               if (w) {
                  match++;                              /* Found something */
                  cw += w;                              /* Cumulative weight */
		  cell[i].x = (float)j;
		  cell[i].y = (float)k;
		  cell[i].z = (float)m;
               } 
            }
         }
      }
   }
   if (match != N) { fprintf(stderr,"NGP missed %d\n", N-match); }
   return cell;
}


float cic3d(float x1, float x2, float y1, float y2, float z1, float z2,
                int M, XYZ *xyz, int N, float ***grid)
{
   int i, j, k, m;
   XYZ cen, diff;
   float dx, dy, dz;
   float wx, wy, wz, w;
   float ddx, ddy, ddz;
   int match = 0;
   int jmin, jmax;
   int kmin, kmax;
   int mmin, mmax;
   int gx, gy, gz;
   float cw = 0;

   dx = (x2-x1)/M;
   dy = (y2-y1)/M;
   dz = (z2-z1)/M;
   for (i=0;i<N;i++) {
      gx = (xyz[i].x - x1)/dx;          /* Guess at which grid cell is centre */
      gy = (xyz[i].y - y1)/dy;          /* Guess at which grid cell is centre */
      gz = (xyz[i].z - z1)/dz;          /* Guess at which grid cell is centre */
      jmin = gx-2; if (jmin < 0) jmin = 0;      /* Only search near this */
      jmax = gx+2; if (jmax > M) jmax = M;
      kmin = gy-2; if (kmin < 0) kmin = 0;
      kmax = gy+2; if (kmax > M) kmax = M;
      mmin = gz-2; if (mmin < 0) mmin = 0;
      mmax = gz+2; if (mmax > M) mmax = M;

      for (j=jmin;j<jmax;j++) {
         cen.x = (j+0.5)*dx + x1;               /* Centre of the current cell: x */
         diff.x = fabs(xyz[i].x - cen.x);       /* Distance to centre of current: x */
         wx = 0.0;                              /* Reset weight */
         ddx = diff.x/dx;
         if (ddx <= 1.0) { wx = 1.0 - ddx; };   /* Weight condition */
         for (k=kmin;k<kmax;k++) {
            cen.y = (k+0.5)*dy + y1;            /* Centre of the current cell: y */
            diff.y = fabs(xyz[i].y - cen.y);    /* Distance to centre of current: x */
            wy = 0;                             /* Reset weight */
            ddy = diff.y/dy;
            if (ddy <= 1.0) { wy = 1.0 - ddy; };        /* Weight condition */
            for (m=mmin;m<mmax;m++) {
               cen.z = (m+0.5)*dz + z1;         /* Centre of the current cell: y */
               diff.z = fabs(xyz[i].z - cen.z); /* Distance to centre of current: x */
               wz = 0;                          /* Reset weight */
               ddz = diff.z/dz;
               if (ddz <= 1.0) { wz = 1.0 - ddz; };     /* Weight condition */
               w = wx*wy*wz;                            /* Weight for this cell */
               if (w) {                         /* Non-zero weight */
                  if (((j-gx) == 0) && ((k-gy) == 0) && ((m-gz) == 0)) {
                     match++;                   /* Check that something was found */
                  }
                  grid[j][k][m] += w;           /* Increase weight of this cell */
                  cw += w;
               }
            }
         }
      }
   }
   if (match != N) { fprintf(stderr,"CIC missed %d\n", N-match); }
   return cw;
}

float tsc3d(float x1, float x2, float y1, float y2, float z1, float z2,
                int M, XYZ *xyz, int N, float ***grid)
{
   int i, j, k, m;
   XYZ cen, diff;
   float dx, dy, dz;
   float wx, wy, wz, w;
   float ddx, ddy, ddz;
   float cw = 0;
   int match;
   int jmin,jmax;
   int kmin,kmax;
   int mmin,mmax;
   int gx,gy,gz;

   dx = (x2-x1)/M;
   dy = (y2-y1)/M;
   dz = (z2-z1)/M;
   match = 0;
   for (i=0;i<N;i++) {
      gx = (xyz[i].x-x1)/dx;            /* Guess at which grid cell is centre */
      gy = (xyz[i].y-y1)/dy;
      gz = (xyz[i].z-z1)/dz;
      jmin = gx-2; if (jmin < 0) jmin = 0;      /* Only search near this */
      jmax = gx+2; if (jmax > M) jmax = M;
      kmin = gy-2; if (kmin < 0) kmin = 0;
      kmax = gy+2; if (kmax > M) kmax = M;
      mmin = gz-2; if (mmin < 0) mmin = 0;
      mmax = gz+2; if (mmax > M) mmax = M;
      for (j=jmin;j<jmax;j++) {
         cen.x = (j+0.5)*dx + x1;               /* Cell centre: x */
         diff.x = fabs(xyz[i].x - cen.x);       /* Distance to cell centre: x */
         wx = 0.0;                              /* Reset weight */
         ddx = diff.x/dx;
         if (ddx <= 0.5) { wx = 0.75 - ddx*ddx; }
         else if ((ddx > 0.5) && (ddx <= 1.5))  /* Weight conditions */
                { wx = 0.5*(1.5 - ddx)*(1.5-ddx); }
         for (k=kmin;k<kmax;k++) {
            cen.y = (k+0.5)*dy + y1;            /* Cell centre: y */
            diff.y = fabs(xyz[i].y - cen.y);    /* Distance to cell centre: y */
            wy = 0;                             /* Reset weight */
            ddy = diff.y/dy;
            if (ddy <= 0.5) { wy = 0.75 - ddy*ddy; }
            if ((ddy > 0.5) && (ddy <= 1.5))    /* Weight condition */
                { wy = 0.5*(1.5 - ddy)*(1.5-ddy); }
            for (m=mmin;m<mmax;m++) {
               cen.z = (m+0.5)*dz + z1;         /* Cell centre: z */
               diff.z = fabs(xyz[i].z - cen.z); /* Distance to cell centre: z */
               wz = 0;                          /* Reset weight */
               ddz = diff.z/dz;
               if (ddz <= 0.5) { wz = 0.75 - ddz*ddz; }
               if ((ddz > 0.5) && (ddz <= 1.5))         /* Weight condition */
                { wz = 0.5*(1.5 - ddz)*(1.5-ddz); }
               w = wx*wy*wz;
               if (w) {
                  if (((j-gx) == 0) && ((k-gy) == 0) & ((m-gz) == 0)) {
                                                /* Was something found */
                      match++;
                  }
                  grid[j][k][m] += w;
                  cw += w;
              }
            }
         }
      }
   }
   if (match != N) { fprintf(stderr,"TSC missed %d\n", N-match); }
   return cw;
}

XYZ *tsc3dRGB(float x1, float x2, float y1, float y2, float z1, float z2,
                int M, XYZ *xyz, int N, float ***grid)
{
   XYZ *cell;
   int i, j, k, m;
   XYZ cen, diff;
   float dx, dy, dz;
   float wx, wy, wz, w;
   float ddx, ddy, ddz;
   float cw = 0;
   int jmin,jmax;
   int kmin,kmax;
   int mmin,mmax;
   int gx,gy,gz;
   int match = 0;

   cell = (XYZ *)calloc(N, sizeof(XYZ));
   

   dx = (x2-x1)/M;
   dy = (y2-y1)/M;
   dz = (z2-z1)/M;
   for (i=0;i<N;i++) {
      gx = (xyz[i].x-x1)/dx;            /* Guess at which grid cell is centre */
      gy = (xyz[i].y-y1)/dy;
      gz = (xyz[i].z-z1)/dz;
      jmin = gx-2; if (jmin < 0) jmin = 0;      /* Only search near this */
      jmax = gx+2; if (jmax > M) jmax = M;
      kmin = gy-2; if (kmin < 0) kmin = 0;
      kmax = gy+2; if (kmax > M) kmax = M;
      mmin = gz-2; if (mmin < 0) mmin = 0;
      mmax = gz+2; if (mmax > M) mmax = M;
      for (j=jmin;j<jmax;j++) {
         cen.x = (j+0.5)*dx + x1;               /* Cell centre: x */
         diff.x = fabs(xyz[i].x - cen.x);       /* Distance to cell centre: x */
         wx = 0.0;                              /* Reset weight */
         ddx = diff.x/dx;
         if (ddx <= 0.5) { wx = 0.75 - ddx*ddx; }
         else if ((ddx > 0.5) && (ddx <= 1.5))  /* Weight conditions */
                { wx = 0.5*(1.5 - ddx)*(1.5-ddx); }
         for (k=kmin;k<kmax;k++) {
            cen.y = (k+0.5)*dy + y1;            /* Cell centre: y */
            diff.y = fabs(xyz[i].y - cen.y);    /* Distance to cell centre: y */
            wy = 0;                             /* Reset weight */
            ddy = diff.y/dy;
            if (ddy <= 0.5) { wy = 0.75 - ddy*ddy; }
            if ((ddy > 0.5) && (ddy <= 1.5))    /* Weight condition */
                { wy = 0.5*(1.5 - ddy)*(1.5-ddy); }
            for (m=mmin;m<mmax;m++) {
               cen.z = (m+0.5)*dz + z1;         /* Cell centre: z */
               diff.z = fabs(xyz[i].z - cen.z); /* Distance to cell centre: z */
               wz = 0;                          /* Reset weight */
               ddz = diff.z/dz;
               if (ddz <= 0.5) { wz = 0.75 - ddz*ddz; }
               if ((ddz > 0.5) && (ddz <= 1.5))         /* Weight condition */
                { wz = 0.5*(1.5 - ddz)*(1.5-ddz); }
               w = wx*wy*wz;
               if (w) {
                  if (((j-gx) == 0) && ((k-gy) == 0) & ((m-gz) == 0)) {
                                                /* Was something found */
                      match++;
                      cell[i].x = (float)j;
                      cell[i].y = (float)k;
                      cell[i].z = (float)m;
                  }
                  grid[j][k][m] += w;
                  cw += w;
              }
            }
         }
      }
   }
   if (match != N) { fprintf(stderr,"TSC missed %d\n", N-match); }
   return cell;
}



