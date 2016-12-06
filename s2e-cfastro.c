/* s2e-cfastro.c
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

#include "s2e-cfastro.h"

#if defined(S2EFITS)

#include "fitsio.h"

/* Type of FITS image to create as a default */
#define TDATA TFLOAT
#define IMG_TYPE FLOAT_IMG
typedef float Imgtype;

void errorFITS(int status)
/* Report an error from dealing with FITS files */
{
   if (status) {
       fits_report_error(stderr, status);
       exit(status);
   }
   return;
}

FITSImage readFITS(char *fname)
/* Read in a FITS file called fname */
{
   FITSImage image;             /* Image which has been read */
   fitsfile *fptr;              /* Pointer to FITS file */
   int status, nfound, anynull;
   long naxes[2], fpixel, npixels, i, j;
   float datamin, datamax, nullval;

   status = 0;                  /* Error condition */

   if (fits_open_file(&fptr, fname, READONLY, &status))
      errorFITS(status);

   if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status))
      errorFITS(status);

   npixels = naxes[0]*naxes[1];         /* Number of pixels */

/* Allocate memory for the image */
   image.array    = (float **)calloc(naxes[1], sizeof(float *));
   image.array[0] = (float *) malloc(naxes[0]*naxes[1]*sizeof(float));
   for (i=1;i<naxes[1];i++) {
      image.array[i] = image.array[i-1] + naxes[0];
   }

   fpixel = 1;
   nullval = 0; /* Don't check for null values in the image */
   datamin =  1.0E30;
   datamax = -1.0E30;

   for (j=0;j<naxes[1];j++) {
      if (fits_read_img(fptr, TFLOAT, fpixel, naxes[0], &nullval,
                image.array[j], &anynull, &status))
         errorFITS(status);

      for (i=0; i<naxes[0]; i++) {
         if (image.array[j][i] < datamin) datamin = image.array[j][i];
         if (image.array[j][i] > datamax) datamax = image.array[j][i];
      }
      fpixel += naxes[0];
   }

   if (fits_close_file(fptr, &status))
      errorFITS(status);

   image.dmin = datamin;
   image.dmax = datamax;
   image.nx   = naxes[0];
   image.ny   = naxes[1];
   image.stdev = 0; 		/* These are not calculated */
   image.mean  = 0;		/* Use statisticsFITS */

   return image;
}


void createFITS(char *fname, int nx, int ny)
/* Create a new FITS file called fanme with dimensions nx * ny */
{
   fitsfile *fptr;
   Imgtype **array;
   int status = 0;

   int bitpix = IMG_TYPE;       /* 16-bit unsigned short pixel values */
   long naxis = 2;              /* 2-dimensional image */
   long naxes[2];               /* width by height */
   int i, j;
   long fpixel, nelements;
   float datamin, datamax;

   naxes[0] = nx;
   naxes[1] = ny;

   array    = (Imgtype **)calloc(ny, sizeof(Imgtype *));
   array[0] = (Imgtype *) malloc(naxes[0]*naxes[1]*sizeof(Imgtype));

   for (i=1;i<naxes[1];i++) {
      array[i] = array[i-1] + naxes[0];
   }
   remove(fname);       /* Delete file if it already exists */

   status = 0;          /* This needs to be done before every fits call */
   if (fits_create_file(&fptr, fname, &status))
      errorFITS(status);

   status = 0;
   if (fits_create_img(fptr, bitpix, naxis, naxes, &status))
      errorFITS(status);

   datamin =  1.0E30;
   datamax = -1.0E30;
   for (j=0;j<naxes[1];j++) {
      for (i=0; i<naxes[0]; i++) {
         array[j][i] = i+j;
         if (array[j][i] < datamin) datamin = array[j][i];
         if (array[j][i] > datamax) datamax = array[j][i];
      }
   }
   fpixel = 1;          /* First pixel */
   nelements = naxes[0] * naxes[1];

   if (fits_write_img(fptr, TDATA, fpixel, nelements, array[0], &status))
      errorFITS(status);

   free(array);

   status = 0;
   if (fits_close_file(fptr, &status) )
      errorFITS(status);

}



void statisticsFITS(FITSImage *image)
{
   float sum, sumsq;            /* Dummy variables for mean and std-dev */
   int i, j;

/* Calculate the mean and std-dev */
   sum   = 0;
   sumsq = 0;
   int count = 0;
   for (j=0;j<image->ny;j++) {
      for (i=0;i<image->nx;i++) {
	if (!isnan(image->array[j][i])) {
	  sum   += image->array[j][i];
	  sumsq += (image->array[j][i]*image->array[j][i]);
	  count++;
	}
      }
   }
   image->mean  = sum/(float)count;
   image->stdev = sqrt((sumsq - image->mean*image->mean)/(float)(count-1));

}

unsigned int us2tidFITS(char *fname, float *asp, char clip)
/* Read in a FITS file and return an s2plot texture ID */
{
   double dd;                   /* Spacing for colour */
   int i, j;                    /* Loop variables */
   unsigned int id;             /* s2plot texture id */
   int width, height;           /* Size of blank texture */
   unsigned char *tex;                   /* Blank texture */
   long idx;                    /* Dummy variable: index into memory */
   unsigned char cgs;           /* Gray scale colour to use */
   FITSImage image;		/* The FITS image */

   image = readFITS(fname);

/* Load the texture, assign an id and copy into a char * memory block */
   id  = ss2ct(image.nx,image.ny);	/* Create a texture big enough for FITS */
   tex = ss2gt(id, &width, &height);	/* Get the memory for this texture */

   dd   = 1.0/(image.dmax - image.dmin);	/* Image data range */

   for (j=0;j<image.ny;j++) {
      for (i=0;i<image.nx;i++) {
         cgs = (char)(255.0*(image.array[j][i]-image.dmin)*dd);

	 if (cgs < clip) cgs = 0;

         idx = (j*width + i) *4;	/* Grey scale */

         tex[ idx     ] = 255;          /* r - cgs*/
         tex[ idx + 1 ] = 255;          /* g - cgs*/
         tex[ idx + 2 ] = 255;          /* b - cgs*/
         tex[ idx + 3 ] = cgs;          /* alpha */

      }
   }

   ss2pt(id);                   /* Restore the texture for future use */

   if (tex != NULL) { free(tex); tex = NULL; }

   *asp = (float)image.ny/(float)image.nx;

   return id;

}


void us2f4tFITSxy(int id, XYZ cen, float scale, float asp, char trans)
{
   XYZ P[4];
   COLOUR col = { 1.0, 1.0, 1.0 };

   float w = scale/2.0;
   float h = w*asp;

   P[0].x = cen.x-w;    P[0].y = cen.y+h;    P[0].z = cen.z;
   P[1].x = cen.x+w;    P[1].y = cen.y+h;    P[1].z = cen.z;
   P[2].x = cen.x+w;    P[2].y = cen.y-h;    P[2].z = cen.z;
   P[3].x = cen.x-w;    P[3].y = cen.y-h;    P[3].z = cen.z;

   ns2vf4x(P, col, id, 1.0, trans);

}


void us2f4tFITSxz(int id, XYZ cen, float scale, float asp, char trans) 
{
   XYZ P[4];
   COLOUR col = { 1.0, 1.0, 1.0 };

   float w = scale/2.0;
   float h = w*asp;

   P[0].x = cen.x-w;    P[0].z = cen.z+h;    P[0].y = cen.y;
   P[1].x = cen.x+w;    P[1].z = cen.z+h;    P[1].y = cen.y;
   P[2].x = cen.x+w;    P[2].z = cen.z-h;    P[2].y = cen.y;
   P[3].x = cen.x-w;    P[3].z = cen.z-h;    P[3].y = cen.y;

   ns2vf4x(P, col, id, 1.0, trans);

}

void us2f4tFITSyz(int id, XYZ cen, float scale, float asp, char trans) 
{
   XYZ P[4];
   COLOUR col = { 1.0, 1.0, 1.0 };

   float w = scale/2.0;
   float h = w*asp;

   P[0].z = cen.z-w;    P[0].y = cen.y+h;    P[0].x = cen.x;
   P[1].z = cen.z+w;    P[1].y = cen.y+h;    P[1].x = cen.x;
   P[2].z = cen.z+w;    P[2].y = cen.y-h;    P[2].x = cen.x;
   P[3].z = cen.z-w;    P[3].y = cen.y-h;    P[3].x = cen.x;

   ns2vf4x(P, col, id, 1.0, trans);
}



void us2f4tFITS(int id, XYZ *P, char trans) 
{
   COLOUR col = { 1.0, 1.0, 1.0 };
   ns2vf4x(P, col, id, 1.0, trans);
}

#endif // #if defined(S2EFITS)


XYZ us2rade2xyz(float ra, float dec, float radius)
/* Converts a RA (0..24) and Dec (-90..+90) to an XYZ point */
{
   XYZ xyz;
   double p, t, rsinp;

   p = (DEMAX - dec)*D2R;             	/* Convert to radians */
   t = ra*RA2R;                   	/* Convert to radians */

   rsinp =  radius * sin(p);           	/* Save one calculation */
   xyz.x =  rsinp  * cos(t);
   xyz.z = -rsinp  * sin(t);
   xyz.y =  radius * cos(p);

   return xyz;
}

void us2radept1(float ra, float dec, float radius, int sym)
{
   XYZ xyz;
   xyz = us2rade2xyz(ra, dec, radius);
   s2pt1(xyz.x, xyz.y, xyz.z, sym);
}

void us2radept(int N, float *ra, float *dec, float radius, int sym)
{
   int i;
   for (i=0;i<N;i++) {
      us2radept1(ra[i], dec[i], radius, sym);
   }
}


#define RASTEPS                  48
#define DESTEPS                  48
#define DEPS                     1.0E-6 /* A small number */
void us2deline(float de, float ramin, float ramax, float radius)
{
   float rr = ramin*RA2R;
   float p1 = radius * sin((DEMAX-de)*D2R);
   float p2 = radius * cos((DEMAX-de)*D2R);
   float d1, d2; 
   float deg = (ramax - ramin)*15.0;
   int nseg = 100;

   XYZ centre = { 0.0, 0.0, 0.0 };
   XYZ normal = { 0.0, radius, 0.0 };
   XYZ start  = { p1, p2, 0.0 };

   /* Rotate the starting vector to ramin */
   d1 =  cos(rr) * start.x; 
   d2 = -sin(rr) * start.x; 
   start.x = d1;
   start.z = d2;

   ns2varc(centre, normal, start, deg, nseg);

}

void us2raline(float ra, float demin, float demax, float radius)
{
   us2radeline(ra, demin, ra, demax, radius);
}

void us2radegrid(int rastep, int destep, float radius)
{
   int i;
   for (i=RAMIN;i<=RAMAX;i+=rastep)
      us2raline((float)i, DEMIN, DEMAX, radius);

   for (i=DEMIN;i<=DEMAX;i+=destep)
      us2deline((float)i, RAMIN, RAMAX, radius);
   
}

void us2radeline(float ra1, float de1, float ra2, float de2, float radius)
{
   XYZ xyz1, xyz2;
   XYZ centre = { 0.0, 0.0, 0.0 };
   XYZ normal;
   float angle;
   int nseg = 100;

   xyz1 = us2rade2xyz(ra1, de1, radius);
   xyz2 = us2rade2xyz(ra2, de2, radius);
   
   normal = cross(xyz1, xyz2);
   angle  = dota(xyz1, xyz2) * R2D;

   ns2varc(centre, normal, xyz1, angle, nseg);
   
} 


void us2radearrow(double ra, double de, double r1, double r2, int ci, int sz)
{  
   double p, t, rsinp; 
   float x[2], y[2], z[2];
   
   p = (90.0 - de)*D2R; 
   t = ra*15.0*D2R;
   
   rsinp = r1 * sin(p);
   x[0] = rsinp * cos(t);
   y[0] = rsinp * sin(t);
   z[0] = r1 * cos(p); 

   rsinp = r2 * sin(p);
   x[1] = rsinp * cos(t);
   y[1] = rsinp * sin(t); 
   z[1] = r2 * cos(p);
   
   s2sci(ci);
   s2sch(sz);
   s2arro(x[0],y[0],z[0],x[1],y[1],z[1]);
   s2sch(1);
}

XYZ cross(XYZ a, XYZ b)
{
   XYZ c;
   c.x = a.y*b.z - a.z*b.y;
   c.y = a.z*b.x - a.x*b.z;
   c.z = a.x*b.y - a.y*b.x;

   return c;
}

//#define EPS 1.0E-6
XYZ unit(XYZ a)
{
   XYZ b;
   float r;
   r = length(a);
   if (fabs(r) < DEPS) {
      b.x = b.y = b.z = 0;
   } else {
      b.x = a.x/r;
      b.y = a.y/r;
      b.z = a.z/r;
   }
   return b;
}

float dotp(XYZ a, XYZ b)
{
   return (a.x*b.x + a.y*b.y + a.z*b.z);
}

float length(XYZ a)
{
   return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

float dota(XYZ a, XYZ b)
{
   XYZ a1 = unit(a);
   XYZ b1 = unit(b);
   float d = dotp(a1, b1);
   return acos(d);
}


int countLines(FILE *fp)
{
   int N = 0;
   long int ptr;
   char string[1024];

   ptr = ftell(fp);

   fgets(string,1024,fp);
   while (!feof(fp)) { fgets(string,1024,fp); N++; }
   fseek(fp, ptr, SEEK_SET);            /* Back to where file was */

   return N;
}



float **initSlice(int nx, int ny)
/* Allocate memory and initialise a data cube */
{
   float **slice;
   int i, j;

   slice = (float **)malloc(nx * sizeof(float *));
   if (slice == NULL) {
      fprintf(stderr,"Failed to allocate %ld bytes\n",nx*sizeof(float *));
      exit(-1);
   }
   for (i=0;i<nx;i++) {
      slice[i] = (float *)malloc(ny * sizeof(float));
      if (slice[i] == NULL) {
         fprintf(stderr,"Failed to allocate %ld bytes\n",nx*sizeof(float));
         exit(-1);
      }
      for (j=0;j<ny;j++) {
         slice[i][j] = 0.0;
      }
   }

   return slice;
}

void equalizeVolume(float ***volume, int M, float vmin, float vmax)
{  
   int i, j, k;
   float *pdf, *cdf;
   int maxn; 
   int idx;

   maxn = (int)(vmax+1.0);
   pdf = (float *)calloc(maxn, sizeof(float));
   cdf = (float *)calloc(maxn, sizeof(float));
   for (i=0;i<maxn;i++) {
      pdf[i] = 0; 
      cdf[i] = 0;
   }        
            
   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         for (k=0;k<M;k++) {
            idx = (int)round(volume[i][j][k]);
            pdf[idx]+=1.0;
         }
      }
   }
   cdf[0] = pdf[0];
   for (i=1;i<maxn;i++) {
      cdf[i] = cdf[i-1] + pdf[i];
   }
   float alpha = maxn/(float)(M*M*M);

   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         for (k=0;k<M;k++) {
            idx = (int)round(volume[i][j][k]);
            volume[i][j][k] = cdf[idx]*alpha;
         }
      }
   }
}

void equalizeSlice(float **slice, int M, float smin, float smax)
{  
   int i, j;
   float *pdf, *cdf;
   int maxn; 
   int idx;

   maxn = (int)(smax+1.0);
   pdf = (float *)calloc(maxn, sizeof(float));
   cdf = (float *)calloc(maxn, sizeof(float));
   for (i=0;i<maxn;i++) {
      pdf[i] = 0; 
      cdf[i] = 0;
   }        
            
   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         idx = (int)round(slice[i][j]);
         pdf[idx]+=1.0;
      }
   }
   cdf[0] = pdf[0];
   for (i=1;i<maxn;i++) {
      cdf[i] = cdf[i-1] + pdf[i];
   }
   float alpha = maxn/(float)(M*M);

   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         idx = (int)round(slice[i][j]);
         slice[i][j] = cdf[idx]*alpha;
      }
   }
}

void maskVolume(float ***volume, int M, float clip, float value)
{
   int i, j, k;

   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         for (k=0;k<M;k++) {
            if (volume[i][j][k] < clip) volume[i][j][k] = value;
         }
      }
   }
}

void maskSlice(float **slice, int M, float clip, float value)
{
   int i, j;

   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         if (slice[i][j] < clip) slice[i][j] = value;
      }
   }
}

void normaliseSlice(float **slice, int M, float smin, float smax, float scale)
{
   int i, j;

   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         slice[i][j] = scale*(slice[i][j]-smin)/(smax-smin);
      }
   }
}

void normaliseVolume(float ***volume, int M, float vmin, float vmax, 
		float scale)
{
   int i, j, k;

   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         for (k=0;k<M;k++) {
            volume[i][j][k] = scale*(volume[i][j][k]-vmin)/(vmax-vmin);
         }
      }
   }
}

void normaliseVolume3(float ***volume, int Mx, int My, int Mz, 
		float vmin, float vmax, float scale)
{
   int i, j, k;

   for (i=0;i<Mx;i++) {
      for (j=0;j<My;j++) {
         for (k=0;k<Mz;k++) {
            volume[i][j][k] = scale*(volume[i][j][k]-vmin)/(vmax-vmin);
         }
      }
   }
}

void log10Volume(float ***volume, int M)
{
   int i, j, k;

   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         for (k=0;k<M;k++) {
            volume[i][j][k] = log10(volume[i][j][k]+1.0);
         }
      }
   }
}

void log10Slice(float **slice, int M)
{
   int i, j;

   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         slice[i][j] = log10(slice[i][j]+1.0);
      }
   }
}


XYZ minmaxVolumeCell(float ***volume, int M, float *rvmin, float *rvmax)
{
   int i, j, k;
   float vmin, vmax;
   XYZ cell;

   vmin = volume[0][0][0];
   vmax = volume[0][0][0];
   cell.x = cell.y = cell.z = 0;
   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         for (k=0;k<M;k++) {
            if (volume[i][j][k] < vmin) { 
               vmin = volume[i][j][k];
            }
            if (volume[i][j][k] > vmax) {
               vmax = volume[i][j][k];
	       cell.x = i;
	       cell.y = j;
	       cell.z = k;
            }
         }
      }
   }
   *rvmin = vmin;
   *rvmax = vmax;
   return cell;
}

void minmaxVolume(float ***volume, int M, float *rvmin, float *rvmax)
{
   int i, j, k;
   float vmin, vmax;

   vmin = volume[0][0][0];
   vmax = volume[0][0][0];
   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         for (k=0;k<M;k++) {
            vmin = (volume[i][j][k] < vmin) ? volume[i][j][k] : vmin;
            vmax = (volume[i][j][k] > vmax) ? volume[i][j][k] : vmax;
         }
      }
   }
   *rvmin = vmin;
   *rvmax = vmax;
}

void minmaxSlice(float **slice, int M, float *rsmin, float *rsmax)
{
   int i, j;
   float smin, smax;

   smin = slice[0][0];
   smax = slice[0][0];
   for (i=0;i<M;i++) {
      for (j=0;j<M;j++) {
         smin = (slice[i][j] < smin) ? slice[i][j] : smin;
         smax = (slice[i][j] > smax) ? slice[i][j] : smax;
      }
   }
   *rsmin = smin;
   *rsmax = smax;
}



float ***initVolume(int nx, int ny, int nz)
/* Allocate memory and initialise a data cube */
{
   float ***volume;
   int i, j, k;

   volume = (float ***)malloc(nx * sizeof(float **));
   if (volume == NULL) {
      fprintf(stderr,"Failed to allocate %ld bytes\n",nx*sizeof(float **));
      exit(-1);
   }
   for (i=0;i<nx;i++) {
      volume[i] = (float **)malloc(ny * sizeof(float *));
      if (volume[i] == NULL) {
         fprintf(stderr,"Failed to allocate %ld bytes\n",nx*sizeof(float *));
         exit(-1);
      }
      for (j=0;j<ny;j++) {
         volume[i][j] = (float *)malloc(nz * sizeof(float));
         if (volume[i][j] == NULL) {
            fprintf(stderr,"Failed to allocate %ld bytes\n",nx*sizeof(float));
            exit(-1);
         }
         for (k=0;k<nz;k++) {
            volume[i][j][k] = 0.0;              /* Initialise */
         }
      }
   }

   return volume;
}



float hsl2rgb(float T, float S, float L)
/* 
   Hue is between 0.0 and 360.0 
   Call as:   
      R = hsl2rgb(hue/360.0 + 1.0/3.0, S, L);
      G = hsl2rgb(hue/360.0          , S, L);
      B = hsl2rgb(hue/360.0 - 1.0/3.0, S, L);
*/
{
   float Q, P;

   if (L< 0.5) Q = L * (1.0 + S);
   else Q = L + S - (L*S);
   P = 2.0*L - Q;
   if (T < 0) T = T + 1.0;
   if (T < (1.0/6.0)) { return (P + ((Q-P)*6.0*T)); }
   if (T < 0.5) { return Q; }
   if (T < (2.0/3.0)) { return (P + ((Q-P)*(2.0/3.0 - T) * 6.0)); }
   return P;
}
