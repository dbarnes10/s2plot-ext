/* s2e-campath.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "paulslib.h"
#include "s2plot.h"
#include "s2e-campath.h"

void ReadViewFile(char *);
void SaveViewFile(CAMERA, char *);
XYZ _s2_world2device(XYZ world);
XYZ _s2_world2device_so(XYZ world);
XYZ _s2_device2world(XYZ device);
extern CAMERA _s2x_camera;

S2E_CAMPATH readCamPath(char *pname) {

  // make new camera path
  S2E_CAMPATH cpath;
  cpath.n = 0;
  cpath.framestamps = NULL;
  cpath.cameras = NULL;
  cpath.strength = NULL;

  // read as many files as there are and store in cpath
  char cname[255];
  int ic = 0;
  sprintf(cname, "%s.%04d.s2v", pname, ic);
  FILE *fptr = fopen(cname, "r");
  while (fptr != NULL) {
    fclose(fptr);
    ReadViewFile(cname);
    cpath.cameras = (CAMERA *)realloc(cpath.cameras, (ic + 1)*sizeof(CAMERA));
    bcopy((char *)&_s2x_camera, (char *)&cpath.cameras[ic], sizeof(CAMERA));
    ic++;
    sprintf(cname, "%s.%04d.s2v", pname, ic);
    fptr = fopen(cname, "r");
  }
  cpath.n = ic;

  cpath.framestamps = (int *)malloc(cpath.n * sizeof(int));
  cpath.strength = (float *)malloc(cpath.n * sizeof(float));
  for (ic = 0; ic < cpath.n; ic++) {
    cpath.framestamps[ic] = ic*10;
    cpath.strength[ic] = 0.33;
  }

  // return the read path
  return cpath;
}

char *writeCamPath(S2E_CAMPATH cpath) {
  static int pathnum = 0;

  char *cname = (char *)calloc(255, sizeof(char *));
  sprintf(cname, "path%d.0000.s2v", pathnum);
  FILE *fptr = fopen(cname, "r");
  while (fptr != NULL) {
    pathnum++;
    fclose(fptr);
    sprintf(cname, "path%d.0000.s2v", pathnum);
    fptr = fopen(cname, "r");
  }
    
  int ic;
  for (ic = 0; ic < cpath.n; ic++) {
    sprintf(cname, "path%d.%04d.s2v", pathnum, ic);
    SaveViewFile(cpath.cameras[ic], cname);
  }

  fprintf(stderr, "Saved %d camera positions to files \"path%d.nnnn.s2v\".\n",
	  cpath.n, pathnum);

  // incr. pathnum for next save
  pathnum++;

  // return the prefix name used
  sprintf(cname, "path%d", pathnum);
  return cname;
}

void freeCamPath(S2E_CAMPATH cpath) {
  if (cpath.framestamps) {
    free(cpath.framestamps);
  }
  cpath.framestamps = NULL;
  if (cpath.cameras) {
    free(cpath.cameras);
  }
  cpath.cameras = NULL;
  if (cpath.strength) {
    free(cpath.strength);
  }
  cpath.strength = NULL;
  cpath.n = 0;
}

// Find a spline-interpolated position between P1 and P2.  Pm is the
// tangent at P1, Pn is the tangent at P2. Compute these prior to
// calling by e.g. using P0 and P3 if they exist.  Parameter t =
// 0 gives P1, t = 1 gives P2.  strength1 and strength2 give the
// tightness of the spline fit at P1 and P2.  Try values around 1 to
// 15.
XYZ splineInterp(XYZ P1, XYZ P2, XYZ Pm, XYZ Pn, float t,
		 float strength1, float strength2) {

  XYZ Pt, P1a, P2a;
  float ll1;
  XYZ result;
  float ax, bx, cx;
  float strengthScale = 1.0;

  // 0. get displacement vector P1->P2
  Pt = VectorSub(P1, P2);
  strengthScale = Modulus(Pt);
  Normalise(&Pt);

  // 1. get P1a from P1 + 1/3 tangent projected on P1->P2
  ll1 = fabs(DotProduct(Pm, Pt));
  ll1 = (ll1 < strength1*strengthScale) ? strength1*strengthScale : ll1;
  SetVectorLength(&Pm, 0.333 * ll1);
  P1a = VectorAdd(P1, Pm);

  // 2. get P2a from P2 + 1/3 tangent projected on P1->P2
  ll1 = fabs(DotProduct(Pn, Pt));
  ll1 = (ll1 < strength2*strengthScale) ? strength2*strengthScale : ll1;
  SetVectorLength(&Pn, 0.333 * ll1);
  P2a = VectorAdd(P2, Pn);

  // 3. calculate new position:
  cx = 3. * (P1a.x - P1.x);
  bx = 3. * (P2a.x - P1a.x) - cx;
  ax = P2.x - P1.x - cx - bx;
  result.x = ax * t * t * t + bx * t * t + cx * t + P1.x;
  
  cx = 3. * (P1a.y - P1.y);
  bx = 3. * (P2a.y - P1a.y) - cx;
  ax = P2.y - P1.y - cx - bx;
  result.y = ax * t * t * t + bx * t * t + cx * t + P1.y;
  
  cx = 3. * (P1a.z - P1.z);
  bx = 3. * (P2a.z - P1a.z) - cx;
  ax = P2.z - P1.z - cx - bx;
  result.z = ax * t * t * t + bx * t * t + cx * t + P1.z;

  return result;
}

// Set a position on the camera path.
int doCamPath(S2E_CAMPATH cpath, int frameno,
	      int *zerocamidx, float *pathfrac) {
  XYZ P0, P3, Pm, Pn;
  float llt;
  XYZ newcampos, ncvd, ncup, cright;
  int wh = 0;
  
  // 1. find previous and next frame
  wh = 0;
  while (wh < cpath.n-1 && frameno >= cpath.framestamps[wh]) {
    wh++;
  }
  wh--;
  
  // 2. get start and end points of bezier
  P0 = cpath.cameras[wh].vp;
  P3 = cpath.cameras[wh+1].vp;
  
  // 3a. get tangent at P0: it is either line P-1 to P3, or P0 to P3
  //     if P-1 doesn't exist
  if (wh > 0) {
    Pm = VectorSub(cpath.cameras[wh-1].vp, P3);
  } else {
    Pm = VectorSub(P0, P3);
  }
  
  // 4a. get tangent at P3: it is either line P4 to P0, or P3 to P0
  //     if P4 doesn't exist
  if (wh+1 < cpath.n-1) {
    Pn = VectorSub(cpath.cameras[wh+2].vp, P0);
  } else {
    Pn = VectorSub(P3, P0);
  }
  
  // 5. get parameter along the curve
  llt = (float)(frameno - cpath.framestamps[wh]) / 
    (float)(cpath.framestamps[wh+1] - cpath.framestamps[wh]);
  
  newcampos = splineInterp(P0, P3, Pm, Pn, llt, //OVERSHOOT, OVERSHOOT);
			   cpath.strength[wh],
			   cpath.strength[wh+1]);
  
  // REPEAT FOR VIEW DIRECTION
  // 2. get start and end points of bezier
  P0 = cpath.cameras[wh].vd;
  P3 = cpath.cameras[wh+1].vd;
  
  // 3a. get tangent at P0: it is either line P-1 to P3, or P0 to P3
  //     if P-1 doesn't exist
  if (wh > 0) {
    Pm = VectorSub(cpath.cameras[wh-1].vd, P3);
  } else {
    Pm = VectorSub(P0, P3);
  }
  
  // 4a. get tangent at P3: it is either line P4 to P0, or P3 to P0
  //     if P4 doesn't exist
  if (wh+1 < cpath.n-1) {
    Pn = VectorSub(cpath.cameras[wh+2].vd, P0);
  } else {
    Pn = VectorSub(P3, P0);
  }
  
  ncvd = splineInterp(P0, P3, Pm, Pn, llt, 0.35, 0.35);
  
  Normalise(&ncvd);
  
  /* determine new up vector */
  ss2qc(NULL, &ncup, NULL, 0);
  cright = CrossProduct(ncvd, ncup);
  ncup = CrossProduct(cright, ncvd);
  Normalise(&ncup);
  
  ss2sc(newcampos, ncup, ncvd, 0);

  *zerocamidx = wh;
  *pathfrac = llt;

  return 1;
}


// Draw a camera path, optionally with handles.
void drawCamPath(S2E_CAMPATH cpath, int addhandles, int handle0id) {

  int frame;
  XYZ P0, P3, Pm, Pn;
  float llt;
  XYZ newcampos;
  
  COLOUR cyan = {0., 1., 1.};
  COLOUR red = {1., 0., 0.};
  COLOUR yellow = {1., 1., 0.};

  XYZ vdv;
  for (frame = 0; frame < cpath.n; frame++) {
    // draw campos in cyan dot, size increasing with frame number
    ns2vthpoint(_s2_device2world(cpath.cameras[frame].vp), cyan, 2+frame);

    if (addhandles) {
      // add handle at campos
      cs2ah(_s2_device2world(cpath.cameras[frame].vp), 200., 
	    cyan, cyan, handle0id + frame, 0);
    }

    // draw line in direction of viewdir in yellow
    vdv = VectorAdd(cpath.cameras[frame].vp,
		    VectorMul(cpath.cameras[frame].vd,
			      0.06));
    ns2vthline(_s2_device2world(cpath.cameras[frame].vp), 
	       _s2_device2world(vdv), 
	       yellow, 2.);
    if (addhandles) { 
      // and a handle
      cs2ah(_s2_device2world(vdv), 200., 
	    yellow, yellow, handle0id + cpath.n +frame, 0);
    }

  }
  
  for (frame = 0; frame <= cpath.framestamps[cpath.n-1]; frame++) {
    
    // 1. find previous and next frame
    int wh = 0;
    while (wh < cpath.n-1 && frame >= cpath.framestamps[wh]) {
      wh++;
    }
    wh--;
    
    // 2. get start and end points of bezier
    P0 = cpath.cameras[wh].vp;
    P3 = cpath.cameras[wh+1].vp;
    
    // 3a. get tangent at P0: it is either line P-1 to P3, or P0 to P3
    //     if P-1 doesn't exist
    if (wh > 0) {
      Pm = VectorSub(cpath.cameras[wh-1].vp, P3);
    } else {
      Pm = VectorSub(P0, P3);
    }
    
    // 4a. get tangent at P3: it is either line P4 to P0, or P3 to P0
    //     if P4 doesn't exist
    if (wh+1 < cpath.n-1) {
      Pn = VectorSub(cpath.cameras[wh+2].vp, P0);
    } else {
      Pn = VectorSub(P3, P0);
    }
    
    // 4. parameter along the curve:
    llt = (float)(frame - cpath.framestamps[wh]) /
      (float)(cpath.framestamps[wh+1] - cpath.framestamps[wh]);
    
    //newcampos = splineInterp(P0, P3, Pm, Pn, llt, OVERSHOOT, OVERSHOOT);
    newcampos = splineInterp(P0, P3, Pm, Pn, llt, 
			     cpath.strength[wh], 
			     cpath.strength[wh+1]);
    
    ns2vpoint(_s2_device2world(newcampos), red);
  } 
}

int modifyCamPath(S2E_CAMPATH cpath, int handle0id, int *id, XYZ *p) {
  int consumed = 0;
  if (*id >= handle0id) {
    if (*id < handle0id + cpath.n) {
      // cam position moved
      cpath.cameras[*id-handle0id].vp = _s2_world2device(*p);
      consumed = 1;
    } else if (*id < handle0id + 2 * cpath.n) {
      // cam view direction changed
      XYZ tmp = VectorSub(cpath.cameras[*id-handle0id-cpath.n].vp, 
			  _s2_world2device(*p));
      Normalise(&tmp);
      cpath.cameras[*id-handle0id-cpath.n].vd = tmp;
      consumed = 1;
    }
  }
  return consumed;
}
