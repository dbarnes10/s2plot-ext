/* s2e-campath.h
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

#include "s2base.h"

/* Camera paths are stored on disk in files like: "base.nnnn.s2v"
 * with "nnnn" zero-padded and running from 0 (first camera position)
 * to (number of camera positions - 1).  A minimum of three positions 
 * are required (ie. 0, 1, 2).  
 *
 * The s2v file format is the exact format that is saved when frames
 * and view files are saved from S2PLOT programs.  So the skeleton
 * of a camera path can be created by saving a number of images from
 * an S2PLOT program then arranging and renaming them to form a 
 * camera path.
 *
 * Ordinarily "base" should be "pathN" with N an integer ... this is 
 * the format that writeCamPath will create when saving camera paths.
 *
 * The camera path starts at the first position, passes through all 
 * subsequent points, ending at the final position.
 *
 * Presently the view position and direction are controlled along the 
 * path; there is no camera roll.
 *
 */

#if defined(__cplusplus)
extern "C" {
#endif

// Structure that holds a camera path.
typedef struct {
  int n;             // number of points
  int *framestamps; // frame numbers
  CAMERA *cameras; // cameras
  float *strength; // strength (for spline fitting at each point)
} S2E_CAMPATH;

// Read a camera path with prefix pname, return the read camera path.
// Note that the framestamps field of the returned camera path is 
// allocated but NOT filled in.  It is up to the caller to do this 
// after reading a path as this information is not stored in the 
// s2v format view files.
S2E_CAMPATH readCamPath(char *pname);

// Write a camera path; return the prefix used for the filenames
char *writeCamPath(S2E_CAMPATH cpath);

// Free a camera path (convenience function)
void freeCamPath(S2E_CAMPATH cpath);

// Find a spline-interpolated position between P1 and P2.  Pm is the
// tangent at P1, Pn is the tangent at P2. Compute these prior to
// calling by e.g. using P0 and P3 if they exist.  Parameter t =
// 0 gives P1, t = 1 gives P2.  strength1 and strength2 give the
// tightness of the spline fit at P1 and P2.  Try values around 1 to
// 15.
XYZ splineInterp(XYZ P1, XYZ P2, XYZ Pm, XYZ Pn, float t,
		 float strength1, float strength2);

// Set a new position (for the given frame number) on the camera path.
// Currently, the UP and
// RIGHT vectors are NOT pathed ... they *are* set, but only following
// this strategy: 
// new right = new viewdir CROSS old updir
// new up = new right CROSS new viewdir
// Because of this, the result of following a camera path will differ
// for differing initial up/right vectors (camera roll).  And it is 
// possible (although untested) that the camera path is not "reversible".
// This will change in the future.
//
// On return, zerocamidx is the index of the camera path element (cpath)
// that marks the start of the segment the camera is currently traversing.
// pathfrac gives the fraction along this segment to the next camera node.
int doCamPath(S2E_CAMPATH cpath, int frameno,
	      int *zerocamidx, float *pathfrac);

// Draw a camera path, optionally with handles.  If handles are used,
// their identifiers will start at handle0id, so set this higher than
// any other handles you are using.
void drawCamPath(S2E_CAMPATH cpath, int addhandles, int handle0id);

// Modify a camera path.  This should really only ever be called from
// a handle drag callback.  It returns 1 if a handle related to this
// camera path was modified ... ie. if the "event" was "consumed".
int modifyCamPath(S2E_CAMPATH cpath, int handle0id, int *id, XYZ *p);

#if defined(__cplusplus)
} // extern "C" {
#endif

