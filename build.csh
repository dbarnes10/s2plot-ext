#!/bin/csh
#
# build.csh
#
# Copyright 2012 David G. Barnes
#
# This file is part of S2PLOT-EXT.
#
# S2PLOT-EXT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# S2PLOT-EXT is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with S2PLOT-EXT.  If not, see <http://www.gnu.org/licenses/>. 
#
# We would appreciate it if research outcomes using S2PLOT-EXT would
# provide the following acknowledgement:
#
# "Three-dimensional visualisation was conducted with the S2PLOT
# progamming library"
#
# and a reference to
#
# D.G.Barnes, C.J.Fluke, P.D.Bourke & O.T.Parry, 2006, Publications
# of the Astronomical Society of Australia, 23(2), 82-93.
#

if (!(${?S2EXTRAINC})) then
  setenv S2EXTRAINC ""
endif
if (!(${?S2EXTRALIB})) then
  setenv S2EXTRALIB ""
endif

clbuild.csh s2ext s2e-campath s2e-cfastro s2e-cfsmooth s2e-cfslider

