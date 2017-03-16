/*
QSTEM - image simulation for TEM/STEM/CBED
    Copyright (C) 2000-2010  Christoph Koch
	Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef STEMTYPES_H
#define STEMTYPES_H

#include "memory_fftw3.h"
//#include "boost/shared_ptr.hpp"

////////////////////////////////////////////////////////////////////////
// define whether to use single or double precision
///////////////////////////////////////////////////////////////////////
#define FLOAT_PRECISION 1


#define BW (2.0F/3.0F)	/* bandwidth limit */
#define DOYLE_TURNER 0
#define WEICK_KOHL 1
#define CUSTOM 2
#define STEM    1
#define CBED    2
#define TEM     3
#define REFINE  4
#define MSCBED  5
#define TOMO    6
#define NBED    7

////////////////////////////////////////////////////////////////////////
// Define physical constants
////////////////////////////////////////////////////////////////////////
#define ELECTRON_CHARGE (1.6021773e-19)
#define PICO_AMPERE (1e-12/ELECTRON_CHARGE)
#define MILLISEC_PICOAMP (1e-3*PICO_AMPERE)

// #include "floatdef.h"
#include "fftw3.h"

////////////////////////////////////////////////////////////////
#if FLOAT_PRECISION == 1
#define fftw_real float
#ifndef float_tt
#define float_tt  float
#endif
#define real      float
#else  // FLOAT_PRECISION
#define fftw_real double
#ifndef float_tt
#define float_tt  float
#endif
#define real      float
#endif  // FLOAT_PRECISION
////////////////////////////////////////////////////////////////

typedef struct atomStruct {
  float z,y,x;
  // float dx,dy,dz;  // thermal displacements
  float dw;      // Debye-Waller factor
  float occ;     // occupancy
  float q;       // charge 
  int Znum;
} atom;

/* Planes will be defined by the standard equation for a plane, i.e.
 * a point (point) and 2 vectors (vect1, vect2)
 */
typedef struct planeStruct {
  double normX,normY,normZ;
  double vect1X,vect1Y,vect1Z;
  double vect2X,vect2Y,vect2Z;
  double pointX,pointY,pointZ;
} plane;

typedef struct grainBoxStruct {
  int amorphFlag;
  double density,rmin, rFactor;  /* density, atomic distance, reduced atomic distance
				  * for amorphous material.  Red. r is for making a hex.
				  * closed packed structure, which will fill all space, 
				  * but will only be sparsely filled, and later relaxed.
				  */ 
  char *name;
  atom *unitCell; /* definition of unit cell */
  int natoms;     /* number of atoms in unit cell */
  double ax,by,cz; /* unit cell parameters */
  double alpha, beta, gamma; /* unit cell parameters */
  double tiltx,tilty,tiltz;
  double shiftx,shifty,shiftz;
  plane *planes;   /* pointer to array of bounding planes */
  double sphereRadius, sphereX,sphereY,sphereZ; /* defines a sphere instead of a grain with straight edges */
  int nplanes; /* number of planes in array planes */
} grainBox;

typedef struct superCellBoxStruct {
  double cmx,cmy,cmz;  /* fractional center of mass coordinates */
  double ax,by,cz;
  int natoms;
  atom *atoms; /* contains all the atoms within the super cell */
} superCellBox;

typedef struct atomBoxStruct {
  int used;   /* indicate here whether this atom is used in the
		 particular problem */
  int nx,ny,nz;
  float_tt dx,dy,dz;
  double B;
#if FLOAT_PRECISION == 1
  fftwf_complex ***potential;   /* 3D array containg 1st quadrant of real space potential */
  float_tt ***rpotential;   /* 3D array containg 1st quadrant of real space potential */
#else
  fftw_complex ***potential;   /* 3D array containg 1st quadrant of real space potential */
  float_tt ***rpotential;   /* 3D array containg 1st quadrant of real space potential */
#endif
} atomBox;

#endif // STEMTYPES_H
