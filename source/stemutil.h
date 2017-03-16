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

#ifndef STEMUTIL_H
#define STEMUTIL_H

#include "data_containers.h"
#include "stemtypes_fftw3.h"

#define pNPIX		0   /* number of pix 1 for real and 2 for complex */
#define	pRMAX		1	/* maximum value of the real part of the image */
#define	pIMAX		2	/* maximum value of the imaginary part of the image */
#define pRMIN		3	/* minimum value of the real part of the image */
#define	pIMIN		4	/* minimum value of the imaginary part of the image */
#define pXBTILT		5	/* x beam tilt in rad */
#define pYBTILT		6	/* y beam tilt in rad */
#define pC		7	/* c unit cell dimension in Angstroms */
#define pRES		8	/* real space resolution in atompot */
#define pXCTILT		9	/* x crystal tilt in rad */
#define pYCTILT		10	/* y crystal tilt in rad */
#define pDEFOCUS	11	/* defocus in Angstroms */
#define pASTIG		12	/* astigmatism in Angstroms */
#define pTHETA		13	/* angle of astigmatism in radians */
#define pDX		14	/* dimension of pixel in x direction in Angstroms */
#define pDY		15	/* dimension of pixel in y direction in Angstroms */
#define pENERGY		16	/* beam energy in keV */
#define pOAPERT		17	/* objective aperture semi-angle in radians */
#define pCS		18	/* spherical aberration in Angstroms */
#define pWAVEL		19	/* electron wavelength in Angstroms */
#define pCAPERT		21	/* condenser (CTEM) illumination angle in radians */
#define pDDF		22	/* defocus spread in Angstroms */
#define pNSLICES	29	/* number of slices */
#define pMINDET		31	/* minimum detector angle (STEM) in radians */
#define pMAXDET		32	/* maximum detector angle (STEM) in radians */

typedef struct timeval timev;
typedef struct timezone timez;


int writeCFG(atom *atoms,int natoms,char *fileName, MULS *muls);
int phononDisplacement(double *u,MULS *muls,int id,int icx,int icy,
		       int icz,int atomCount,double dw,int maxAtom, int Znum);

void *memcopy(void *dest, const void *src, size_t n);
// void saveSTEMimages(MULS *muls);
// atom *readUnitCell(int *natom,char *fileName,MULS *muls,int handleVacancies);
atom *readCFGUnitCell(int *natom,char *fileName,MULS *muls);

double v3DatomLUT(int iz,double r,int tdsFlag,int scatFlag);
double vzatomLUT( int Z, double r ,int tdsFlag,int scatFlag);
double v3DzatomLUT(int Znum, real r2D, real z);  // real

double wavelength( double kev );
double v3Datom(int Z, double r,int tdsFlag,int scatFlag);
double vzatom( int Z, double radius,int tdsFlag,int scatFlag);

double gasdev(long *idum);
double rangauss( unsigned long *iseed );
int ReadfeTable(int scatFlag );
int ReadLine( FILE* fpRead, char* cRead, int cMax, const char *mesg );
int getZNumber(char *element);
double sigma( double kev );
void splinh( double x[], double y[],
	     double b[], double c[], double d[], int n);
double seval( double *x, double *y, double *b, double *c,
	     double *d, int n, double x0 );
double ranflat( unsigned long *iseed );
int parlay( const char c[], int islice[], int nsmax, int lmax,
	    int *nslice, int fperr );
/* long powerof2( long n ); */
double fe3D(int Z, double q2,int tdsFlag, double scale,int scatFlag);
double sfLUT(double s,int atKind, MULS *muls);
double bicubic(double **ff,int Nz, int Nx,double z,double x);
int atomCompare(const void *atom1,const void *atom2);


double getTime();
double cputim();



#endif
