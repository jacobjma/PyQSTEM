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

/* file: STEMutil.c */

#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#ifndef WIN32
//#include <sys/time.h>
#endif

#include "stemlib.h"
#include "stemutil.h"
#include "memory_fftw3.h"	/* memory allocation routines */
// #include "tiffsubs.h"
#include "matrixlib.h"
#include "readparams.h"
#include "fileio_fftw3.h"
// #include "fileio.h"
// #include "floatdef.h"

// #define WRITEK

#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif


#define K_B 1.38062e-23      /* Boltzman constant */
#define PID 3.14159265358979 /* pi */
#define SQRT_PI 1.77245385090552 /* sqrt(pi) */
#define HPLANCK 6.6262e-34   /* planck's constant */
#define AMU 1.66053e-27      /* atomic mass unit */
#define THZ_AMU_PI_HBAR 0.05012012918415 /*  A�^2*THz*amu/(pi*hbar) */
#define THZ_AMU_HBAR 0.15745702964189    /*   A�^2*THz*amu/(hbar)   */
// 4.46677327584453 /* 1e10/sqrt(THz*amu/(pi*hbar)) */
#define THZ_HBAR_2KB  3.81927135604119     /* THz*hbar/(2*kB) */
#define THZ_HBAR_KB   1.90963567802059     /* THz*hbar/kB */
#define AMU_THZ2_A2_KB   1.20274224623720     /* AMU*THz^2*A^2/kB */

#define NCINMAX  500	/* max number of characers in stacking spec */
#define NRMAX	50	/* number of values in look-up-table in vzatomLUT */
#define RMIN	0.01	/* min r (in Ang) range of LUT for vzatomLUT() */
#define RMAX	5

#define NPDTMAX 8       /* number of parameters for doyle turner sfacts */
#define NPMAX	12	/* number of parameters for each Z */
#define NZMIN	1	/* min Z in featom.tab */
#define NZMAX	98      /* max Z (in featom.dat ZMAX=103) */
#define	NCMAX	132	/* characters per line to read */
#define NPARAM	64	/* number of parameters in tiff files */
int feTableRead=0;	/* flag to remember if the param file has been read */
double **fparams=NULL;	/* to get fe(k) parameters */
const int nl=3, ng=3;	/* number of Lorenzians and Gaussians */


/******************************************************
 * memcopy replaces memcpy, because memcpy seems to have a bug!
 */
void *memcopy(void *dest, const void *src, size_t n)
{
  int i;

  for (i=0;i<n;i++)
    ((char *)dest)[i] = ((char *)src)[i];

  return dest;
  /* memcpy(dest,src,n) */
}

/*****************************************************************
 * v3DatomLUT(int Z, double r)
 * returns 3D potential (not projected) at radius r
 * (using Lookup table)
 ****************************************************************/
double v3DatomLUT(int Z,double r,int tdsFlag,int scatFlag)
{
  int i,iz;
  double dlnr;
  static double *splinr = NULL;
  static double **splinv = NULL;
  static double **splinb,**splinc,**splind;
  static int *nspline = NULL;

  iz = Z-1;

  if(splinr == NULL) {
    splinr = double1D(NRMAX,"splinr");
    splinv = (double **)malloc(NZMAX*sizeof(double*));
    splinb = (double **)malloc(NZMAX*sizeof(double*));
    splinc = (double **)malloc(NZMAX*sizeof(double*));
    splind = (double **)malloc(NZMAX*sizeof(double*));
    nspline = (int *)malloc( NZMAX*sizeof(int));
    for( i=0; i<NZMAX; i++)
      nspline[i] = 0;

    /*  generate a set of logarithmic r values */
    dlnr = log(RMAX/RMIN)/(NRMAX-1);
    for( i=0; i<NRMAX; i++)
      splinr[i] = RMIN * exp( i * dlnr );
    /*
      printf( "fit from r= %g to r= %g\n", splinr[0], splinr[NRMAX-1] );
    */

    nspline = (int *)malloc( NZMAX*sizeof(int));
    for( i=0; i<NZMAX; i++)
      nspline[i] = 0;
  }

  /* if this atomic number has not been called before
     generate the spline coefficients */
  if( nspline[iz] == 0 ) {
    /*
    printf("generating 3D spline %d \n",Z);
    */
    splinv[iz] = double1D( NRMAX, "splinv" );
    splinb[iz] = double1D( NRMAX, "splinb" );
    splinc[iz] = double1D( NRMAX, "splinc" );
    splind[iz] = double1D( NRMAX, "splind" );

    for( i=0; i<NRMAX; i++)
      splinv[iz][i] = v3Datom(Z,splinr[i],tdsFlag,scatFlag);
    nspline[iz] = NRMAX;
    splinh( splinr, splinv[iz], splinb[iz],
	    splinc[iz], splind[iz], NRMAX);
  }

  /***************************************************
   * now that everything is set up find the
   * scattering factor by interpolation in the table
   * We can either do the real interpolation (lookup <0)
   * or simply look up the potential in the table
   **************************************************/
  return seval(splinr, splinv[iz], splinb[iz],
	       splinc[iz], splind[iz], nspline[iz], r );
}  /* end v3DatomLUT() */


/*--------------------- vzatomLUT() -----------------------------------*/
/*
	return the (real space) projected atomic potential for atomic
	number Z at radius r (in Angstroms)

	this mimics vzatom() in slicelib.c but uses a look-up-table
	with cubic spline interpolatoin to make it run about 2X-4X faster

	started 23-may-1997 E. Kirkland
	fix Z range to allow Hydrogen 1-jan-1998 ejk

	Z = atomic number 1 <= iz <= 98
	r = radius in Angstroms
*/

double vzatomLUT(int Z, double r,int tdsFlag,int scatFlag)
{
  int i,iz;
  double dlnr,result;
  static double *splinr = NULL;
  static double **splinv = NULL;
  static double **splinb,**splinc,**splind;
  static int *nspline = NULL;

  iz = Z-1;

  if(splinr == NULL) {
    splinr = double1D(NRMAX,"splinr");
    splinv = (double **)malloc(NZMAX*sizeof(double*));
    splinb = (double **)malloc(NZMAX*sizeof(double*));
    splinc = (double **)malloc(NZMAX*sizeof(double*));
    splind = (double **)malloc(NZMAX*sizeof(double*));
    nspline = (int *)malloc( NZMAX*sizeof(int));
    for( i=0; i<NZMAX; i++)
      nspline[i] = 0;

    /*  generate a set of logarithmic r values */
    dlnr = log(RMAX/RMIN)/(NRMAX-1);
    for( i=0; i<NRMAX; i++)
      splinr[i] = RMIN * exp( i * dlnr );
    printf( "fit from r= %gA to r= %gA\n", splinr[0], splinr[NRMAX-1] );

    nspline = (int *)malloc( NZMAX*sizeof(int));
    for( i=0; i<NZMAX; i++)
      nspline[i] = 0;
  }

  /* if this atomic number has not been called before
     generate the spline coefficients */
  if( nspline[iz] == 0 ) {
    /*
    printf("generating spline %d\n",Z);
    */
    splinv[iz] = double1D( NRMAX, "splinv" );
    splinb[iz] = double1D( NRMAX, "splinb" );
    splinc[iz] = double1D( NRMAX, "splinc" );
    splind[iz] = double1D( NRMAX, "splind" );

    for( i=0; i<NRMAX; i++)
      splinv[iz][i] = vzatom(Z,splinr[i],tdsFlag,scatFlag);
    nspline[iz] = NRMAX;
    splinh( splinr, splinv[iz], splinb[iz],
	    splinc[iz], splind[iz], NRMAX);
  }

  /***************************************************
   * now that everything is set up find the
   * scattering factor by interpolation in the table
   * We can either do the real interpolation (lookup <0)
   * or simply look up the potential in the table
   **************************************************/
  result = seval(splinr, splinv[iz], splinb[iz],
	       splinc[iz], splind[iz], nspline[iz], r );

  return result;
}  /* end vzatomLUT() */




/* This function saves the STEM images for all detectors in float tiff files
 * the name will have the sample thickness appended to it.
 */
/*
void saveSTEMimages(MULS *muls) {
  int i,ix,iy,nx,ny,npix;
  real rmin,rmax,x;
  float param[NPARAM];
  char fileName[32],datetime[32];
  static float **image=NULL;
  static int nx_old=0,ny_old=0;
  double t;

  nx = (*muls).scanXN;
  ny = (*muls).scanYN;

  if ((nx !=nx_old) || (ny!=ny_old) || (image == NULL)) {
    if (image != NULL) {
      free(image[0]);
      free(image);
    }
    image = (float **)malloc(nx*sizeof(float *));
    image[0] = (float *)malloc(nx*ny*sizeof(float));
    for (ix=0;ix<nx;ix++) image[ix] = &(image[0][ix*ny]);
    nx_old = nx;
    ny_old = ny;
  }

  //  printf("will save STEM images for %d detectors\n",(*muls).detectorNum);
  for (i=0;i<(*muls).detectorNum;i++) {
    sprintf(fileName,"./%s/%s_%d.tif",(*muls).folder,(*muls).detectors[i].name,
	    (int)((*muls).thickness+0.5));

    /**********************************************************
     * if this is not the first averaging run, we want to read
     * the average of previous runs first and include this one
     **********************************************************
    if ((*muls).avgCount > 0) {
      // printf("will open %s now (avgCount=%d)\n",fileName,(*muls).avgCount);
      if( topenFloat(fileName) != 1 ) {
	printf("Cannot open floating point image %s (change to .img format)\n",fileName);
	for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) image[ix][iy] = 0.0;
      }
      else {
	if(treadFloatPix(image, (long)nx,(long)ny,
			  &npix, datetime,param) != 1) {
	  printf("Cannot read input file %s\n",fileName);
	  for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) image[ix][iy] = 0.0;
	}
	//printf("will be closing now\n");
	tclose();
	// printf("did close\n");
      }
    } // end of if muls.avgCount > 0..
    else
      for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) image[ix][iy] = 0.0;

    (*muls).detectors[i].error = 0.0;
    for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
      t = ((double)((*muls).avgCount)*image[ix][iy] +
	       (double)((*muls).detectors[i].image[ix][iy])) /
	((double)((*muls).avgCount+1.0));
      (*muls).detectors[i].error += (image[ix][iy]-t)*(image[ix][iy]-t);
      image[ix][iy] = t;
    }
    (*muls).detectors[i].error = sqrt((*muls).detectors[i].error/(nx*ny));

    rmin = image[0][0];
    rmax = rmin;
    for( ix=0; ix<nx; ix++)  for( iy=0; iy<ny; iy++) {
      x =  image[ix][iy];
      if( x < rmin ) rmin = x;
      if( x > rmax ) rmax = x;
    }

    param[pRMAX]  = rmax;
    param[pIMAX]  = 0.0f;
    param[pRMIN]  = rmin;
    param[pIMIN]  = 0.0f;
    param[pXCTILT] = (*muls).ctiltx;
    param[pYCTILT] = (*muls).ctiltx;
    param[pENERGY] = (*muls).v0;
    param[pDX] = ((*muls).scanXStop-(*muls).scanXStart)/(real)nx;
    param[pDY] = ((*muls).scanYStop-(*muls).scanYStart)/(real)ny;
    param[pWAVEL] = (real)wavelength((*muls).v0);
    param[pNSLICES] = (real)((*muls).nslic0);
    param[pDEFOCUS] = (*muls).df0;
    param[pOAPERT] = 0.0;
    param[pCS] = (*muls).Cs;
    param[pCAPERT] = (*muls).alpha;
    param[pDDF] = 0.0;
    param[pC]= (*muls).thickness;

    // printf("will create file %s now,nx:%d, ny:%d\n",fileName,nx,ny);
    tcreateFloatPixFile(fileName,image, (long)nx,(long)ny, 1,param );
    // printf("Wrote %s (%g ... %g)\n",fileName,rmin,rmax);
  }
}
*/

/*-------------------- bessi0() ---------------*/
/*
    modified Bessel function I0(x)
    see Abramowitz and Stegun page 379

    x = (double) real arguments

    12-feb-1997 E. Kirkland
 */
 double bessi0( double x )
 {
 	int i;
 	double ax, sum, t;

 	double i0a[] = { 1.0, 3.5156229, 3.0899424, 1.2067492,
		0.2659732, 0.0360768, 0.0045813 };

 	double i0b[] = { 0.39894228, 0.01328592, 0.00225319,
 		-0.00157565, 0.00916281, -0.02057706, 0.02635537,
 		-0.01647633, 0.00392377};

	ax = fabs( x );
	if( ax <= 3.75 ) {
		t = x / 3.75;
		t = t * t;
		sum = i0a[6];
		for( i=5; i>=0; i--) sum = sum*t + i0a[i];
	} else {
		t = 3.75 / ax;
		sum = i0b[8];
		for( i=7; i>=0; i--) sum = sum*t + i0b[i];
		sum = exp( ax ) * sum / sqrt( ax );
	}
	return( sum );

}  /* end bessi0() */

/*-------------------- bessk0() ---------------*/
/*
    modified Bessel function K0(x)
    see Abramowitz and Stegun page 380

    Note: K0(0) is not define and this function
        returns 1E20

    x = (double) real arguments

    this routine calls bessi0() = Bessel function I0(x)

    12-feb-1997 E. Kirkland
 */
 double bessk0( double x )
 {
 	double bessi0(double);

 	int i;
 	double ax, x2, sum;
 	double k0a[] = { -0.57721566, 0.42278420, 0.23069756,
 		 0.03488590, 0.00262698, 0.00010750, 0.00000740};

 	double k0b[] = { 1.25331414, -0.07832358, 0.02189568,
 		 -0.01062446, 0.00587872, -0.00251540, 0.00053208};

	ax = fabs( x );
	if( (ax > 0.0)  && ( ax <=  2.0 ) ){
		x2 = ax/2.0;
		x2 = x2 * x2;
		sum = k0a[6];
		for( i=5; i>=0; i--) sum = sum*x2 + k0a[i];
		sum = -log(ax/2.0) * bessi0(x) + sum;
	} else if( ax > 2.0 ) {
		x2 = 2.0/ax;
		sum = k0b[6];
		for( i=5; i>=0; i--) sum = sum*x2 + k0b[i];
		sum = exp( -ax ) * sum / sqrt( ax );
	} else sum = 1.0e20;
	return ( sum );

}  /* end bessk0() */



/*--------------------- vzatom() -----------------------------------*/
/*
	return the real space projected atomic potential
	in volt-Angstroms for atomic number Z at radius r

	Z = atomic number 1 <= Z <= 103
	radius  = radius in Angstroms (MUST be > 0)
	tdsFlag: indicated wether we calculate thermal diffuse scattering
                 or use debye-waller factors
	scatFlag: indicates whether we use Doyle-Turner, or Weickenmeier and Kohl
	          scattering factors

  assumed global vars:

#define NZMIN	1	= min Z in featom.tab
#define NZMAX	103	= max Z in featom.tab

int feTableRead=0; = flag to remember if the param file has been read
int nl=3, ng=3; = number of Lorenzians and Gaussians
double fparams[][] = fe parameters

  al and ag calculated using physical constants from:
	H. L. Anderson, editor "A Physicist's Desk Reference",
		2nd edition, Amer. Instit. Physics, 1989

  Formula: Kirkland, "Advanced Computing in Electron Microscopy"
			p. 68/69 (5.9) and (5.10)

  in 2D: vz(x,y)  = 4*pi^2*a0*e*sum_i=1^3{a_i*K_0(2*pi*r*sqrt(b_i)}+
				   +2*pi^2*a0*e*sum_i=1^3{c_i/d_i*exp(-pi^2*r^2/d_i)}
  (r^2 = x^2+y^2)
*/

double vzatom( int Z, double radius,int tdsFlag,int scatFlag)
{
   int i,j;
   double r,sum,t;
   double suml,sumg,x;
   static double **f2par = NULL;
   const double pc1 = 150.4121417;


   /* Lorenzian, Gaussian constants */
   const double al=300.8242834, ag=150.4121417;
   const double pi=3.141592654;

   if( (Z<NZMIN) || (Z>NZMAX) ) return( 0.0 );

   /* read in the table from a file if this is the
	first time this is called */
   if( feTableRead == 0 ) ReadfeTable(scatFlag);

   r = fabs( radius );

   if (scatFlag == WEICK_KOHL) {
     /*     if( r < RMIN ) r = RMIN; */ /* avoid singularity at r=0 */
      if( r < 1.0e-10) r = 1.0e-10;  /* avoid singularity at r=0 */
    suml = sumg = 0.0;

     /* Lorenztians  nl = 3 = # of lorentzians*/
     x = 2.0*pi*r;
     for( i=0; i<2*nl; i+=2 )
       suml += fparams[Z][i]* bessk0( x*sqrt(fparams[Z][i+1]) );

     /* Gaussians  ng = 3 = # of gaussians */
     x = pi*r;
     x = x*x;
     for( i=2*nl; i<2*(nl+ng); i+=2 )
       sumg += fparams[Z][i]*exp(-x/fparams[Z][i+1]) / fparams[Z][i+1];

     return( al*suml + ag*sumg );
   }
   else {
     /*********************************************************************
      * determining V2D(r) using Doyle & Turner expansion
      * Scattering factors from :
      * L.-M. Peng Micron 30 (1999) p. 625-648
      *
      * We also need the Debye-Waller factor B=8*pi^2*avg(u^2), which is
      * expected to be stored in the last parameter in fparams
      * ReadFeTable will assign a default value of 0.5 as the Debye-Waller
      * factor, but one may want to change that by setting this last parameter
      * manally before running this routine.
      *
      * All length units are in Angstroem
      *
      * The potential is calculated according to the following formula:
      *
      * V(r) = 8*pi^2*a0*e*sum{a_i/(b_i+B)*exp(-4*pi^2*r^2/(b_i+B))}
      *
      *******************************************************************/
     if (f2par == NULL) {
       /* make calculation more efficient by combining some of the
	  parameters to new ones:
       */
       f2par = double2D( NZMAX+1,NPDTMAX, "f2par" );
       for (i=1;i<=NZMAX;i++) for (j=0;j<NPDTMAX;j+=2) {

	 if (tdsFlag)
	   t = fparams[i][j+1]/4.0; /* t = (b_i)/4 */
	 else
	   t = (fparams[i][j+1]+fparams[i][NPDTMAX])/4.0; /* t = (b_i+B)/4 */

	 f2par[i][j] = pc1*fparams[i][j]/t;
	 f2par[i][j+1] = -pi*pi/t;
	 /*
	   if (i==17)
	   printf("%g %g, t=%g, pi=%g, B=%g\n",f2par[i][j],f2par[i][j+1],t,pi,fparams[i][NPDTMAX]);
	 */
       }
     }

     sum = 0;
     for (i=0;i<NPDTMAX;i+=2)
       sum += f2par[Z][i]*exp(f2par[Z][i+1]*r*r);
     return(sum);
   }
}  /* end vzatom() */

/*********************************************************************/
/*--------------------- v3Datom() -----------------------------------*/
/* written by Christoph Koch 5/19/00

	return the real space 3D atomic potential
	in Volt for atomic number Z at radius r

	Z = atomic number 1 <= Z <= 103
	r = radius in A
	tdsFlag: indicated wether we calculate thermal diffuse scattering
                 or use debye-waller factors
	scatFlag: indicates whether we use Doyle-Turner, or Weickenmeier and Kohl
	          scattering factors

	assumed global vars:

	#define NZMIN	1	= min Z in featom.tab
	#define NZMAX	103	= max Z in featom.tab

	int nl=3, ng=3; = number of Lorenzians and Gaussians
	double fparams[][] = fe parameters

	ag calculated using physical constants from:
	H. L. Anderson, editor "A Physicist's Desk Reference",
	2nd edition, Amer. Instit. Physics, 1989

	Formula: Kirkland, "Advanced Computing in Electron Microscopy"
	p. 68/69 (5.9) and (5.10)

v3(x,y,z)= 2*pi^2*a0*e*sum{a_i/r*exp(-2*pi*r*sqrt(b_i))+
	  +2*pi^(5/2)*a0*e*sum{c_i*d_i^(-3/2)*exp(-pi^2*r^2/d_i)}

    (r^2 = x^2+y^2+z^2)

  simplify the 3D formula to
  v3 = pc1*[sum{p1_i/r*exp(p2_i*r)}+pc2*sum{p3_i*exp(p4_i*r^2)}]

  pc1   = 2*pi^2*a0*e = ag
  pc2   = sqrt(pi)
  p1_i = a_i
  p2_i = -2*pi*sqrt(b_i)
  p3_i = c_i*d_i^(-3/2)
  p4_i = -pi^2*r^2/d_i

*/

double v3Datom(int Z, double r,int tdsFlag,int scatFlag)
{
   int i,j,iz;
   double sum1,t;
   double sum2,r2;
   static double **f2par = NULL;

   /* Gaussian constants */
   const double pi=3.141592654;
   const double pc1 = 150.4121417;
   const double pc2 = sqrt(pi);


   if( (Z<NZMIN) || (Z>NZMAX) ) return( 0.0 );
   iz = Z-1;

   /*********************************************************
    * read in the table from a file if this is the
    * first time this is called
    ********************************************************/
   if(fparams == NULL)
     ReadfeTable(scatFlag);

   if (scatFlag == WEICK_KOHL) {
     if (f2par == NULL) {
       /* make calculation more efficient by combining some of the
	  parameters to new ones:
       */
       f2par = double2D( NZMAX+1,2*(nl+nl), "f2par" );
       for (i=0;i<NZMAX;i++) {
	 for (j=0;j<2*nl;j+=2) {
	   f2par[i][j]      = fparams[i+1][j];                /* p1_i = a_i */
	   f2par[i][j+1]    = -2*pi*sqrt(fparams[i+1][j+1]);  /* p2_i = -2*pi*sqrt(b_i) */
	   t = fparams[i][2*nl+j+1];
	   f2par[i][j+1+2*nl] =  -pi*pi/t;                  /* p4_i = -pi^2*r^2/d_i */
	   t = sqrt(t*t*t);
	   f2par[i][j+2*nl] =  fparams[i+1][2*nl+j]/t;        /* p3_i = c_i*d_i^(-3/2) */
	 }
       }
     }

     /* avoid singularity at r=0 */
     /*     if( r < RMIN )  r = RMIN;*/
     if( r < 1.0e-10 )  r = 1.0e-10;
     sum1 = sum2 = 0.0;

     /*****************************************************************************
      * Formula: v3D = pc1*[1/r*sum{p1_i*exp(p2_i*r)}+pc2*sum{p3_i*exp(p4_i*r^2)}]
      ****************************************************************************/
     /* first sum: */
     for( i=0; i<2*nl; i+=2 )
       sum1 += f2par[iz][i]*exp(f2par[iz][i+1]*r);

     /* second sum: */
     r2 = r*r;
     for( i=2*nl; i<2*(nl+ng); i+=2 )
       sum2 += f2par[Z][i]*exp(f2par[Z][i+1]*r2);

     return(pc1*(sum1/r+pc2*sum2));
   }
   else {  /* if scatFlag==DOYLE_TURNER */
     /*********************************************************************
      * determining V(r) using Doyle & Turner expansion
      * Scattering factors from :
      * L.-M. Peng Micron 30 (1999) p. 625-648
      *
      * We also need the Debye-Waller factor B=8*pi^2*avg(u^2), which is
      * expected to be stored in the last parameter in fparams
      * ReadFeTable will assign a default value of 0.5 as the Debye-Waller
      * factor, but one may want to change that by setting this last parameter
      * manally before running this routine.
      *
      * All length units are in Angstroem
      *
      * The potential is calculated according to the following formula:
      *
      * V(r) = 8*pi^5/2*a0*e*sum{a_i*((b_i+B)/4)^-3/2*exp(-pi^2*r^2*4/(b_i+B))}
      *
      *******************************************************************/
     if (f2par == NULL) {
       /* make calculation more efficient by combining some of the
	  parameters to new ones:
       */
       f2par = double2D( NZMAX+1,NPDTMAX, "f2par" );
       for (i=1;i<=NZMAX;i++) for (j=0;j<NPDTMAX;j+=2) {
	 if (tdsFlag)
	   t=fparams[i][j+1]/4.0;   /* t = (b_i)/4 */
	 else
	   t = (fparams[i][j+1]+fparams[i][NPDTMAX])/4.0; /* t = (b_i+B)/4 */

	 f2par[i][j] = pc1*pc2*fparams[i][j]/sqrt(t*t*t);
	 f2par[i][j+1] = -pi*pi/t;
	 /*
	   if (i==17)
	   printf("%g %g, t=%g, pi=%g, B=%g\n",f2par[i][j],f2par[i][j+1],t,pi,fparams[i][NPDTMAX]);
	 */
       }
     }

     sum1 = 0;
     for (i=0;i<NPDTMAX;i+=2)
       sum1 += f2par[Z][i]*exp(f2par[Z][i+1]*r*r);
     return(sum1);
   }
}  /* end v3Datom() */




/*--------------------- wavelength() -----------------------------------*/
/*
	return the electron wavelength (in Angstroms)
	keep this is one place so I don't have to keep typing in these
	constants (that I can never remember anyhow)

	ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
		(The American Institute of Physics, New York) 1989
		page 4.

	kev = electron energy in keV

*/

double wavelength( double kev )
{
  double w;
  const double emass=510.99906; /* electron rest mass in keV */
  const double hc=12.3984244; /* Planck's const x speed of light*/

  /* electron wavelength in Angstroms */
  w = hc/sqrt( kev * ( 2*emass + kev ) );

  return( w );

}  /* end wavelength() */


double getTime() {
/*
#ifdef WIN32
	return (double)time(NULL);
#else
  timev tv;
  timez tz;

  gettimeofday(&tv, &tz);
  // printf("Time: %ldsec, %ld usec\n",tv.tv_sec,tv.tv_usec);
  */
  return 0;//(double)tv.tv_sec+0.000001*(double)tv.tv_usec;
//#endif
}


/*--------------------- cputim() -----------------------------------*/
/*
  retrieve current CPU time in seconds
*/

double cputim()
{
//  return ( ( (double)clock() ) / ( (double)CLOCKS_PER_SEC) );
return 0;
}  /* end cputim() */



/*-------------------- rangauss() -------------------------- */
/*
  Return a normally distributed random number with
  zero mean and unit variance using Box-Muller method

  ranflat() is the source of uniform deviates

  ref.  Numerical Recipes, 2nd edit. page 289

  added log(0) test 10-jan-1998 E. Kirkland
*/
double rangauss( unsigned long *iseed )
{
  double ranflat( unsigned long* );
  double x1, x2, y;
  static double tpi=2.0*3.141592654;

  /* be careful to avoid taking log(0) */
  do{
    x1 = ranflat( iseed );
    x2 = ranflat( iseed );

  } while ( (x1 < 1.0e-30) || (x2 < 1.0e-30) );

  y = sqrt( - 2.0 * log(x1) ) * cos( tpi * x2 );

  return( y );

}  /* end rangauss() */

/*************************************************************/
/*--------------------- ReadfeTable() -----------------------*/
/*
   read electron scattering factors parameters
	   from file fparam.dat

  the constants that must be defined above are

#define NPMAX	12	= number of parameters for each Z
#define NZMIN	1	= min Z in featom.tab
#define NZMAX	103	= max Z in featom.tab
#define	NCMAX	132	= characters per line to read

  assumed global vars:

int feTableRead=0; = flag to remember if the param file has been read
int nl=3, ng=3; = number of Lorenzians and Gaussians
double fparams[][] = fe parameters

********************************************************************/
int ReadfeTable(int scatFlag)
{
   static char cline[NCMAX],dummy[32];
   static char fileName[32] = "sfact_peng.dat";
   int na,i,j;
   double w;
   char *cstatus;
   int n, zi, z;
   FILE *fpTable;

   /* if the file has been read already then just return */
   if( feTableRead == 1 ) return(0);

   if (scatFlag == DOYLE_TURNER)
     sprintf(fileName,"sfact_peng.dat");
   else
     sprintf(fileName,"fparams.dat");

   if ( (fpTable = fopen( fileName, "r" )) == NULL )  {
	printf("ReadfeTable() can't open file %s\n",fileName);
	exit( 0 );
   }

   if (scatFlag == WEICK_KOHL) {
     printf("Reading Weickenmeier & Kohl parameters from %s\n",fileName);
     na = 2*( nl + ng );	/* number of params to fit */
     fparams = double2D( NZMAX+1, na+1, "fparams" );
     n = 0;

     for( zi=NZMIN; zi<=NZMAX; zi++) {

       /* define some (random) debye waller factors */
       if (zi == 14)
	 fparams[zi][na] = 0.45; /*  for Si */
       else if (zi == 79)
	 fparams[zi][na] = 0.186;  /* for Au */
       else if (zi == 29)
	 fparams[zi][na] = 0.21;  /* for Cu (liquid N2) */
       else
	 fparams[zi][na] = 0.5;  /* for all other elements */

       /* find Z delimiter */
       do {
		   cstatus = fgets( cline, NCMAX, fpTable );
	 if( cstatus == NULL ) break;
       } while ( strncmp( cline, "Z=", 2 ) != 0 );
       if( cstatus == NULL ) break;

       n += 1;
       sscanf( cline, "Z=%d,  chisq=%lf\n", &z, &w);
       for( j=0; j<na; j+=4 ) {
		   fgets( cline, NCMAX, fpTable );
	 for( i=0; i<4; i++) {
	   sscanf(&cline[i*17],"%le", &fparams[z][i+j] );
	 }
       }
       if( z != zi ) {  /* test integrity of the file */
	 printf( "Warning, Z= %d read when expecting"
		 " Z= %d, in ReadfeTable().\n",
		 z, zi);
       }
     }  /* end for(zi=2.. */
   }
   else {  /* DOYLE_TURNER */
     printf("Reading Doyle & Turner parameters from %s\n",fileName);
     fparams = double2D( NZMAX+1, NPDTMAX+1, "fparams" );
     n = 0;
     for( zi=NZMIN; zi<=NZMAX; zi++) {
       /* define some (random) debye waller factors */
       if (zi == 14)
	 fparams[zi][NPDTMAX] = 0.45; /*  for Si */
       else if (zi == 79)
	 fparams[zi][NPDTMAX] = 0.186;  /* for Au */
       else if (zi == 6)
	 fparams[zi][NPDTMAX] = 0.4;  /* for C */
       else
	 fparams[zi][NPDTMAX] = 0.5;  /* for all other elements */
       /*     printf("Using Debye-Waller Factor DW=%f for Z=%d\n",
	      fparams[zi][NPDTMAX],zi);
       */

	   if ( fgets( cline, NCMAX, fpTable ) == NULL )
	 break;
       sscanf(cline,"%s %d %le %le %le %le %le %le %le %le",
	      dummy,&z,&fparams[zi][0],&fparams[zi][2],
	      &fparams[zi][4],&fparams[zi][6],
	      &fparams[zi][1],&fparams[zi][3],&fparams[zi][5],&fparams[zi][7]);
       if (z != zi)
	 printf("Warning: corrupt data file %s!\n",fileName);
       n++;
     }  /* end for(zi=1.. */
   } /* DOYLE_TURNER */

   if( n != (NZMAX-NZMIN + 1) ) {
     printf("Warning, only %d elements read in "
	    "in feTableRead() (too small).\n", n );
   }
   fclose( fpTable );

   feTableRead = 1;	/* remember that table has been read */
   return( n );

} /* end ReadfeTable() */





/*--------------------- sigma() -----------------------------------*/
/*
	return the interaction parameter sigma in radians/(kv-Angstroms)
	keep this in one place so I don't have to keep typing in these
	constants (that I can never remember anyhow)

	ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
		(The American Institute of Physics, New York) 1989
		page 4.

	kev = electron energy in keV

*/

double sigma( double kev )
{
  double s, pi, wavl, x;
  const double emass=510.99906; /* electron rest mass in keV */
  double wavelength( double kev );  /*  get electron wavelength */

  x = ( emass + kev ) / ( 2.0*emass + kev);
  wavl = wavelength( kev );
  pi = 4.0 * atan( 1.0 );

  s = 2.0 * pi * x / (wavl*kev);  // 2*pi*kz*(1+kev/emaxx)/(2*emass+kev)

  return( s );

}  /* end sigma() */


/*------------------ splinh() -----------------------------*/
/*
	fit a quasi-Hermite  cubic spline

	[1] Spline fit as in H.Akima, J. ACM 17(1970)p.589-602
		'A New Method of Interpolation and Smooth
		Curve Fitting Based on Local Procedures'

	[2] H.Akima, Comm. ACM, 15(1972)p.914-918

	E. Kirkland 4-JUL-85
	changed zero test to be a small nonzero number 8-jul-85 ejk
	converted to C 24-jun-1995 ejk

	The inputs are:
		x[n] = array of x values in ascending order, each X(I) must
			be unique
		y[n] = array of y values corresponding to X(N)
		n  = number of data points must be 2 or greater

	The outputs are (with z=x-x(i)):
		b[n] = array of spline coeficients for (x-x[i])
		c[n] = array of spline coeficients for (x-x[i])**2
		d[n] = array of spline coeficients for (x-x[i])**3
		( x[i] <= x <= x[i+1] )
	To interpolate y(x) = yi + bi*z + c*z*z + d*z*z*z

	The coeficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
	interval. NOTE that the last set of coefficients,
	b[n-1], c[n-1], d[n-1] are meaningless.
*/
void splinh( double x[], double y[],
	     double b[], double c[], double d[], int n)
{
#define SMALL 1.0e-25

  int i, nm1, nm4;
  double m1, m2, m3, m4, m5, t1, t2, m54, m43, m32, m21, x43;

  if( n < 4) return;

  /* Do the first end point (special case),
     and get starting values */

  m5 = ( y[3] - y[2] ) / ( x[3] - x[2] );	/* mx = slope at pt x */
  m4 = ( y[2] - y[1] ) / ( x[2] - x[1] );
  m3 = ( y[1] - y[0] ) / ( x[1] - x[0] );

  m2 = m3 + m3 - m4;	/* eq. (9) of reference [1] */
  m1 = m2 + m2 - m3;

  m54 = fabs( m5 - m4);
  m43 = fabs( m4 - m3);
  m32 = fabs( m3 - m2);
  m21 = fabs( m2 - m1);

  if ( (m43+m21) > SMALL )
    t1 = ( m43*m2 + m21*m3 ) / ( m43 + m21 );
  else
    t1 = 0.5 * ( m2 + m3 );

  /*  Do everything up to the last end points */

  nm1 = n-1;
  nm4 = n-4;

  for( i=0; i<nm1; i++) {

    if( (m54+m32) > SMALL )
      t2= (m54*m3 + m32*m4) / (m54 + m32);
    else
      t2 = 0.5* ( m3 + m4 );

    x43 = x[i+1] - x[i];
    b[i] = t1;
    c[i] = ( 3.0*m3 - t1 - t1 - t2 ) /x43;
    d[i] = ( t1 + t2 - m3 - m3 ) / ( x43*x43 );

    m1 = m2;
    m2 = m3;
    m3 = m4;
    m4 = m5;
    if( i < nm4 ) {
      m5 = ( y[i+4] - y[i+3] ) / ( x[i+4] - x[i+3] );
    } else {
      m5 = m4 + m4 - m3;
    }

    m21 = m32;
    m32 = m43;
    m43 = m54;
    m54 = fabs( m5 - m4 );
    t1 = t2;
  }

  return;

} /* end splinh() */

/*----------------------- seval() ----------------------*/
/*
	Interpolate from cubic spline coefficients

	E. Kirkland 4-JUL-85
	modified to do a binary search for efficiency 13-Oct-1994 ejk
	converted to C 26-jun-1995 ejk
	fixed problem on end-of-range 16-July-1995 ejk

	The inputs are:
		x[n] = array of x values in ascending order, each x[i] must
			be unique
		y[n] = array of y values corresponding to x[n]
		b[n] = array of spline coeficients for (x-x[i])
		c[n] = array of spline coeficients for (x-x[i])**2
		d[n] = array of spline coeficients for (x-x[i])**3
		n  = number of data points
		x0  = the x value to interpolate at
		(x[i] <= x <= x[i+1]) and all inputs remain unchanged

	The value returned is the interpolated y value.

	The coeficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
	interval. NOTE that the last set of coefficients,
	b[n-1], c[n-1], d[n-1] are meaningless.
*/
double seval( double *x, double *y, double *b, double *c,
	     double *d, int n, double x0 )
{
  int i, j, k;
  double z, seval1;

  /*  exit if x0 is outside the spline range */
  if( x0 <= x[0] ) i = 0;
  else if( x0 >= x[n-2] ) i = n-2;
  else {
    i = 0;
    j = n;
    do{ k = ( i + j ) / 2 ;
    if( x0 < x[k] )  j = k;
    else if( x0 >= x[k] ) i = k;
    } while ( (j-i) > 1 );
  }

  z = x0 - x[i];
  seval1 = y[i] + ( b[i] + ( c[i] + d[i] *z ) *z) *z;

  return( seval1 );

} /* end seval() */

/*---------------------------- ranflat -------------------------------*/
/*
  return a random number in the range 0.0->1.0
  with uniform distribution

  the 'Magic Numbers' are from
  Numerical Recipes 2nd edit pg. 285
*/
double ranflat( unsigned long *iseed )
{
  static unsigned long a=1366, c=150889L, m=714025L;

  *iseed = ( a * (*iseed) + c ) % m;

  return( ((double) *iseed)/m);

}  /* end ranflat() */

/*--------------------- parlay() -----------------------------------*/
/*
  subroutine to parse the atomic layer stacking sequence
  for use with multislice programs.

  This converts layer structure definition of the form:
      2(abc)d
  into a sequence of numerical indices where a=1, b=2, c=3, ...
  The main attraction is the repeat operator n(...) where
  the structure inside the parenthesis (...) is repeated n times.
  for instance the above structure is translated into:
      1 2 3 1 2 3 4
  The parenthesis may be nested up to 100 levels (determined by
  nlmax). For instance  5( 2(abc) 3(cba) ) is also a valid structure
  definition. This is a compact way of specifying the layer structure
  for a multislice calculation. tab's may be present in the structure
  anywhere (they are ignored).

  This is done by pushing the position of each '(' and its repeat
  count onto a stack. Each ')' pops the last entry from the stack
  and invokes a duplication process.

  Layers refer to the distinquishable subset of different
  types of layers, whereas slices refer to the way these
  layers are sequentially stacked.

  fortran version started 22-dec-1989 earl j. kirkland
  added nested parenthesis by stacking entry points
     31-jan-1990 ejk
  added tab handling (in ascii and ebcdic) 1-feb-1990 ejk
  minor changes to error messages to make rs6000's happy
      7-july-1992 ejk
  converted to ANSI-C 26-jun-1995 ejk
  added include of stdlib.h 6-july-1995 ejk

   c             = input character string
   islice[nsmax] = integer array to get stacking sequence indicies
   nsmax         = (integer) size of layer
   lmax          = (integer) maximum allowed layer index
   nslice        = (integer) number of layers
   returned value= (integer) success/error code
                       0 : success
                      -1 : layer out of range
                      -2 : missing left parenthesis
                      -3 : parenthesis nested too deep
                      -4 : bad repeat code
                      -5 : unmatched right parenthesis
                      -6 : invalid character
                      -7 : too many layers
                      -8 : incomplete stacking sequence

   fperr       = (int) if this is not a NULL then write
                   error messages

  NOTE: islice and nslice may be modified by this routine

  cname determines the mapping of characters into numbers.

*/

#define NSTKMAX 100	/* maximum stack depth (i.e. maximum level
				 of nesting the parenthesis) */
#define NCHARS 52	/* maximum number of character symbols */

int parlay( const char c[], int islice[], int nsmax, int lmax,
	    int *nslice, int fperr )
{
  /* define our own symbol sequence */
  const char cname[] =
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

  int ic, lenc, name, i, j, istack, n,
    ipoint[NSTKMAX], repeat[NSTKMAX];

  /*  initialize misc constants  */

  *nslice = 0;
  istack = 0;
  lenc = (int) strlen( c );  /* length of c */

  for( ic=0; ic<lenc; ic++ ) {	/* loop over all characters in c */

    /*  skip over embedded spaces, tabs and wierd characters */
    while( isspace( c[ic] )  ) ic++;

    /*  if its a character do this  */
    if( isalpha( c[ic] ) ) {
      for(name=0; name<NCHARS; name++)
	if( c[ic] == cname[name] ) break;
      if( name > lmax ) {
	if( fperr != 0 )
	  printf("Layer /%c/ out of range in the stacking sequence"
		 " in subroutine parlay.\n", c[ic] );
	return( -1 );
      } else {
	if( *nslice >= nsmax ) {
	  if( fperr != 0 )
	    printf("Too many layers generated in the stacking"
		   " sequence in parlay.\n");
	  return( -7 );
	}
	islice[ (*nslice)++ ] = name;
      }

      /*  if its a number then extract repeat count up to the '('
	  and save it until a ')' appears (i.e. push on the stack)
      */
    } else if ( isdigit( c[ic] ) ) {
      if( istack >= NSTKMAX ) {
	if( fperr != 0 )
	  printf("Parenthesis nested too deep in "
		 "the stacking sequence in subroutine parlay.\n");
	return( -3 );
      }
      repeat[istack] = atoi( &c[ic] ) - 1;
      while( isdigit(c[ic]) || isspace(c[ic]) ) ic++;
      if( c[ic] != '(' ) {
	if( fperr != 0 )
	  printf("Missing left parenthesis in "
		 "the stacking sequence in subroutine parlay.'\n");
	for(i=0; i<=ic; i++) printf("%c", c[i]);
	printf("?\n");
	return( -2 );
      }
      ipoint[istack++] = *nslice;

      /*  if its a ')' then repeat back to the last '('
	  (i.e. pop the stack) */

    } else if( c[ic] == ')' ) {
      if( istack <= 0 ) {
	if( fperr != 0 )
	  printf("Unmatched right parenthesis in "
		 "the stacking sequence in subroutine parlay.\n");
	return( -5 );
      }
      n = *nslice;
      istack--;
      for( j=0; j<repeat[istack]; j++)
	for(i=ipoint[istack]; i<n; i++){
	  if( *nslice >= nsmax ) {
	    if( fperr != 0 )
	      printf("Too many layers generated in the stacking"
		     " sequence in parlay.\n");
	    return( -7 );
	  }
	  islice[ (*nslice)++ ] = islice[ i ];
	}
    } else {
      if( fperr != 0 )
	printf("Invalid character /%c/ encountered in "
	       "the stacking sequence in subroutine parlay.\n",
	       c[ic]);
      return( -6 );
    }

  } /* end for( ic... */

  if( istack != 0 ) {
    if( fperr != 0 )
      printf("incomplete stacking sequence in parlay.\n");
    return( -8 );
  } else return( 0 );

#undef NSTKMAX
#undef NCHARS

}  /* end parlay() */

/*------------------------- powerof2() ------------------------*/
/*
	return the nearest power of 2 greater than or equal to
	the argument
*/
/*
  long powerof2( long n )
  {
  int ln;
  long n2;

  ln = 1;
  n2 = 2;
  while( (n2 < n) && (ln<31) ) { n2*=2; ln++; }

  return( n2 );
  }
*/



/*********************************************************************/
/*--------------------- fe3D() -----------------------------------*/
/* written by Christoph Koch 5/22/01


	return the reciprocal space 3D atomic potential
	in Volt for atomic number Z at wavevector of length q

	Z = atomic number 1 <= Z <= 103
	q2 = rec. radius in 1/A^2
	tdsFlag: indicates whether we simulate temperature using
	    Debye-Waller factors, or random displacement of atoms
	scale: overall scaling factor, e.g. scaling for FFT
	scatFlag: whether we use doyle-turner or weickenmeier-kohl parameters

	assumed global constants:

	#define NZMIN	1	= min Z in sfact_peng.dat
	#define NZMAX	98	= max Z in sfact_peng.dat

*******************************************************************/

double fe3D(int Z, double q2,int tdsFlag,double scale,int scatFlag)
{
  int i,j; // iz;
  double sum=0.0; // t;
  static double **f2par = NULL;

   /* Gaussian constants */
   const double pi=3.141592654;
   const double a0 = .529;  /* A */
   const double echarge = 14.39;  /* units: V*A */

   if( (Z<NZMIN) || (Z>NZMAX) ) return( 0.0 );

   /*********************************************************
    * read in the table from a file if this is the
    * first time this is called
    ********************************************************/
   if(fparams == NULL)
     ReadfeTable(scatFlag);

   if (scatFlag == WEICK_KOHL) {
     if (f2par == NULL) {
       printf("Will use Weickenmeier & Kohl electron scattering "
	      "factor parameterization\n");
       f2par = fparams;
       if (!tdsFlag)
	 printf("Will use DW-factor [B(Si)=%g]\n",fparams[14][2*(nl+ng)]);
       else
	 printf("Will not use DW-factors\n");

     }
     sum = 0.0;
     for (j=0;j<2*nl;j+=2)
       sum += fparams[Z][j]/(fparams[Z][j+1]+q2);
     for (;j<2*(nl+ng);j+=2)
       sum += fparams[Z][j]*exp(-fparams[Z][j+1]*q2);
     if (!tdsFlag)
       sum *= exp(-fparams[Z][2*(nl+ng)]*q2/(16*pi*pi));
     sum *= 2*pi*a0*echarge;
   }
   else {  /* ifdef DOYLE_TURNER */
     /*********************************************************************
      * determining V(r) using Doyle & Turner expansion
      * Scattering factors from :
      * L.-M. Peng, Micron 30, p. 625-648 (1999)
      *
      * We also need the Debye-Waller factor B=8*pi^2*avg(u^2), which is
      * expected to be stored in the last parameter in fparams
      * ReadFeTable will assign a default value of 0.5 as the Debye-Waller
      * factor, but one may want to change that by setting this last parameter
      * manally before running this routine.
      *
      * q = 4*pi*s
      * fe(s) = K*integral{phi(r)*exp(-i*q*r)}dr^3
      *       = K*integral{phi(r)*exp(-2*pi*i*(2s)*r)}dr^3
      * => fe(q/2) = K*integral{phi(r)*exp(-2*pi*i*q*r)}dr^3
      *            = FT(phi(r))
      * fe3D(s) = sum_{i=1}^NPDTMAX [a_i*exp(-b_i*s^2)]
      * fe3D(q) = sum_{i=1}^NPDTMAX [a_i*exp(-b_i*q^2/(16*pi^2))]
      *
      *
      * for Kirklands scattering factor expansion:
      * Weickenmeier and Kohl, Acta Cryst. A47, p. 590-597 (1991)
      * fe(q) = sum_{i=1}^3 [a_i/(q^2+b_i)]+sum_{i=1}^3 [c_i*exp(-d_i*q^2)]
      *
      *******************************************************************/
     if (f2par == NULL) {
       /* make calculation more efficient by combining some of the
	  parameters to new ones:
       */
       printf("Will use Doyle-Turner electron scattering factors\n");
       if (!tdsFlag)
	 printf("Will use DW-factor [B(Si)=%g]\n",fparams[14][NPDTMAX]);
       else
	 printf("Will not use DW-factors\n");

       f2par = double2D( NZMAX+1,NPDTMAX, "f2par" );
       for (i=1;i<=NZMAX;i++) for (j=0;j<NPDTMAX;j+=2) {
	 if (tdsFlag)  /* t = -(b_i)/(16*pi^2) */
	   f2par[i][j+1] =-fparams[i][j+1]/(16.0*pi*pi);
	 /* f2par[i][j+1] =-fparams[i][j+1];  */
	 else   /* t = -(b_i+B)/(16*pi^2) */
	   f2par[i][j+1] = -(fparams[i][j+1]+fparams[i][NPDTMAX])/(16*pi*pi);
	 /*f2par[i][j+1] = -(fparams[i][j+1]+fparams[i][NPDTMAX]); */

	 f2par[i][j] = scale*fparams[i][j];
	 /*
	   if (i==14)
	   printf("f2par[%d][%d] = %g, f2par[%d][%d] = %g\n",
	   i,j,f2par[i][j],i,j+1,f2par[i][j+1]);
	 */
       }
     }

     sum = 0.0;
     for (j=0;j<NPDTMAX;j+=2)
       sum += f2par[Z][j]*exp(f2par[Z][j+1]*q2);

   }
   return(sum);
}  /* end fe3D() */


double sfLUT(double s,int atKind, MULS *muls)
{
   int i;
   double sf;
   static double *splinx=NULL;
   static double **spliny=NULL;
   static double **splinb=NULL;
   static double **splinc=NULL;
   static double **splind=NULL;
   static int sfSize = 0;
   static int atKinds = 0;
   static double maxK = 0;

   if(splinx == NULL) {
     // splinx = s;
     // spliny = sfC;
     sfSize = muls->sfNk;
     splinx = muls->sfkArray;
     spliny = muls->sfTable;
     atKinds = muls->atomKinds;
     splinb = double2D(atKinds,sfSize, "splinb" );
     splinc = double2D(atKinds,sfSize, "splinc" );
     splind = double2D(atKinds,sfSize, "splind" );
     maxK = splinx[sfSize-1];

     for (i=0;i<atKinds;i++)
       splinh(splinx,spliny[i],splinb[i],splinc[i],splind[i],sfSize);
   }


   /* now that everything is set up find the
      scattering factor by interpolation in the table
   */
   if (s > maxK) return 0.0;
   if (atKind < atKinds) {
     sf = seval(splinx,spliny[atKind],splinb[atKind],splinc[atKind],splind[atKind],sfSize,s);
     if (sf < 0) return 0.0;
     return(sf);
   }
   printf("sfLUT: invalid atom kind (%d) - exit!\n",atKind);
   exit(0);
}  /* end sfLUT() */


/* function bicubic from Matlab toolbox
 * zz = values of function F(z,x) F[iz][ix] = zz[iz*Nx+ix];
 * s = z-coordinate
 * t = x-coordinate
 */
#define FX (ptr[xi-1]*x0 + ptr[xi]*x1 + ptr[xi+1]*x2 + ptr[xi+2]*x)
double bicubic(double **ff,int Nz, int Nx,double z,double x) {
  static double x0,x1,x2,f;
  static double *ptr;
  static int xi,zi;

  // Now interpolate using computationally efficient algorithm.
  // s and t are the x- and y- coordinates of the points we are looking for
  // it is assumed that the function ff is defined at equally spaced intervals
  // 0,1,...N-1 in both directions.

  xi = (int)x;
  zi = (int)z;
  x   = x-(double)xi;
  z   = z-(double)zi;


  x0  = ((2.0-x)*x-1.0)*x;
  x1  = (3.0*x-5.0)*x*x+2.0;
  x2  = ((4.0-3.0*x)*x+1.0)*x;
  x   = (x-1.0)*x*x;


  if (Nz > 2) {
    if ((zi > Nz-3) || (xi > Nx-3)) return 0.0;
    if ((zi < 1) && (xi < 1)) return ff[0][0];
    if (zi < 1) return ff[0][xi];
    if (xi < 1) return ff[zi][0];

    ptr = ff[zi-1];
    f   = FX * (((2.0-z)*z-1)*z);
    ptr = ff[zi];
    f   += FX * ((3.0*z-5.0)*z*z+2.0);
    ptr = ff[zi+1];
    f   += FX * (((4.0-3.0*z)*z+1.0)*z);
    ptr = ff[zi+2];
    f   += FX * ((z-1.0)*z*z);
    f   *= 0.25;
  }
  else {
    ptr = ff[0];
    f = 0.5 * FX;
  }
  return f;
}
#undef FX


int atomCompare(const void *atom1,const void *atom2) {
  /*  return (((*(atom *)atom1).z == (*(atom *)atom2).z) ? 0 :
	  (((*(atom *)atom1).z > (*(atom *)atom2).z) ? 1 : -1));
  */
  /* Use the fact that z is the first element in the atom struct */
  return ((*(real *)atom1 == *(real *)atom2) ? 0 :
	  ((*(real *)atom1 > *(real *)atom2) ? 1 : -1));
}
