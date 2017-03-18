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

#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "stemlib.h"
#include "memory_fftw3.h"	/* memory allocation routines */
#include "stemutil.h"
// #include "tiffsubs.h"
#include "imagelib_fftw3.h"
#include "fileio_fftw3.h"
// #include "floatdef.h"
// #include "imagelib.h"


#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

#define OVERSAMP_X 2
#define OVERSAMP_Z 18

#define NSMAX 1000	/* max number of slices */
#define NLMAX	52	/* maximum number of layers */
#define NCINMAX  500	/* max number of characers in stacking spec */
#define NCMAX 256	/* max characters in file names */
#define NPARAM	64	/* number of parameters */
#define NZMIN	1	/* min Z in featom.tab */
#define NZMAX	103	/* max Z in featom.tab */
#define EXTRA_LAYERS 3  /* number of extra layers for potential overlap */

#define SUB_SLICES  5    /* number of sub slices per slice (for integration) */
#define POTENTIAL_3D
#define INTEGRAL_TOL 1e-5
#define MAX_INTEGRAL_STEPS 15
#define MIN_INTEGRAL_STEPS 2
#define OVERSAMPLING 3
#define OVERSAMPLINGZ (3*OVERSAMPLING)
/*#define USE_VATOM_LUT */ /* set if you want to use vzatomLUT/v3DatomLUT */
/*#define USE_VZATOM_IN_CENTER */
/////////////////////////////////////////////////
// for debugging:
#define SHOW_SINGLE_POTENTIAL 0
/////////////////////////////////////////////////



#define BUF_LEN 256
#define PI 3.14159265358979
#define USE_REZ_SFACTS    1  // used in getAtomPotential3D and getAtomPotentialOffset3D
#define Z_INTERPOLATION   0  // used in make3DSlices (central function for producing atom potential slices)
#define USE_Q_POT_OFFSETS 1  // used in make3DSlices (central function for producing atom potential slices)


const char cname[] = "abcdefghijklmnopqrstuvwxyz"
"ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const double pid=PI;
const double pid2=2*PI;

// define all the scat facts (or just those for Sr, Ti, O, In, P, He,Cl, Si,Ca,Ba,Fe) from Rez et al / Doyle and Turner
// the final 3 values must be zero to achieve a nice tapering off
#if USE_REZ_SFACTS
// provide array variable names
// - scatPar[N_ELEM][N_SF] and
// - scatParOffs[N_ELEM][N_SF]
// and also define N_SF and N_ELEM:
#include "scatfactsRez.h"
#else

#define N_SF 30
#define N_ELEM 14
// tabulated scattering factors from Doyle and Turner (copied by hand!)
double scatPar[N_ELEM][N_SF] = {{0.0000,0.0500,0.1000,0.1500,0.2000,0.2500,0.3000,0.3500,0.4000,0.4500,0.5000,
0.6000,0.7000,0.8000,0.9000,1.0000,1.2000,1.4000,1.6000,1.8000,2.0000,2.5000,
3.0000,3.5000,4.0000,5.0000,6.0000,80,90,100},
{13.1090,11.4760,8.4780,6.2000,4.7940,3.8820,3.2240,2.7180,2.3150,1.9930,1.7330,
1.3500,1.0890,0.9020,0.7620,0.6510,0.4880,0.3750,0.2970,0.2390,0.1940,0.1290,
0.0920,0.0690,0.0540,0.0350,0.0240,0,0,0},
{8.7760,7.9370,6.1990,4.6430,3.5640,2.8440,2.3410,1.9640,1.6680,1.4280,1.2300,
0.9300,0.7210,0.5730,0.4670,0.3890,0.2850,0.2190,0.1750,0.1430,0.1170,0.0780,
0.0550,0.0410,0.0310,0.0200,0.0140,0,0,0},
{1.9830,1.9370,1.8080,1.6250,1.4220,1.2220,1.0400,0.8810,0.7470,0.6350,0.5420,
0.4030,0.3070,0.2410,0.1930,0.1590,0.1130,0.0850,0.0660,0.0530,0.0440,0.0290,
0.0200,0.0150,0.0120,0.0080,0.0050,0,0,0},
{10.434,9.7680,8.2970,6.8050,5.6010,4.6820,3.9730,3.4120,2.9550,2.5760,2.2570,
1.7580,1.3970,1.1320,0.9360,0.7890,0.5870,0.4580,0.3680,0.3020,0.2500,0.1660,
0.1180,0.0880,0.0680,0.0345,0.0310,0,0,0},
{5.4880,5.1920,4.4570,3.5860,2.7960,2.1690,1.7020,1.3620,1.1150,0.9330,0.7970,
0.6100,0.4870,0.4010,0.3350,0.2840,0.2100,0.1600,0.1250,0.1000,0.0820,0.0530,
0.0370,0.0280,0.0210,0.0140,0.0100,0,0,0},
{0.4180,0.4100,0.3900,0.3590,0.3230,0.2860,0.2500,0.2170,0.1890,0.1640,0.1430,
0.1100,0.0860,0.0680,0.0550,0.0460,0.0320,0.0240,0.0190,0.0150,0.0120,0.0080,
0.0050,0.0040,0.0030,0.0020,0.0010,0,0,0},
{4.8570,4.6850,4.2270,3.6200,2.9970,2.4380,1.9740,1.6060,1.3190,1.0980,0.8280,
0.6920,0.5410,0.4400,0.3660,0.3110,0.2320,0.1780,0.1410,0.1130,0.0930,0.0600,
0.0420,0.0310,0.0240,0.0160,0.0110,0,0,0},
{5.8280,5.4210,4.4670,3.4370,2.5890,1.9690,1.5340,1.2310,1.0170,0.8610,0.7430,
0.5780,0.4650,0.3830,0.3200,0.2700,0.1980,0.1500,0.1170,0.0930,0.0760,0.0500,
0.0350,0.0260,0.0200,0.0130,0.0090,0,0,0},
{9.9130,8.7030,6.3880,4.5500,3.4080,2.6950,2.2060,1.8380,1.5480,1.3140,1.1230,
0.8380,0.6470,0.5150,0.4220,0.3540,0.2620,0.2020,0.1620,0.1320,0.1070,0.0710,
0.0500,0.0370,0.0280,0.0180,0.0130,0,0,0},
{18.267,15.854,11.675,8.6820,6.8290,5.5700,4.6280,3.8950,3.3180,2.8610,2.4940,
1.9510,1.5700,1.2880,1.0730,0.9040,0.6660,0.5110,0.4110,0.3370,0.2770,0.1890,
0.1340,0.1000,0.0780,0.0510,0.0360,0,0,0},
{7.1650,6.6690,5.5580,4.4360,3.5620,2.9280,2.4610,2.1040,1.8180,1.5840,1.3880,
1.0800,0.8540,0.6860,0.5610,0.4660,0.3360,0.2550,0.2020,0.1650,0.1360,0.0910,
0.0650,0.0480,0.0370,0.0240,0.0170,0,0,0},
// Al:
{5.8990,5.3710,4.2370,3.1280,2.2990,1.7370,1.3630,1.1110,0.9320,0.8010,0.7000,
0.5510,0.4450,0.3660,0.3040,0.2550,0.1850,0.1390,0.1090,0.0870,0.0700,0.0460,
0.0320,0.0240,0.0190,0.0120,0.0090,0,0,0},
// Y:
{12.307,10.968,8.3983,6.3131,4.9361,4.0087,3.3336,2.8186,2.4107,2.0840,1.8121,
1.4173,1.1500,0.9536,0.8049,0.6849,0.5114,0.3917,0.3098,0.2480,0.2031,0.1363,
0.0957,0.0727,0.0569,0.0369,0.0258,0,0,0}};
#endif  // USE_REZ_SFACTS
/****************************************************************************
* function: atomBoxLookUp
*
* Znum = element
* x,y,z = real space position (in A)
* B = Debye-Waller factor, B=8 pi^2 <u^2>
***************************************************************************/
void atomBoxLookUp(fftw_complex *vlu,MULS *muls,int Znum,double x,double y,double z,double B) {
	static int boxNx,boxNy,boxNz;
	static double dx,dy,dz,ddx,ddy,ddz;
	static atomBox *aBox = NULL;
	static int ix,iy,iz; // idz, intSteps;
	// static double x2,y2,z2,r2;
	// static int avgSteps,maxSteps,stepCount,maxStepCount;
	static double maxRadius2;
	static char fileName[256],systStr[256];
	static fftw_complex sum;
	static int tZ, tnx, tny, tnz, tzOversample;
	static double tdx, tdy, tdz, tv0, tB;
	FILE *fpBox;
	int numRead = 0,dummy;


	(*vlu)[0] = 0.0;
	(*vlu)[1] = 0.0;

	/* initialize all the atoms to non-used */
	if (aBox == NULL) {
		aBox = (atomBox *)malloc(sizeof(atomBox)*(NZMAX+1));
		for (ix=0;ix<=NZMAX;ix++) {
			aBox[ix].potential = NULL;
			aBox[ix].rpotential = NULL;
			aBox[ix].B = -1.0;
		}
		dx = (*muls).resolutionX;    ddx = dx/(double)OVERSAMPLING;
		dy = (*muls).resolutionY;    ddy = dy/(double)OVERSAMPLING;
		dz = (*muls).sliceThickness; ddz = dz/(double)OVERSAMPLINGZ;
		maxRadius2 = (*muls).atomRadius*(*muls).atomRadius;
		/* For now we don't care, if the box has only small
		* prime factors, because we will not fourier transform it
		* especially not very often.
		*/
		boxNx = (int)((*muls).atomRadius/ddx+2.0);
		boxNy = (int)((*muls).atomRadius/ddy+2.0);
		boxNz = (int)((*muls).atomRadius/ddz+2.0);
		if ((*muls).potential3D == 0)
			boxNz = 1;

		if (muls->printLevel > 2)
			printf("Atombox has real space resolution of %g x %g x %gA (%d x %d x %d pixels)\n",
			ddx,ddy,ddz,boxNx,boxNy,boxNz);
	}
	// printf("Debugging: %d %g %g: %g\n",Znum,aBox[Znum].B,B,fabs(aBox[Znum].B - B));

	/* Creating/Reading a atombox for every new kind of atom, but only as needed */
	if (fabs(aBox[Znum].B - B) > 1e-6) {
		//  printf("Debugging 1 (%d: %.7g-%.7g= %.7g), %d\n",
		//	   Znum,aBox[Znum].B,B,fabs(aBox[Znum].B - B),fabs(aBox[Znum].B - B) > 1e-6);
		aBox[Znum].B = B;
		/* Open the file with the projected potential for this particular element
		*/
		sprintf(fileName,"potential_%d_B%d.prj",Znum,(int)(100.0*B));
		if ( (fpBox = fopen( fileName, "r" )) == NULL ) {
			sprintf(systStr,"scatpot %s %d %g %d %d %d %g %g %g %d %g",
				fileName,Znum,B,boxNx,boxNy,boxNz,ddx,ddy,ddz,OVERSAMPLINGZ,(*muls).v0);
			if (muls->printLevel > 2) {
				printf("Could not find precalculated potential for Z=%d,"
					" will calculate now.\n",Znum);
				printf("Calling: %s\n",systStr);
			}
			system(systStr);
			for (dummy=0;dummy < 10000;dummy++);
			if ( (fpBox = fopen( fileName, "r" )) == NULL ) {
				if (muls->printLevel >0)
					printf("cannot calculate projected potential using scatpot - exit!\n");
				exit(0);
			}

		}
		fgets( systStr, 250, fpBox );
		sscanf(systStr,"%d %le %d %d %d %le %le %le %d %le\n",
			&tZ, &tB, &tnx, &tny, &tnz, &tdx, &tdy, &tdz, &tzOversample, &tv0);
		/* If the parameters in the file don't match the current ones,
		* we need to create a new potential file
		*/
		if ((tZ != Znum) || (fabs(tB-B)>1e-6) || (tnx != boxNx) || (tny != boxNy) || (tnz != boxNz) ||
			(fabs(tdx-ddx) > 1e-5) || (fabs(tdy-ddy) > 1e-5) || (fabs(tdz-ddz) > 1e-5) ||
			(tzOversample != OVERSAMPLINGZ) || (tv0 != muls->v0)) {
				if (muls->printLevel > 2) {
					printf("Potential input file %s has the wrong parameters\n",fileName);
					printf("Parameters:\n"
						"file:    Z=%d, B=%.3f A^2 (%d, %d, %d) (%.7f, %.7f %.7f) nsz=%d V=%g\n"
						"program: Z=%d, B=%.3f A^2 (%d, %d, %d) (%.7f, %.7f %.7f) nsz=%d V=%g\n"
						"will create new potential file, please wait ...\n",
						tZ,tB,tnx,tny,tnz,tdx,tdy,tdz,tzOversample,tv0,
						Znum,B,boxNx,boxNy,boxNz,ddx,ddy,ddz,OVERSAMPLINGZ,(*muls).v0);
					/* printf("%d %d %d %d %d %d %d %d %d %d\n",
					(tZ != Znum),(tB != B),(tnx != boxNx),(tny != boxNy),(tnz != boxNz),
					(fabs(tdx-ddx) > 1e-5),(fabs(tdy-ddy) > 1e-5),(fabs(tdz-ddz) > 1e-5),
					(tzOversample != OVERSAMPLINGZ),(tv0 != muls->v0));
					*/
				}
				/* Close the old file, Create a new potential file now
				*/
				fclose( fpBox );
				sprintf(systStr,"scatpot %s %d %g %d %d %d %g %g %g %d %g",
					fileName,Znum,B,boxNx,boxNy,boxNz,ddx,ddy,ddz,OVERSAMPLINGZ,(*muls).v0);
				system(systStr);
				if ( (fpBox = fopen( fileName, "r" )) == NULL ) {
					if (muls->printLevel >0)
						printf("cannot calculate projected potential using scatpot - exit!\n");
					exit(0);
				}
				fgets( systStr, 250, fpBox );
		}

		/* Finally we can read in the projected potential
		*/
		if (B == 0) {
			aBox[Znum].rpotential = float3D(boxNz,boxNx,boxNy,"atomBox");
			numRead = fread(aBox[Znum].rpotential[0][0],sizeof(real),
				(size_t)(boxNx*boxNy*boxNz), fpBox );
		}
		else {
#if FLOAT_PRECISION == 1
			aBox[Znum].potential = complex3Df(boxNz,boxNx,boxNy,"atomBox");
			numRead = fread(aBox[Znum].potential[0][0],sizeof(fftwf_complex),
				(size_t)(boxNx*boxNy*boxNz), fpBox );
#else
			aBox[Znum].potential = complex3D(boxNz,boxNx,boxNy,"atomBox");
			numRead = fread(aBox[Znum].potential[0][0],sizeof(fftw_complex),
				(size_t)(boxNx*boxNy*boxNz),fpBox);
#endif
		}

		/* writeImage_old(aBox[Znum].potential[0],boxNx,boxNy, 0.0,"potential.img");
		system("showimage potential.img");
		*/
		fclose( fpBox );

		if (numRead == boxNx*boxNy*boxNz) {
			if (muls->printLevel > 1)
				printf("Sucessfully read in the projected potential\n");
		}
		else {
			if (muls->printLevel > 0)
				printf("error while reading potential file %s: read %d of %d values\n",
				fileName,numRead,boxNx*boxNy*boxNz);
			exit(0);
		}
	}

	/***************************************************************
	* Do the trilinear interpolation
	*/
	sum[0] = 0.0;
	sum[1] = 0.0;
	if (x*x+y*y+z*z > maxRadius2) {
		return;
	}
	x = fabs(x);
	y = fabs(y);
	z = fabs(z);
	ix = (int)(x/ddx);
	iy = (int)(y/ddy);
	iz = (int)(z/ddz);
	dx = x-(double)ix*ddx;
	dy = y-(double)iy*ddy;
	dz = z-(double)iz*ddz;
	if ((dx < 0) || (dy<0) || (dz<0)) {
		/* printf("Warning, dx(%g), dy(%g), dz(%g) < 0, (x=%g, y=%g, z=%g)\n",dx,dy,dz,x,y,z);
		*/
		if (dx < 0) dx = 0.0;
		if (dy < 0) dy = 0.0;
		if (dz < 0) dz = 0.0;
	}


	if ((*muls).potential3D) {
		if (aBox[Znum].B > 0) {
			sum[0] = (1.0-dz)*((1.0-dy)*((1.0-dx)*aBox[Znum].potential[iz][ix][iy][0]+
				dx*aBox[Znum].potential[iz][ix+1][iy][0])+
				dy*((1.0-dx)*aBox[Znum].potential[iz][ix][iy+1][0]+
				dx*aBox[Znum].potential[iz][ix+1][iy+1][0]))+
				dz*((1.0-dy)*((1.0-dx)*aBox[Znum].potential[iz+1][ix][iy][0]+
				dx*aBox[Znum].potential[iz+1][ix+1][iy][0])+
				dy*((1.0-dx)*aBox[Znum].potential[iz+1][ix][iy+1][0]+
				dx*aBox[Znum].potential[iz+1][ix+1][iy+1][0]));
			sum[1] = (1.0-dz)*((1.0-dy)*((1.0-dx)*aBox[Znum].potential[iz][ix][iy][1]+
				dx*aBox[Znum].potential[iz][ix+1][iy][1])+
				dy*((1.0-dx)*aBox[Znum].potential[iz][ix][iy+1][1]+
				dx*aBox[Znum].potential[iz][ix+1][iy+1][1]))+
				dz*((1.0-dy)*((1.0-dx)*aBox[Znum].potential[iz+1][ix][iy][1]+
				dx*aBox[Znum].potential[iz+1][ix+1][iy][1])+
				dy*((1.0-dx)*aBox[Znum].potential[iz+1][ix][iy+1][1]+
				dx*aBox[Znum].potential[iz+1][ix+1][iy+1][1]));
		}
		else {
			sum[0] = (1.0-dz)*((1.0-dy)*((1.0-dx)*aBox[Znum].rpotential[iz][ix][iy]+
				dx*aBox[Znum].rpotential[iz][ix+1][iy])+
				dy*((1.0-dx)*aBox[Znum].rpotential[iz][ix][iy+1]+
				dx*aBox[Znum].rpotential[iz][ix+1][iy+1]))+
				dz*((1.0-dy)*((1.0-dx)*aBox[Znum].rpotential[iz+1][ix][iy]+
				dx*aBox[Znum].rpotential[iz+1][ix+1][iy])+
				dy*((1.0-dx)*aBox[Znum].rpotential[iz+1][ix][iy+1]+
				dx*aBox[Znum].rpotential[iz+1][ix+1][iy+1]));
		}
	}
	else {
		if (aBox[Znum].B > 0) {
			sum[0] = (1.0-dy)*((1.0-dx)*aBox[Znum].potential[0][ix][iy][0]+
				dx*aBox[Znum].potential[0][ix+1][iy][0])+
				dy*((1.0-dx)*aBox[Znum].potential[0][ix][iy+1][0]+
				dx*aBox[Znum].potential[0][ix+1][iy+1][0]);
			sum[1] = (1.0-dy)*((1.0-dx)*aBox[Znum].potential[0][ix][iy][1]+
				dx*aBox[Znum].potential[0][ix+1][iy][1])+
				dy*((1.0-dx)*aBox[Znum].potential[0][ix][iy+1][1]+
				dx*aBox[Znum].potential[0][ix+1][iy+1][1]);
		}
		else {
			sum[0] = (1.0-dy)*((1.0-dx)*aBox[Znum].rpotential[0][ix][iy]+
				dx*aBox[Znum].rpotential[0][ix+1][iy])+
				dy*((1.0-dx)*aBox[Znum].rpotential[0][ix][iy+1]+
				dx*aBox[Znum].rpotential[0][ix+1][iy+1]);
		}
	}
	(*vlu)[0] = sum[0];
	(*vlu)[1] = sum[1];
}



/*****************************************************
* void make3DSlices()
*
* This function will create a 3D potential from whole
* unit cell, slice it, and make transr/i and propr/i
* Call this function with center = NULL, if you don't
* want the array to be shifted.
****************************************************/
void make3DSlices(MULS *muls,int nlayer,char *fileIn,atom *center) {
	// FILE *fpu2;
	char fileOut[512]; // RAM: this is terrible, why is fileName a function argument and here we have filename?  FIXED: rename function argument to fileIn and this to fileOut
	int natom,iatom,iz;  /* number of atoms */
	atom *atoms;
	real dx,dy,dz;
	real c,atomX,atomY,atomZ;
	int i=0,j,nx,ny,ix,iy,iax,iay,iaz,sliceStep;
	int iAtomX,iAtomY,iAtomZ,iRadX,iRadY,iRadZ,iRad2;
	int iax0,iax1,iay0,iay1,iaz0,iaz1,nyAtBox,nyAtBox2,nxyAtBox,nxyAtBox2,iOffsX,iOffsY,iOffsZ;
	int nzSub,Nr,ir,Nz_lut;
	int iOffsLimHi,iOffsLimLo,iOffsStep;

	real *slicePos;
	double z,x,y,r,ddx,ddy,ddr,dr,r2sqr,x2,y2,potVal,dOffsZ;
	// char *sliceFile = "slices.dat";
	char buf[BUF_LEN];
	FILE *sliceFp;
	real minX,maxX,minY,maxY,minZ,maxZ;
	double atomRadius2;
	time_t time0,time1;
	float s11,s12,s21,s22;
	fftwf_complex	*atPotPtr;
	float *potPtr=NULL, *ptr;
	static int divCount = 0;
	static real **tempPot = NULL;
	/*
#if FLOAT_PRECISION == 1
	static fftwf_complex ***oldTrans = NULL;
	static fftwf_complex ***oldTrans0 = NULL;
#else
	static fftw_complex ***oldTrans = NULL;
	static fftw_complex ***oldTrans0 = NULL;
#endif
*/
	ImageIOPtr imageIO = ImageIOPtr(new CImageIO(muls->potNx,muls->potNy,
				muls->sliceThickness,muls->resolutionX,muls->resolutionY));
	fftw_complex dPot;
#if Z_INTERPOLATION
	double ddz;
#endif
#if USE_Q_POT_OFFSETS
	fftwf_complex	*atPotOffsPtr;
#endif

	if (muls->trans == NULL) {
		printf("Severe error: trans-array not allocated - exit!\n");
		exit(0);
	}

/*
	if (oldTrans0 == NULL) {
#if FLOAT_PRECISION == 1
		oldTrans0 = (fftwf_complex ***)fftw_malloc(nlayer * sizeof(fftwf_complex**));
#else
		oldTrans0 = (fftw_complex ***)fftw_malloc(nlayer * sizeof(fftw_complex**));
#endif
		for (i=0;i<nlayer;i++) {
			// printf("%d %d\n",i,(int)(muls->trans));
			oldTrans0[i] = muls->trans[i];
		}
		oldTrans = muls->trans;
	}
	if (oldTrans != muls->trans)
		printf("Warning: Transmission function pointer has changed!\n");
*/
	/* return, if there is nothing to do */
	if (nlayer <1)
		return;

	nx = muls->potNx;
	ny = muls->potNy;



	/* we need to keep track of which subdivision of the unit cell we are in
	* If the cell is not subdivided, then muls.cellDiv-1 = 0.
	*/
	if ((divCount == 0) || (muls->equalDivs))
		divCount = muls->cellDiv;
	divCount--;

	/* we only want to reread and shake the atoms, if we have finished the
	* current unit cell.
	*/
	if (divCount == muls->cellDiv-1) {
		if (muls->avgCount == 0) {
			// if this is the first run, the atoms have already been
			// read during initialization
			natom = (*muls).natom;
			atoms = (*muls).atoms;
		}
		else {
			/*
			the following function makes an array of natom atoms from
			the input file (with x,y,z,dw,occ);
			*/
			// the last parameter is handleVacancies.  If it is set to 1 vacancies

			// and multiple occupancies will be handled.
			atoms = readUnitCell(&natom,fileIn,muls,1);
			if (muls->printLevel>=3)
				printf("Read %d atoms from %s, tds: %d\n",natom,fileIn,muls->tds);
			muls->natom = natom;
			muls->atoms = atoms;
		}
		minX = maxX = atoms[0].x;
		minY = maxY = atoms[0].y;
		minZ = maxZ = atoms[0].z;

		for (i=0;i<natom;i++) {
			if (atoms[i].x < minX) minX = atoms[i].x;
			if (atoms[i].x > maxX) maxX = atoms[i].x;
			if (atoms[i].y < minY) minY = atoms[i].y;
			if (atoms[i].y > maxY) maxY = atoms[i].y;
			if (atoms[i].z < minZ) minZ = atoms[i].z;
			if (atoms[i].z > maxZ) maxZ = atoms[i].z;
		}
		/*
		printf("Root of mean square TDS displacement: %f A (wobble=%g at %gK) %g %g %g\n",
		sqrt(u2/natom),wobble,(*muls).tds_temp,ux,uy,uz);
		*/
		if (muls->printLevel >= 2) {
			printf("range of thermally displaced atoms (%d atoms): \n",natom);
			printf("X: %g .. %g\n",minX,maxX);
			printf("Y: %g .. %g\n",minY,maxY);
			printf("Z: %g .. %g\n",minZ,maxZ);
		}
		/*
		define the center of our unit cell by moving the atom specified
		by "center" at position (0.5,0.5,0.0)
		*/

		if (center != NULL) {
			dx = (*muls).ax/2.0f - (*center).x;
			dy = (*muls).by/2.0f - (*center).y;
			dz = -(*center).z;
			for (i=0;i<natom;i++) {
				atoms[i].x += dx;
				if (atoms[i].x < 0.0f) atoms[i].x += (*muls).ax;
				else if (atoms[i].x > (*muls).ax) atoms[i].x -= (*muls).ax;
				atoms[i].y += dy;
				if (atoms[i].y < 0.0f) atoms[i].y += (*muls).by;
				else if (atoms[i].y > (*muls).by) atoms[i].y -= (*muls).by;
				atoms[i].z += dz;
				if (atoms[i].z < 0.0f) atoms[i].z += (*muls).c;
				else if (atoms[i].z > (*muls).c) atoms[i].z -= (*muls).c;
			}
		}

		/**********************************************************
		* Sort the atoms in z.
		*********************************************************/
		// RAM I think this is not working correctly to output the file
		// Muls passed in as *muls, so by pointer rather than by-value.
		// Ok, and cfgFile is not set before here...  Look at Muls and see if there's been some variable confusion


		//printf( "DEBUG: stemlib::make3Dslices : muls.cfgFile = %s \n", muls->cfgFile );

		qsort(atoms,natom,sizeof(atom),atomCompare);

		/*
		if ((*muls).cfgFile != NULL)
		{
			sprintf(buf,"%s/%s",muls->folder,muls->cfgFile);
			// append the TDS run number
			if (strcmp(buf+strlen(buf)-4,".cfg") == 0) *(buf+strlen(buf)-4) = '\0';
			if (muls->tds) sprintf(buf+strlen(buf),"_%d.cfg",muls->avgCount);
			else sprintf(buf+strlen(buf),".cfg");

			// printf("Will write CFG file <%s> (%d)\n",buf,muls->tds)
			writeCFG(atoms,natom,buf,muls);

			if (muls->readPotential)
			{
				sprintf(buf,"nanopot %s/%s %d %d %d %s",muls->folder,muls->cfgFile,
					ny,nx,muls->slices*muls->cellDiv,muls->folder);
				system(buf);
			}
		}
		*/
	} /* end of if divCount==cellDiv-1 ... */
	else {
		natom = muls->natom;
		atoms = muls->atoms;
	}
	/**************************************************************
	*	setup the slices with their start and end positions
	*	then loop through all the atoms and add their potential to
	*	the slice that their potential reaches into (up to RMAX)
	*************************************************************/
	// c = (*muls).c/(real)((*muls).cellDiv);
	c = muls->sliceThickness * muls->slices;
	dx = (*muls).resolutionX;
	dy = (*muls).resolutionY;
	dr   = muls->resolutionX/OVERSAMP_X;  // define step width in which radial V(r,z) is defined
	iRadX = (int)ceil((*muls).atomRadius/dx);
	iRadY = (int)ceil((*muls).atomRadius/dy);
	iRadZ = (int)ceil((*muls).atomRadius/muls->sliceThickness);
	iRad2 = iRadX*iRadX+iRadY*iRadY;
	atomRadius2 = (*muls).atomRadius * (*muls).atomRadius;

	if (muls->printLevel >= 3) {
		printf("Slab thickness: %gA z-offset: %gA (cellDiv=%d)\n",
			c,c*(real)(muls->cellDiv-divCount-1),divCount);
	}
	/*******************************************************
	* initializing slicPos, cz, and transr
	*************************************************************/
	/*
	if (oldTrans0[0] != muls->trans[0]) {
		printf("Warning: transmision array pointers have changed!\n");
		for (i=0;i<nlayer;i++)
			muls->trans[i] = oldTrans0[i];
	}
	*/
	/*
	if (oldTrans0[0][0] != muls->trans[0][0]) {
	printf("Warning: transmision array pointers have changed!\n");
	for (i=0;i<nlayer;i++)
	muls->trans[i] = oldTrans0[i];
	}
	*/
	if ((*muls).cz == NULL) {
		(*muls).cz = float1D(nlayer,"cz");
	}
	// sliceFp = fopen(sliceFile,"r");
	sliceFp = NULL;
	slicePos = float1D(nlayer,"slicePos");


	if (muls->sliceThickness == 0)
		(*muls).cz[0] = c/(real)nlayer;
	else
		(*muls).cz[0] = muls->sliceThickness;
	slicePos[0] = (*muls).czOffset;
	/*
	************************************************************/


	for (i=1;i<nlayer;i++) {
		if (sliceFp == NULL) (*muls).cz[i] = (*muls).cz[0];
		/* don't need to all be the same, yes they do for fast 3D-FFT method! */
		else {
			fgets(buf,BUF_LEN,sliceFp);
			(*muls).cz[i] = atof(buf);
		}
		slicePos[i] = slicePos[i-1]+(*muls).cz[i-1]/2.0+(*muls).cz[i]/2.0;
	}

	memset(muls->trans[0][0],0,nlayer*nx*ny*sizeof(fftwf_complex));
	/* check whether we have constant slice thickness */

	if (muls->fftpotential) {
		for (i = 0;i<nlayer;i++)  if ((*muls).cz[0] != (*muls).cz[i]) break;
		if (i<nlayer) printf("Warning: slice thickness not constant, will give wrong results (iz=%d)!\n",i);

	}

	/*************************************************************************
	* read the potential that has been created externally!
	*/
	if (muls->readPotential) {
		for (i=(divCount+1)*muls->slices-1,j=0;i>=(divCount)*muls->slices;i--,j++) {
			sprintf(buf,"%s/potential_%d.img",muls->folder,i);
			imageIO->ReadImage((void **)tempPot,nx,ny,buf);
			for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
				(*muls).trans[j][ix][iy][0] = tempPot[ix][iy];
				(*muls).trans[j][ix][iy][1] = 0.0;
			}
		}
		return;
	}

	// reset the potential to zero:
#if FLOAT_PRECISION == 1
	memset((void *)&(muls->trans[0][0][0][0]),0,
		muls->slices*muls->potNx*muls->potNy*sizeof(fftwf_complex));
#else
	memset((void *)&(muls->trans[0][0][0][0]),0,
		muls->slices*muls->potNx*muls->potNy*sizeof(fftw_complex));
#endif
	nyAtBox   = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY);
	nxyAtBox  = nyAtBox*(2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX));
	nyAtBox2  = 2*nyAtBox;
	nxyAtBox2 = 2*nxyAtBox;
	sliceStep = 2*muls->potNx*muls->potNy;

	/*
	for (i=0;i<nlayer;i++)
	printf("slice center: %g width: %g\n",slicePos[i],(*muls).cz[i]);
	*/

	/****************************************************************
	* Loop through all the atoms and add their potential
	* to the slices:
	***************************************************************/

	time(&time0);
	for (iatom = 0;iatom<natom;iatom++) {
		// make sure we skip vacancies:
		while (atoms[iatom].Znum == 0) iatom++;
		if (iatom >=natom) break;

		if ((muls->printLevel >= 4) && (muls->displayPotCalcInterval > 0)) {
			if (((iatom+1) % (muls->displayPotCalcInterval)) == 0) {
				printf("Adding potential for atom %d (Z=%d, pos=[%.1f, %.1f, %.1f])\n",iatom+1,atoms[iatom].Znum,atoms[iatom].x,atoms[iatom].y,atoms[iatom].z);
			}
		}
		// printf("c=%g, slice thickness=%g, slices=%d, %d\n",c,muls->sliceThickness,muls->slices,muls->displayPotCalcInterval);
		/*
		* c = the thickness of the current slab.
		*
		* if the z-position of this atom is outside the potential slab
		* we won't consider it and skip to the next
		*/
		/* cellDiv = number of times that the big super cell is divided into
		* less big ones, yet often still bigger than a single unit cell
		* (for saving memory)
		* divCount = counter of how many such semi-super cells we already
		* passed through.
		* c = height in A of one such semi-super cell.
		* Since cellDiv is always >=1, and divCount starts at 0, the actual position
		* of this atom with the super cell is given by:
		*/
		/* c*(real)((*muls).cellDiv-divCount-1) will pick the right super-cell
		* division in the big super-cell
		* The z-offset 0.5*cz[0] will position atoms at z=0 into the middle of the first
		* slice.
		*/
		atomZ = atoms[iatom].z-c*(real)(muls->cellDiv-divCount-1) + muls->czOffset
			-(0.5*muls->sliceThickness*(1-muls->centerSlices));
		// make sure that slices are centered for 2D and 3D differently:
		if (muls->potential3D==0)	atomZ += 0.5*muls->sliceThickness;
		else atomZ -= muls->sliceThickness;

		/* Now we need to find the first atom that contributes to this slice */
		/* Make use of the fact that we sorted the atoms in z */
		if ((*muls).nonPeriodZ) {
			if (((*muls).potential3D) && (atomZ -(*muls).atomRadius > c)) break;
			if (((*muls).potential3D==0) && (atomZ > c)) break;
			do {
				// printf("z: %g c: %g\n",atomZ,c);
				if (((*muls).potential3D) && (atomZ+(*muls).atomRadius+muls->sliceThickness >=0)) break;
				if (((*muls).potential3D==0) && (atomZ >=0)) break;
				// atomZ = atoms[++iatom].z-c*(real)((*muls).cellDiv-divCount-1)+ (0.5*(*muls).cz[0]*muls->centerSlices);
				atomZ = atoms[++iatom].z-c*(real)(muls->cellDiv-divCount-1) + muls->czOffset
					-(0.5*muls->sliceThickness*(1-muls->centerSlices));
				if (muls->potential3D==0)	atomZ += 0.5*muls->sliceThickness;
				else atomZ -= muls->sliceThickness;
			}
			while (iatom < natom-1);
		}
		/* atom coordinates in cartesian coords
		* The x- and y-position will be offset by the starting point
		* of the actually needed array of projected potential
		*/
		atomX = atoms[iatom].x -(*muls).potOffsetX;
		atomY = atoms[iatom].y -(*muls).potOffsetY;

		/* so far we need periodicity in z-direction.
		* This requirement can later be removed, if we
		* use some sort of residue slice which will contain the
		* proj. potential that we need to add to the first slice of
		* the next stack of slices
		*
		*
		*/

		/*************************************************************
		* real space potential lookup table summation
		************************************************************/
		if (!muls->fftpotential) {
			/* Warning: will assume constant slice thickness ! */
			/* do not round here: atomX=0..dx -> iAtomX=0 */
			/*
			iAtomX = (int)(atomX/dx);
			if (atomX/dx < (float)iAtomX) iAtomX--; // in case iAtomX is negative
			iAtomY = (int)(atomY/dy);
			if (atomY/dy < (float)iAtomY) iAtomY--;
			iAtomZ = (int)(atomZ/(*muls).cz[0]);
			if (atomZ/(*muls).cz[0] < (float)iAtomZ) iAtomZ--;
			*/
			iAtomX = (int)floor(atomX/dx);
			iAtomY = (int)floor(atomY/dy);
			iAtomZ = (int)floor(atomZ/muls->cz[0]);

			// printf("atomZ(%d)=%g(%d)\t",iatom,atomZ,iAtomZ);

			if (muls->displayPotCalcInterval > 0) {
				if ((muls->printLevel>=3) && ((iatom+1) % muls->displayPotCalcInterval == 0)) {
					printf("adding atom %d [%.3f %.3f %.3f (%.3f)], Z=%d\n",
						iatom+1,atomX+(*muls).potOffsetX,atomY+(*muls).potOffsetY,
						atoms[iatom].z,atomZ,atoms[iatom].Znum);
					/*    (*muls).ax,(*muls).by,(*muls).potOffsetX,(*muls).potOffsetY); */
				}
			}

			for (iax = -iRadX;iax<=iRadX;iax++) {
				if ((*muls).nonPeriod) {
					if (iax+iAtomX < 0) {
						iax = -iAtomX;
						if (abs(iax)>iRadX) break;
					}
					if (iax+iAtomX >= nx)	break;
				}
				x = (double)(iAtomX+iax)*dx-atomX;
				ix = (iax+iAtomX+16*nx) % nx;	/* shift into the positive range */
				for (iay=-iRadY;iay<=iRadY;iay++) {
					if ((*muls).nonPeriod) {
						if (iay+iAtomY < 0) {
							iay = -iAtomY;
							if (abs(iay)>iRadY) break;
						}
						if (iay+iAtomY >= ny)	break;
					}
					y = (double)(iAtomY+iay)*dy-atomY;
					iy = (iay+iAtomY+16*ny) % ny;	  /* shift into the positive range */
					r2sqr = x*x + y*y;

					if (r2sqr <= atomRadius2) {

						if (muls->potential3D) {
							/* calculate the range which we have left to cover with z-variation */
							/* iRadZ is the number of slices (rounded up) that this atom
							* will contribute to, given its current x,y-radius
							*/
							iRadZ = (int)(sqrt(atomRadius2-r2sqr)/(*muls).cz[0]+1.0);
							/* loop through the slices that this atoms contributes to */
							for (iaz=-iRadZ;iaz <=iRadZ;iaz++) {
								if ((*muls).nonPeriodZ) {
									if (iaz+iAtomZ < 0)  {
										if (-iAtomZ <= iRadZ) iaz = -iAtomZ;
										else break;
										if (abs(iaz)>nlayer) break;
									}
									if (iaz+iAtomZ >= nlayer)	break;
								}
								z = (double)(iAtomZ+iaz+0.5)*(*muls).cz[0]-atomZ;
								/* shift into the positive range */
								iz = (iaz+iAtomZ+32*nlayer) % nlayer;
								/* x,y,z is the true vector from the atom center
								* We can look up the proj potential at that spot
								* using trilinear extrapolation.
								*/
								atomBoxLookUp(&dPot,muls,atoms[iatom].Znum,x,y,z,
									muls->tds ? 0 : atoms[iatom].dw);
								//    printf("access: %d %d %d\n",iz,ix,iy);
								muls->trans[iz][ix][iy][0] += dPot[0];
								muls->trans[iz][ix][iy][1] += dPot[1];
							} /* end of for iaz=-iRadZ .. iRadZ */
						} /* end of if potential3D */

						/********************************************************************/

						else { /* if 2D potential */
							if ((*muls).nonPeriodZ) {
								if (iAtomZ < 0)  break;
								if (iAtomZ >= nlayer)	break;
							}
							iz = (iAtomZ+32*nlayer) % nlayer;	  /* shift into the positive range */
							atomBoxLookUp(&dPot,muls,atoms[iatom].Znum,x,y,0,
								muls->tds ? 0 : atoms[iatom].dw);
							z = (double)(iAtomZ+1)*(*muls).cz[0]-atomZ;

							/*
							printf("iz=%d,%d, z=%g, atomZ=%g(%g), iAtomZ=%d, c=%g (%d, %d))\n",
							iz,iatom,z,atomZ,atoms[iatom].z,iAtomZ,c, divCount,(*muls).cellDiv);
							*/
							/* split the atom if it is close to the top edge of the slice */
							//    printf("access: %d %d %d\n",iz,ix,iy);
							if ((z<0.15*(*muls).cz[0]) && (iz >0)) {
								muls->trans[iz][ix][iy][0] += 0.5*dPot[0];
								muls->trans[iz][ix][iy][1] += 0.5*dPot[1];
								muls->trans[iz-1][ix][iy][0] += 0.5*dPot[0];
								muls->trans[iz-1][ix][iy][1] += 0.5*dPot[1];
							}
							/* split the atom if it is close to the bottom edge of the slice */
							else {
								if ((z>0.85*(*muls).cz[0]) && (iz < nlayer-1)) {
									(*muls).trans[iz][ix][iy][0] += 0.5*dPot[0];
									(*muls).trans[iz][ix][iy][1] += 0.5*dPot[1];
									(*muls).trans[iz+1][ix][iy][0] += 0.5*dPot[0];
									(*muls).trans[iz+1][ix][iy][1] += 0.5*dPot[1];
								}
								else {
									(*muls).trans[iz][ix][iy][0] += dPot[0];
									(*muls).trans[iz][ix][iy][1] += dPot[1];
								}
							}
							/*
							if ((iatom % 100 == 0) && (fabs(x) <0.04) && (fabs(y) < 0.04))
							printf("i: %d iz: %d atomZ: %g (%g)\n",iatom,iz,atomZ,z);
							*/
						}
					}
				}
			}
		}


		/**************************************************************************
		* Newer, even faster method based on FFT of tabulated scattering factors
		**************************************************************************/
		else {	/* fftpotential */
			iAtomX = (int)floor(atomX/dx);
			iAtomY = (int)floor(atomY/dy);
			if (muls->potential3D) atomZ+=muls->sliceThickness;  // why ??? !!!!!
			// printf("%d: pos=[%d, %d, %.1f]\n",iatom,iAtomX,iAtomY,atomZ);
			// atomZ is z-distance with respect to the start of the current stack of slices.
			// ddz = atomZ-dz*iAtomZ;


			/////////////////////////////////////////////////////////////////////////////
			// if we need to cut away at the edges, i.e. non-periodic potential arrays:
			if (muls->nonPeriod) {
				if (muls->potential3D) {
					iAtomZ = (int)floor(atomZ/muls->sliceThickness+0.5);
					// printf("iAtomZ: %d\n",iAtomZ);

					iax0 = iAtomX-iRadX <  0 ? 0 : iAtomX-iRadX;
					iax1 = iAtomX+iRadX >= muls->potNx ? muls->potNx-1 : iAtomX+iRadX;
					iay0 = iAtomY-iRadY <  0 ? 0 : iAtomY-iRadY;
					iay1 = iAtomY+iRadY >= muls->potNy ? muls->potNy-1 : iAtomY+iRadY;
					// if within the potential map range:
					if ((iax0 <  muls->potNx) && (iax1 >= 0) && (iay0 <  muls->potNy) && (iay1 >= 0)) {


						// define range of sampling from atomZ-/+atomRadius
						iaz0 = iAtomZ-iRadZ <  0 ? -iAtomZ : -iRadZ;
						iaz1 = iAtomZ+iRadZ >= muls->slices ?  muls->slices-iAtomZ-1 : iRadZ;
						// iaz0 = 0;  iaz1 = 0;
						// printf("iatomZ: %d, %d..%d cz=%g, %g (%d), dOffsZ=%g (%d)\n",iAtomZ,iaz0,iaz1,muls->sliceThickness,atomZ,(int)atomZ,(iAtomZ+iaz0-atomZ/muls->sliceThickness)*nzSub,(int)(iAtomZ+iaz0-atomZ/muls->sliceThickness)*nzSub+0.5);
						if ((iAtomZ+iaz0 <	muls->slices) && (iAtomZ+iaz1 >= 0)) {
							// retrieve the pointer for this atom
							atPotPtr     = getAtomPotential3D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut);
#if USE_Q_POT_OFFSETS
							// retrieve the pointer to the array of charge-dependent potential offset
							// This function will return NULL; if the charge of this atom is zero:
							atPotOffsPtr = getAtomPotentialOffset3D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut,atoms[iatom].q);
#endif // USE_Q_POT_OFFSETS
							iOffsLimHi   =  Nr*(Nz_lut-1);
							iOffsLimLo   = -Nr*(Nz_lut-1);
							iOffsStep    = nzSub*Nr;

							// Slices around the slice that this atom is located in must be affected by this atom:
							// iaz must be relative to the first slice of the atom potential box.
							for (iax=iax0; iax <= iax1; iax++) {
								potPtr = &(muls->trans[iAtomZ+iaz0][iax][iay0][0]);
								// potPtr = &(muls->trans[iAtomZ-iaz0+iaz][iax][iay0][0]);
								// printf("access: %d %d %d (%d)\n",iAtomZ+iaz0,iax,iay0,(int)potPtr);

								//////////////////////////////////////////////////////////////////////
								// Computation of Radius must be made faster by using pre-calculated ddx
								// and LUT for sqrt:
								// use of sqrt slows down from 120sec to 180 sec.
								x2 = iax*dx - atomX;  x2 *= x2;
								for (iay=iay0; iay <= iay1; iay++) {
									// printf("iax=%d, iay=%d\n",iax,iay);
									y2 = iay*dy - atomY;  y2 *= y2;
									r = sqrt(x2+y2);
									// r = (x2+y2);
									ddr = r/dr;
									ir	= (int)floor(ddr);
									// add in different slices, once r has been defined
									if (ir < Nr-1) {
										ddr = ddr-(double)ir;
										ptr = potPtr;

										dOffsZ = (iAtomZ+iaz0-atomZ/muls->sliceThickness)*nzSub;
#if Z_INTERPOLATION
										iOffsZ = (int)dOffsZ;
										ddz    = fabs(dOffsZ - (double)iOffsZ);
#else
										iOffsZ = (int)(dOffsZ+0.5);
#endif
										iOffsZ *= Nr;

										for (iaz=iaz0; iaz <= iaz1; iaz++) {
											potVal = 0;
											// iOffsZ = (int)(fabs(iAtomZ+iaz-atomZ/muls->sliceThickness)*nzSub+0.5);
											if (iOffsZ < 0) {
												if (iOffsZ > iOffsLimLo) {
													// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
													potVal = (1-ddz)*((1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0])+
														ddz *((1-ddr)*atPotPtr[ir-iOffsZ   ][0]+ddr*atPotPtr[ir+1-iOffsZ   ][0]);
#if USE_Q_POT_OFFSETS
													// add the charge-dependent potential offset
													if (atPotOffsPtr != NULL) {
														potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0])+
															ddz *((1-ddr)*atPotOffsPtr[ir-iOffsZ   ][0]+ddr*atPotOffsPtr[ir+1-iOffsZ   ][0]));
													}
#endif  // USE_Q_POT_OFFSETS
#else   // Z_INTERPOLATION
													potVal = (1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0];
#if USE_Q_POT_OFFSETS
													// add the charge-dependent potential offset
													if (atPotOffsPtr != NULL) {
														potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0]);
													}
#endif  // USE_Q_POT_OFFSETS
#endif  // Z_INTERPOLATION
													// if (r < dr) printf("%d: iaz=%d,pot=%g,ddr=%g, iOffsZ=%d\n",iatom,iaz,potVal,ddr,iOffsZ);
												}
											} // if iOffZ < 0
											else {
												// select the pointer to the right layer in the lookup box
												// printf("%4d: iOffsZ: %d, iaz: %d (slice: %d, pos: %g [%d .. %d])\n",iatom,iOffsZ,iaz,iAtomZ+iaz,atomZ,iaz0,iaz1);
												if (iOffsZ < iOffsLimHi) {
													// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
													potVal = (1-ddz)*((1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0])+
														ddz *((1-ddr)*atPotPtr[ir+iOffsZ+Nr][0]+ddr*atPotPtr[ir+1+iOffsZ+Nr][0]);
#if USE_Q_POT_OFFSETS
													// add the charge-dependent potential offset
													if (atPotOffsPtr != NULL) {
														potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0])+
															ddz *((1-ddr)*atPotOffsPtr[ir+iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1+iOffsZ+Nr][0]));
													}
#endif  // USE_Q_POT_OFFSETS
#else   // Z_INTERPOLATION
													potVal = (1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0];
#if USE_Q_POT_OFFSETS
													// add the charge-dependent potential offset
													if (atPotOffsPtr != NULL) {
														potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0]);
													}
#endif  // USE_Q_POT_OFFSETS
#endif  // Z_INTERPOLATION
													// if (r < dr) printf("%d: iaz=%d,pot=%g,ddr=%g, iOffsZ=%d\n",iatom,iaz,potVal,ddr,iOffsZ);

												}
											} // if iOffsZ >=0
											*ptr += potVal;  // ptr = potPtr = muls->trans[...]

											ptr  += sliceStep;	// advance to the next slice
											// add the remaining potential to the next slice:
											// if (iaz < iaz1)	*ptr += (1-ddz)*potVal;
											iOffsZ += iOffsStep;
										} // end of iaz-loop
									} // if ir < Nr
									// advance pointer to next complex potential point:
									potPtr+=2;
									// make imaginary part zero for now
									// *potPtr = 0;
									// potPtr++;

									// wrap around when end of y-line is reached:
									// if (++iay % muls->potNy == 0) potPtr -= 2*muls->potNy;
								} // iay=iay0 .. iay1
							} // iax=iax0 .. iax1
						} // iaz0+iAtomZ < muls->slices
						// dOffsZ = (iAtomZ-atomZ/muls->sliceThickness)*nzSub;
						// printf("%5d (%2d): iAtomZ=%2d, offsZ=%g, diff=%g, (%g)\n",
						//	  iatom,atoms[iatom].Znum,iAtomZ,dOffsZ,dOffsZ - (int)(dOffsZ+0.5),iAtomZ-atomZ/muls->sliceThickness);
					} // if within bounds
				}  // muls->potential3D and non-periodic in x-y
				////////////////////////////////////////////////////////////////////
				// 2D potential calculation already seems to work!
				// However, the potential is calculated wrongly, it is not integrated
				// in z-direction.	This must be done in the potential slice initialization
				// procedure.
				else {
					iAtomZ = (int)floor(atomZ/muls->sliceThickness);
					iax0 = iAtomX-iRadX <  0 ? 0 : iAtomX-iRadX;
					iax1 = iAtomX+iRadX >= muls->potNx ? muls->potNx-1 : iAtomX+iRadX;
					iay0 = iAtomY-iRadY <  0 ? 0 : iAtomY-iRadY;
					iay1 = iAtomY+iRadY >= muls->potNy ? muls->potNy-1 : iAtomY+iRadY;
					// if within the potential map range:
					if ((iax0 <  muls->potNx) && (iax1 >= 0) && (iay0 <  muls->potNy) && (iay1 >= 0)) {
						ddx = (-(double)iax0+(atomX/dx-(double)iRadX))*(double)OVERSAMP_X;
						ddy = (-(double)iay0+(atomY/dy-(double)iRadY))*(double)OVERSAMP_X;
						iOffsX = (int)floor(ddx);
						iOffsY = (int)floor(ddy);
						ddx -= (double)iOffsX;
						ddy -= (double)iOffsY;
						s11 = (1-ddx)*(1-ddy);
						s12 = (1-ddx)*ddy;
						s21 = ddx*(1-ddy);
						s22 = ddx*ddy;
						atPotPtr = getAtomPotential2D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw);

						for (iax=iax0; iax < iax1; iax++) {
							// printf("(%d, %d): %d,%d\n",iax,nyAtBox,(iOffsX+OVERSAMP_X*(iax-iax0)),iOffsY+iay1-iay0);
							// potPtr and ptr are of type (float *)
							potPtr = &(muls->trans[iAtomZ][iax][iay0][0]);
							ptr = &(atPotPtr[(iOffsX+OVERSAMP_X*(iax-iax0))*nyAtBox+iOffsY][0]);
							for (iay=iay0; iay < iay1; iay++) {
								*potPtr += s11*(*ptr)+s12*(*(ptr+2))+s21*(*(ptr+nyAtBox2))+s22*(*(ptr+nyAtBox2+2));

								potPtr++;
								// *potPtr = 0;
								potPtr++;
								ptr += 2*OVERSAMP_X;
							}
						}

					}  // not muls->potential3D
				}
			} // end of if not periodic

			/////////////////////////////////////////////////////////////////////////////
			// if the potential array is periodic:
			else {


				// printf("Z=%d (z=%d), iOffs: %d, %d (%g, %g) %g, %g, %g %g, atom: %g, %g\n",
				//		atoms[iatom].Znum,iAtomZ,iOffsX,iOffsY,ddx,ddy,s11,s12,s21,s22,atomX,atomY);
				////////////////////////////////////////////////////////////////////
				// add code here!
				if (muls->potential3D) {
					iAtomZ = (int)floor(atomZ/muls->sliceThickness+0.5);
					iax0 = iAtomX-iRadX;
					iax1 = iAtomX+iRadX;
					iay0 = iAtomY-iRadY;
					iay1 = iAtomY+iRadY;

					// define range of sampling from atomZ-/+atomRadius
					iaz0 = iAtomZ-iRadZ <  0 ? -iAtomZ : -iRadZ;
					iaz1 = iAtomZ+iRadZ >= muls->slices ?  muls->slices-iAtomZ-1 : iRadZ;
					// if (iatom < 2) printf("iatomZ: %d, %d cz=%g, %g: %d, %d\n",iAtomZ,iaz0,muls->sliceThickness,atomZ,(int)(-2.5-atomZ),(int)(atomZ+2.5));
					// printf("%d: iatomZ: %d, %d cz=%g, %g\n",iatom,iAtomZ,iaz0,muls->sliceThickness,atomZ);
					if ((iAtomZ+iaz0 <  muls->slices) && (iAtomZ+iaz1 >= 0)) {
						// retrieve the pointer for this atom
						atPotPtr = getAtomPotential3D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut);
#if USE_Q_POT_OFFSETS
						// retrieve the pointer to the array of charge-dependent potential offset
						// This function will return NULL; if the charge of this atom is zero:
						atPotOffsPtr = getAtomPotentialOffset3D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut,atoms[iatom].q);
#endif // USE_Q_POT_OFFSETS
						iOffsLimHi =	Nr*(Nz_lut-1);
						iOffsLimLo = -Nr*(Nz_lut-1);
						iOffsStep  = nzSub*Nr;

						// Slices around the slice that this atom is located in must be affected by this atom:
						// iaz must be relative to the first slice of the atom potential box.
						for (iax=iax0; iax < iax1; iax++) {
							potPtr = &(muls->trans[iAtomZ+iaz0][(iax+2*muls->potNx) % muls->potNx][(iay0+2*muls->potNy) % muls->potNy][0]);
							// potPtr = &(muls->trans[iAtomZ-iaz0+iaz][iax][iay0][0]);
							x2 = iax*dx - atomX;	x2 *= x2;
							for (iay=iay0; iay < iay1; ) {
								// printf("iax=%d, iay=%d\n",iax,iay);
								y2 = iay*dy - atomY;	y2 *= y2;
								r = sqrt(x2+y2);
								ddr = r/dr;
								ir  = (int)floor(ddr);
								// add in different slices, once r has been defined
								if (ir < Nr-1) {
									ddr = ddr-(double)ir;
									ptr = potPtr;
									// Include interpolation in z-direction as well (may do it in a very smart way here !):

									dOffsZ = (iAtomZ+iaz0-atomZ/muls->sliceThickness)*nzSub;
#if Z_INTERPOLATION
									iOffsZ = (int)dOffsZ;
									ddz	 = fabs(dOffsZ - (double)iOffsZ);
#else  // Z_INTERPOLATION
									iOffsZ = (int)(dOffsZ+0.5);
#endif  // Z_INTERPOLATION
									iOffsZ *= Nr;

									for (iaz=iaz0; iaz <= iaz1; iaz++) {
										potVal = 0;
										// iOffsZ = (int)(fabs(iAtomZ+iaz-atomZ/muls->sliceThickness)*nzSub+0.5);
										if (iOffsZ < 0) {
											if (iOffsZ > iOffsLimLo) {
												// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
												potVal = (1-ddz)*((1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0])+
													ddz *((1-ddr)*atPotPtr[ir-iOffsZ	 ][0]+ddr*atPotPtr[ir+1-iOffsZ	 ][0]);
#if USE_Q_POT_OFFSETS
												// add the charge-dependent potential offset
												if (atPotOffsPtr != NULL) {
													potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0])+
														ddz *((1-ddr)*atPotOffsPtr[ir-iOffsZ   ][0]+ddr*atPotOffsPtr[ir+1-iOffsZ   ][0]));
												}
#endif  // USE_Q_POT_OFFSETS
#else  // Z_INTERPOLATION
												potVal = (1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0];
#if USE_Q_POT_OFFSETS
												// add the charge-dependent potential offset
												if (atPotOffsPtr != NULL) {
													potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0]);
												}
#endif  // USE_Q_POT_OFFSETS
#endif  // Z_INTERPOLATION
											}
										}
										else {
											// select the pointer to the right layer in the lookup box
											// printf("%4d: iOffsZ: %d, iaz: %d (slice: %d, pos: %g [%d .. %d])\n",iatom,iOffsZ,iaz,iAtomZ+iaz,atomZ,iaz0,iaz1);
											if (iOffsZ < iOffsLimHi) {
												// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
												potVal = (1-ddz)*((1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0])+
													ddz *((1-ddr)*atPotPtr[ir+iOffsZ+Nr][0]+ddr*atPotPtr[ir+1+iOffsZ+Nr][0]);
#if USE_Q_POT_OFFSETS
												// add the charge-dependent potential offset
												if (atPotOffsPtr != NULL) {
													potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0])+
														ddz *((1-ddr)*atPotOffsPtr[ir+iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1+iOffsZ+Nr][0]));
												}
#endif  // USE_Q_POT_OFFSETS
#else  // Z_INTERPOLATION
												potVal = (1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0];
#if USE_Q_POT_OFFSETS
												// add the charge-dependent potential offset
												if (atPotOffsPtr != NULL) {
													potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0]);
												}
#endif  // USE_Q_POT_OFFSETS
#endif  // Z_INTERPOLATION
											}
										}
										*ptr += potVal;

										ptr  += sliceStep;  // advance to the next slice
										// add the remaining potential to the next slice:
										// if (iaz < iaz1)  *ptr += (1-ddz)*potVal;
										iOffsZ += iOffsStep;
									}
								} // if ir < Nr-1
								// advance pointer to next complex potential point:
								potPtr+=2;
								// make imaginary part zero for now
								// *potPtr = 0;
								// potPtr++;

								// wrap around when end of y-line is reached:
								if (++iay % muls->potNy == 0) potPtr -= 2*muls->potNy;
							}
						}
					} // iaz0+iAtomZ < muls->slices
				}  // muls->potential3D
				////////////////////////////////////////////////////////////////////
				// 2D potential (periodic) calculation already seems to work!
				else {
					iAtomZ = (int)floor(atomZ/muls->sliceThickness);
					iax0 = iAtomX-iRadX+2*muls->potNx;
					iax1 = iAtomX+iRadX+2*muls->potNx;
					iay0 = iAtomY-iRadY+2*muls->potNy;
					iay1 = iAtomY+iRadY+2*muls->potNy;

					ddx = (-(double)iax0+(atomX/dx-(double)(iRadX-2*muls->potNx)))*(double)OVERSAMP_X;
					ddy = (-(double)iay0+(atomY/dy-(double)(iRadY-2*muls->potNy)))*(double)OVERSAMP_X;
					iOffsX = (int)floor(ddx);
					iOffsY = (int)floor(ddy);
					ddx -= (double)iOffsX;
					ddy -= (double)iOffsY;
					/*
					s11 = (1-ddx)*(1-ddy);
					s12 = (1-ddx)*ddy;
					s21 = ddx*(1-ddy);
					s22 = ddx*ddy;
					*/
				    s22 = (1-ddx)*(1-ddy);
					s21 = (1-ddx)*ddy;
					s12 = ddx*(1-ddy);
					s11 = ddx*ddy;

					atPotPtr = getAtomPotential2D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw);

					// if (iatom < 3) printf("atom #%d: ddx=%g, ddy=%g iatomZ=%d, atomZ=%g, %g\n",iatom,ddx,ddy,iAtomZ,atomZ,atoms[iatom].z);
					for (iax=iax0; iax < iax1; iax++) {  // TODO: should use ix += OVERSAMP_X
						// printf("(%d, %d): %d,%d\n",iax,nyAtBox,(iOffsX+OVERSAMP_X*(iax-iax0)),iOffsY+iay1-iay0);
						// potPtr and ptr are of type (float *)
						//////////////////
						// Only the exact slice that this atom is located in is affected by this atom:
						int atPosX = (OVERSAMP_X*(iax-iax0)-iOffsX);
						if ((atPosX >= 0) && (atPosX < nyAtBox-1)) {
							ptr = &(atPotPtr[atPosX*nyAtBox-iOffsY][0]);
						for (iay=iay0; iay < iay1; iay++) {
							// wrap around when end of y-line is reached:
								int atPosY = (iay-iay0)*OVERSAMP_X-iOffsY;
								if ((atPosY < nyAtBox-1) && (atPosY >=0)) {
							// do the real part
									muls->trans[iAtomZ][iax % muls->potNx][iay % muls->potNy][0] +=
									     s11*(*ptr)+s12*(*(ptr+2))+s21*(*(ptr+nyAtBox2))+s22*(*(ptr+nyAtBox2+2));
								}
							// make imaginary part zero for now
								// *potPtr = 0;  potPtr++;
							ptr += 2*OVERSAMP_X;
						}
						} // if atPosX within limits
					}
				}  // muls->potential3D	== 0
			}
			////////////////////////////////////////////////////////////////////
		} /* end of if (fftpotential) */
	} /* for iatom =0 ... */
	time(&time1);
	if (iatom > 0)
	if (muls->printLevel) printf("%g sec used for real space potential calculation (%g sec per atom)\n",difftime(time1,time0),difftime(time1,time0)/iatom);
	else
	if (muls->printLevel) printf("%g sec used for real space potential calculation\n",difftime(time1,time0));


	/*************************************************/
	/* Save the potential slices					   */

	if (muls->savePotential) {
		for (iz = 0;iz<nlayer;iz++){
			/*
			muls->thickness = iz;
			showCrossSection(muls,(*muls).transr[iz],nx,1,0);
			*/
			// find the maximum value of each layer:
			potVal = muls->trans[iz][0][0][0];
			for (ddx=potVal,ddy = potVal,ix=0;ix<muls->potNy*muls->potNx;potVal = muls->trans[iz][0][++ix][0]) {
				if (ddy<potVal) ddy = potVal;
				if (ddx>potVal) ddx = potVal;
			}

#ifndef WIN32
			sprintf(fileOut,"%s/%s%d.img",muls->folder,muls->fileBase,iz);
#else
			sprintf(fileOut,"%s\\%s%d.img",muls->folder,muls->fileBase,iz);
#endif
			if (muls->printLevel >= 3)
				printf( "Saving (complex) potential layer %d to file %s (r: %g..%g)\n", iz, fileOut, ddx, ddy );

			imageIO->SetThickness(muls->sliceThickness);
			sprintf(buf,"Projected Potential (slice %d)",iz);
			imageIO->SetComment(buf);
			imageIO->WriteComplexImage( (void **)muls->trans[iz], fileOut );
		} // loop through all slices
	} /* end of if savePotential ... */
	if (muls->saveTotalPotential) {
		if (tempPot == NULL) tempPot = float2D(muls->potNx,muls->potNy,"total projected potential");

		for (ix=0;ix<muls->potNx;ix++) for (iy=0;iy<muls->potNy;iy++) {
			tempPot[ix][iy] = 0;
			for (iz=0;iz<nlayer;iz++) tempPot[ix][iy] += muls->trans[iz][ix][iy][0];
		}

		for (ddx=tempPot[0][0],ddy = potVal,ix=0;ix<muls->potNy*muls->potNx;potVal = tempPot[0][++ix]) {
			if (ddy<potVal) ddy = potVal;
			if (ddx>potVal) ddx = potVal;
		}
#ifndef WIN32
		sprintf(fileOut,"%s/%sProj.img",muls->folder,muls->fileBase);
#else
		// RAM DEBUG : this is overwriting the original config file, why?  is filename something else?
		// RAM RESOLVED: both fileName and filename were defined above...
		sprintf( fileOut, "%s\\%sProj.img", muls->folder, muls->fileBase );
#endif
		if (muls->printLevel >= 2)
			printf( "Saving total projected potential to file %s (r: %g..%g)\n", fileOut, ddx, ddy );
		imageIO->SetThickness(nlayer*muls->sliceThickness);
		sprintf(buf,"Projected Potential (sum of %d slices)",muls->slices);
		imageIO->SetComment(buf);
		imageIO->WriteRealImage( (void **)tempPot, fileOut );
	}

} // end of make3DSlices



/********************************************************************************
* Create Lookup table for 3D potential due to neutral atoms
********************************************************************************/
#define PHI_SCALE 47.87658
fftwf_complex *getAtomPotential3D(int Znum, MULS *muls,double B,int *nzSub,int *Nr,int *Nz_lut) {
	int ix,iy,iz,iiz,ind3d,iKind,izOffset;
	double zScale,kzmax,zPos,xPos;
	fftwf_plan plan;
	static double f,phase,s2,s3,kmax2,smax2,kx,kz,dkx,dky,dkz; // ,dx2,dy2,dz2;
	static int nx,ny,nz,nzPerSlice;
	static fftwf_complex **atPot = NULL;
	static fftwf_complex *temp = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
	ImageIOPtr imageio = ImageIOPtr();
	static fftwf_complex *ptr = NULL;
	static char fileName[256];
#endif
	static double *splinb=NULL;
	static double *splinc=NULL;
	static double *splind=NULL;


	// scattering factors in:
	// float scatPar[4][30]
	if (atPot == NULL) {
		splinb = double1D(N_SF, "splinb" );
		splinc = double1D(N_SF, "splinc" );
		splind = double1D(N_SF, "splind" );


		nx = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX);
		ny = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY);
		// The FFT-resolution in the z-direction must be high enough to avoid
		// artifacts due to premature cutoff of the rec. space scattering factor
		// we will therefore make it roughly the same as the x-resolution
		// However, we will make sure that a single slice contains an integer number
		// of sampling points.
		nzPerSlice = (int)floor(OVERSAMP_X*muls->sliceThickness/muls->resolutionX);
		// make nzPerSlice odd:
		if (2.0*(nzPerSlice >> 1) == nzPerSlice) nzPerSlice += 1;
		// Total number of z-positions should be twice that of atomRadius/sliceThickness
		nz = (2*(int)ceil(muls->atomRadius/muls->sliceThickness))*nzPerSlice;
		if (muls->printLevel > 1) printf("Will use %d sampling points per slice, total nz=%d (%d)\n",nzPerSlice,nz,nzPerSlice >> 1);

		dkx = 0.5*OVERSAMP_X/(nx*muls->resolutionX);  // nx*muls->resolutionX is roughly 2*muls->atomRadius
		dky = dkx;
		dkz = nzPerSlice/(double)(nz*muls->sliceThickness);
		kmax2 = 0.5*nx*dkx/(double)OVERSAMP_X;  // largest k that we'll admit
		smax2 = kmax2;

		//printf("dkx = %g, nx = %d, kmax2 = %g\n",dkx,nx,kmax2);
		if (muls->printLevel > 1) printf("Cutoff scattering angle: kmax=%g (1/A), dk=(%g,%g %g)\n",kmax2,dkx,dky,dkz);
		scatPar[0][N_SF-1] = 1.2*smax2;
		scatPar[0][N_SF-2] = 1.1*smax2;
		scatPar[0][N_SF-3] = smax2;
		// adjust the resolution of the lookup table if necessary
		if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3]) {

			if (1) {
				// set additional scattering parameters to zero:
				for (ix = 0;ix < N_SF-10;ix++) {
					if (scatPar[0][N_SF-4-ix] < scatPar[0][N_SF-3]-0.001*(ix+1)) break;
					scatPar[0][N_SF-4-ix] = scatPar[0][N_SF-3]-0.001*(ix+1);
					for (iy=1; iy<N_ELEM;iy++) scatPar[iy][N_SF-4-ix] = 0;
				}
			}
			else {
				for (ix = 0;ix < 20;ix++) {
					if (scatPar[0][N_SF-4-ix] < scatPar[0][N_SF-3]) break;
					scatPar[0][N_SF-4-ix] = scatPar[0][N_SF-3];
					for (iy=1; iy<N_ELEM;iy++) scatPar[iy][N_SF-4-ix] = scatPar[iy][N_SF-3];
				}
			}
			if (muls->printLevel > 1) printf("getAtomPotential3D: set resolution of scattering factor to %g/A!\n",
				scatPar[0][N_SF-4-ix]);
		}	// end of if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3])
		smax2 *= smax2;
		kmax2 *= kmax2;
		// allocate a list of pointers for the element-specific potential lookup table
		atPot = (fftwf_complex **)malloc((NZMAX+1)*sizeof(fftwf_complex *));
		for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
		temp  = (fftwf_complex*) fftwf_malloc(nx*nz*sizeof(fftwf_complex));
	}
	// initialize this atom, if it has not been done yet:
	if (atPot[Znum] == NULL) {
		iKind = Znum;

		// setup cubic spline interpolation:
		splinh(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF);

		// allocate a 3D array:
		atPot[Znum] = (fftwf_complex*) fftwf_malloc(nx*nz/4*sizeof(fftwf_complex));
		memset(temp,0,nx*nz*sizeof(fftwf_complex));
		kzmax	  = dkz*nz/2.0;
		// define x-and z-position of atom center:
		// The atom should sit in the top-left corner,
		// however (nzPerSlice+1)/2 above zero in z-direction
		xPos = -2.0*PI*0.0;  // or muls->resolutionX*nx/(OVERSAMP_X), if in center
		izOffset = (nzPerSlice-1)/2;
		zPos = -2.0*PI*(muls->sliceThickness/nzPerSlice*(izOffset));

		// What this look-up procedure will do is to supply V(r,z) computed from fe(q).
		// Since V(r,z) is rotationally symmetric we might as well compute
		// V(x,y,z) at y=0, i.e. V(x,z).
		// In order to do the proper 3D inverse FT without having to do a complete 3D FFT
		// we will pre-compute the qy-integral for y=0.

		// kzborder = dkz*(nz/(2*OVERSAMP_Z) -1);
		for (iz=0;iz<nz;iz++) {
			kz = dkz*(iz<nz/2 ? iz : iz-nz);
			// We also need to taper off the potential in z-direction
			// in order to avoid cutoff artifacts.
			// zScale = fabs(kz) <= kzborder ? 1.0 : 0.5+0.5*cos(M_PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
			// printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
			for (ix=0;ix<nx;ix++) {
				kx = dkx*(ix<nx/2 ? ix : ix-nx);
				s2 = (kx*kx+kz*kz);
				// if this is within the allowed circle:
				if (s2<smax2) {
					ind3d = ix+iz*nx;
					// f = fe3D(Znum,k2,muls->tds,1.0,muls->scatFactor);
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					f = seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s2))*exp(-s2*B*0.25);
					// perform the qy-integration for qy <> 0:
					for (iy=1;iy<nx;iy++) {
						s3 = dkx*iy;
						s3 = s3*s3+s2;
						if (s3<smax2) {
							f += 2*seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s3))*exp(-s3*B*0.25);
						}
						else break;
					}
					f *= dkx;
					// note that the factor 2 is missing in the phase (2pi k*r)
					// this places the atoms in the center of the box.
					phase	= kx*xPos + kz*zPos;
					temp[ind3d][0] = f*cos(phase);  // *zScale
					temp[ind3d][1] = f*sin(phase);  // *zScale
					// if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
				}
			}
		} // for iz ...


#if SHOW_SINGLE_POTENTIAL
		// 0 thickness
		imageio = ImageIOPtr(new CImageIO(nz, nx, 0, dkz, dkx));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		sprintf(fileName,"pot_rec_%d.img",Znum);
		imageio->SetThickness(muls->sliceThickness);
		imageio->WriteComplexImage((void**)temp, fileName);
#endif
		// This converts the 2D kx-kz  map of the scattering factor to a 2D real space map.
		plan = fftwf_plan_dft_2d(nz,nx,temp,temp,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
		// dx2 = muls->resolutionX*muls->resolutionX/(OVERSAMP_X*OVERSAMP_X);
		// dy2 = muls->resolutionY*muls->resolutionY/(OVERSAMP_X*OVERSAMP_X);
		// dz2 = muls->sliceThickness*muls->sliceThickness/(OVERSAMP_Z*OVERSAMP_Z);
		// We also make sure that the potential touches zero at least somewhere.  This will avoid
		// sharp edges that could produce ringing artifacts.
		// It is certainly debatable whether this is a good apprach, or not.
		// printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
		// min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
		for (ix=0;ix<nx/2;ix++)  for (iz=0;iz<nz/2;iz++) {
			ind3d = ix+iz*nx/2;
			// Integrate over nzPerSlice neighboring layers here:::::::::
			for (zScale=0,iiz=-izOffset;iiz<=izOffset;iiz++) {
				if (iz+izOffset+iiz < nz/2) zScale += temp[ix+(iz+izOffset+iiz)*nx][0];
			}
			if (zScale < 0) zScale = 0;
			// assign the iz-th slice the sum of the 3 other slices:
			// and divide by unit cell volume (since this is in 3D):
			// Where does the '1/2' come from???  OVERSAMP_X*OVERSAMP_Y/8 = 1/2
			// if nothing has changed, then OVERSAMP_X=2 OVERSAMP_Z=18.
			// remember, that s=0.5*k;
			// This potential will later again be scaled by lambda*gamma (=0.025*1.39139)
			atPot[Znum][ind3d][0] = 47.8658*dkx*dkz/(nz)*zScale;

			// *8*14.4*0.529=4*a0*e (s. Kirkland's book, p. 207)
			// 2*pi*14.4*0.529 = 7.6176;
			// if (atPot[Znum][ind3d][0] < min) min = atPot[Znum][ind3d][0];
			atPot[Znum][ind3d][1]= 0;
		}
		// make sure we don't produce negative potential:
		// if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny][0] -= min;
#if SHOW_SINGLE_POTENTIAL
		imageio = ImageIOPtr(new CImageIO(nz/2, nx/2, 0, muls->sliceThickness/nzPerSlice,
			muls->resolutionX/OVERSAMP_X));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		imageio->SetThickness(nz*muls->sliceThickness/nzPerSlice);
		sprintf(fileName,"potential_rz_%d.img",Znum);
		ptr = atPot[Znum];
		imageio->WriteComplexImage((void**)ptr, fileName);
#endif
		if (muls->printLevel > 1) printf("Created 3D (r-z) %d x %d potential array for Z=%d (%d, B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)\n",
			nx/2,nz/2,Znum,iKind,B,dkx,dky,dkz,izOffset);
	}
	*Nz_lut = nz/2;
	*nzSub  = nzPerSlice;
	*Nr	  = nx/2;
	return atPot[Znum];
}



/********************************************************************************
* Lookup function for 3D potential offset due to charged atoms (ions)
********************************************************************************/
fftwf_complex *getAtomPotentialOffset3D(int Znum, MULS *muls,double B,int *nzSub,int *Nr,int *Nz_lut,float q) {
	int ix,iy,iz,iiz,ind3d,iKind,izOffset;
	double zScale,kzmax,zPos,xPos;
	fftwf_plan plan;
	static double f,phase,s2,s3,kmax2,kx,kz,dkx,dky,dkz; // ,dx2,dy2,dz2;
	static int nx,ny,nz,nzPerSlice;
	static fftwf_complex **atPot = NULL;
	static fftwf_complex *temp = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
	ImageIOPtr imageio = ImageIOPtr();
	static fftwf_complex *ptr = NULL;
	static char fileName[256];
#endif
	static double *splinb=NULL;
	static double *splinc=NULL;
	static double *splind=NULL;


	// if there is no charge to this atom, return NULL:
	if (q == 0) return NULL;


	// scattering factors in:
	// float scatPar[4][30]
	if (atPot == NULL) {
		splinb = double1D(N_SF, "splinb" );
		splinc = double1D(N_SF, "splinc" );
		splind = double1D(N_SF, "splind" );


		nx = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX);
		ny = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY);
		// The FFT-resolution in the z-direction must be high enough to avoid
		// artifacts due to premature cutoff of the rec. space scattering factor
		// we will therefore make it roughly the same as the x-resolution
		// However, we will make sure that a single slice contains an integer number
		// of sampling points.
		nzPerSlice = (int)floor(OVERSAMP_X*muls->sliceThickness/muls->resolutionX);
		// make nzPerSlice odd:
		if (2.0*floor((double)(nzPerSlice >> 1)) == nzPerSlice) nzPerSlice += 1;
		// Total number of z-positions should be twice that of atomRadius/sliceThickness
		nz = (2*(int)ceil(muls->atomRadius/muls->sliceThickness))*nzPerSlice;
		if (muls->printLevel > 1) printf("Potential offset: will use %d sampling points per slice, total nz=%d (%d)\n",nzPerSlice,nz,nzPerSlice >> 1);

		dkx = 0.5*OVERSAMP_X/(nx*muls->resolutionX);
		dky = 0.5*OVERSAMP_X/(ny*muls->resolutionY);
		dkz = nzPerSlice/(double)(nz*muls->sliceThickness);
		kmax2 = 0.5*nx*dkx/(double)OVERSAMP_X;	// largest k that we'll admit

		// printf("Cutoff scattering angle:kmax=%g, smax=%g (1/A), dk=(%g,%g %g)\n",kmax2,S_SCALE*kmax2,dkx,dky,dkz);
		scatParOffs[0][N_SF-1] = 1.2*kmax2;
		scatParOffs[0][N_SF-2] = 1.1*kmax2;
		scatParOffs[0][N_SF-3] = kmax2;
		if (scatParOffs[0][N_SF-4] > scatParOffs[0][N_SF-3]) {

			if (1) {
				// set additional scattering parameters to zero:
				for (ix = 0;ix < N_SF-10;ix++) {
					if (scatParOffs[0][N_SF-4-ix] < scatParOffs[0][N_SF-3]-0.001*(ix+1)) break;
					scatParOffs[0][N_SF-4-ix] = scatParOffs[0][N_SF-3]-0.001*(ix+1);
					for (iy=1; iy<N_ELEM;iy++) scatParOffs[iy][N_SF-4-ix] = 0;
				}
			}
			else {
				for (ix = 0;ix < 20;ix++) {
					if (scatParOffs[0][N_SF-4-ix] < scatParOffs[0][N_SF-3]) break;
					scatParOffs[0][N_SF-4-ix] = scatParOffs[0][N_SF-3];
					for (iy=1; iy<N_ELEM;iy++) scatParOffs[iy][N_SF-4-ix] = scatParOffs[iy][N_SF-3];
				}
			}
			if (muls->printLevel > 1) printf("getAtomPotentialOffset3D: reduced angular range of scattering factor to %g/A!\n",
				scatParOffs[0][N_SF-4-ix]);
		}  // end of if (scatParOffs[0][N_SF-4] > scatParOffs[0][N_SF-3])
		kmax2 *= kmax2;

		atPot = (fftwf_complex **)malloc((NZMAX+1)*sizeof(fftwf_complex *));
		for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
		temp  = (fftwf_complex*)fftwf_malloc(nx*nz*sizeof(fftwf_complex));
	}
	// initialize this atom, if it has not been done yet:
	if (atPot[Znum] == NULL) {
#if USE_REZ_SFACTS
		iKind = Znum;
#else
		printf("Using charged atoms only works with scattering factors by Rez et al!\n",Znum);
		exit(0);
#endif


		// setup cubic spline interpolation:
		splinh(scatParOffs[0],scatParOffs[iKind],splinb,splinc,splind,N_SF);

		atPot[Znum] = (fftwf_complex*)fftwf_malloc(nx*nz/4*sizeof(fftwf_complex));
		memset(temp,0,nx*nz*sizeof(fftwf_complex));
		kzmax	 = dkz*nz/2.0;
		// define x-and z-position of atom center:
		// The atom should sit in the top-left corner,
		// however (nzPerSlice+1)/2 above zero in z-direction
		xPos = -2.0*PI*0.0;  // or muls->resolutionX*nx/(OVERSAMP_X), if in center
		izOffset = (nzPerSlice-1)/2;
		zPos = -2.0*PI*(muls->sliceThickness/nzPerSlice*(izOffset));

		// kzborder = dkz*(nz/(2*OVERSAMP_Z) -1);
		for (iz=0;iz<nz;iz++) {
			kz = dkz*(iz<nz/2 ? iz : iz-nz);
			// We also need to taper off the potential in z-direction
			// in order to avoid cutoff artifacts.
			// zScale = fabs(kz) <= kzborder ? 1.0 : 0.5+0.5*cos(M_PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
			// printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
			for (ix=0;ix<nx;ix++) {
				kx = dkx*(ix<nx/2 ? ix : ix-nx);
				s2 = (kx*kx+kz*kz);
				// if this is within the allowed circle:
				if (s2<kmax2) {
					ind3d = ix+iz*nx;
					// f = fe3D(Znum,k2,muls->tds,1.0,muls->scatFactor);
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					f = seval(scatParOffs[0],scatParOffs[iKind],splinb,splinc,splind,N_SF,sqrt(s2))*exp(-s2*B*0.25);
					// perform the qy-integration for qy <> 0:
					for (iy=1;iy<nx;iy++) {
						s3 = dky*iy;
						s3 = s3*s3+s2;
						if (s3<kmax2) {
							f += 2*seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s3))*exp(-s3*B*0.25);
						}
						else break;
					}
					f *= dkx;
					// note that the factor 2 is missing in the phase (2pi k*r)
					// this places the atoms in the center of the box.
					phase  = kx*xPos + kz*zPos;
					temp[ind3d][0] = f*cos(phase);	// *zScale
					temp[ind3d][1] = f*sin(phase);	// *zScale
					// if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
				}
			}
		} // for iz ...




#if SHOW_SINGLE_POTENTIAL
		imageio = ImageIOPtr(new CImageIO(nz, nx, 0, dkx, dkz, std::vector<double>(),
			"rec. space potential"));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		imageio->SetThickness(muls->sliceThickness);
		imageio->WriteComplexImage((void**)temp, fileName);
#endif

		plan = fftwf_plan_dft_2d(nz,nx,temp,temp,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
		// dx2 = muls->resolutionX*muls->resolutionX/(OVERSAMP_X*OVERSAMP_X);
		// dy2 = muls->resolutionY*muls->resolutionY/(OVERSAMP_X*OVERSAMP_X);
		// dz2 = muls->sliceThickness*muls->sliceThickness/(OVERSAMP_Z*OVERSAMP_Z);
		// We also make sure that the potential touches zero at least somewhere.  This will avoid
		// sharp edges that could produce ringing artifacts.
		// It is certainly debatable whether this is a good apprach, or not.
		// printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
		// min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
		for (ix=0;ix<nx/2;ix++)  for (iz=0;iz<nz/2;iz++) {
			ind3d = ix+iz*nx/2;
			// Integrate over nzPerSlice neighboring layers here:::::::::
			for (zScale=0,iiz=-izOffset;iiz<=izOffset;iiz++) {
				if (iz+izOffset+iiz < nz/2) zScale += temp[ix+(iz+izOffset+iiz)*nx][0];
			}
			if (zScale < 0) zScale = 0;
			// assign the iz-th slice the sum of the 3 other slices:
			// and divide by unit cell volume (since this is in 3D):
			// Where does the '1/2' come from???  OVERSAMP_X*OVERSAMP_Y/8 = 1/2
			// if nothing has changed, then OVERSAMP_X=2 OVERSAMP_Z=18.
			// remember, that s=0.5*k;
			// This potential will later again be scaled by lambda*gamma (=0.025*1.39139)
			atPot[Znum][ind3d][0] = 47.8658*dkx*dkz/(nz)*zScale;
			atPot[Znum][ind3d][1] = 0;
		}
#if SHOW_SINGLE_POTENTIAL
		imageio = ImageIOPtr(new CImageIO(nz/2, nx/2, 0, muls->resolutionX/OVERSAMP_X,
		muls->sliceThickness/nzPerSlice));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		imageio->SetThickness(nz*muls->sliceThickness/nzPerSlice);
		sprintf(fileName,"potentialOffs_rz_%d.img",Znum);
		ptr = atPot[Znum];
		imageio->WriteComplexImage((void**)ptr, fileName);
#endif
		if (muls->printLevel > 1) printf("Created 3D (r-z) %d x %d potential offset array for Z=%d (%d, B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)\n",
			nx/2,nz/2,Znum,iKind,B,dkx,dky,dkz,izOffset);
	}
	*Nz_lut = nz/2;
	*nzSub = nzPerSlice;
	*Nr	 = nx/2;
	return atPot[Znum];
}
// #undef SHOW_SINGLE_POTENTIAL
// end of fftwf_complex *getAtomPotential3D(...)

#define PHI_SCALE 47.87658
// #define SHOW_SINGLE_POTENTIAL 0
////////////////////////////////////////////////////////////////////////////
// This function should be used yet, because it computes the projected
// potential wrongly, since it doe not yet perform the projection!!!
fftwf_complex *getAtomPotential2D(int Znum, MULS *muls,double B) {
	int ix,iy,iz,ind,iKind;
	double min;
	fftwf_plan plan;
	static double f,phase,s2,s3,kmax2,kx,ky,dkx,dky;
	static int nx,ny;
	static fftwf_complex **atPot = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
	ImageIOPtr imageio = ImageIOPtr();
	static char fileName[256];
#endif
	static double *splinb=NULL;
	static double *splinc=NULL;
	static double *splind=NULL;


	// scattering factors in:
	// float scatPar[4][30]
	if (atPot == NULL) {
		splinb = double1D(N_SF, "splinb" );
		splinc = double1D(N_SF, "splinc" );
		splind = double1D(N_SF, "splind" );

		nx = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX);
		ny = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY);
		dkx = 0.5*OVERSAMP_X/((nx)*muls->resolutionX);
		dky = 0.5*OVERSAMP_X/((ny)*muls->resolutionY);
		kmax2 = 0.5*nx*dkx/(double)OVERSAMP_X;  // largest k that we'll admit

		printf("Cutoff scattering angle:kmax=%g (1/A)\n",kmax2);
		scatPar[0][N_SF-1] = 1.2*kmax2;
		scatPar[0][N_SF-2] = 1.1*kmax2;
		scatPar[0][N_SF-3] = kmax2;
		if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3]) {

			if (1) {
				// set additional scattering parameters to zero:
				for (ix = 0;ix < N_SF-10;ix++) {
					if (scatPar[0][N_SF-4-ix] < scatPar[0][N_SF-3]-0.001*(ix+1)) break;
					scatPar[0][N_SF-4-ix] = scatPar[0][N_SF-3]-0.001*(ix+1);
					for (iy=1; iy<N_ELEM;iy++) scatPar[iy][N_SF-4-ix] = 0;
				}
			}

			if (muls->printLevel > 1) printf("getAtomPotential2D: reduced angular range of scattering factor to %g/A!\n",
				scatPar[0][N_SF-4-ix]);
		}  // end of if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3])
		kmax2 *= kmax2;



		atPot = (fftwf_complex **)malloc((NZMAX+1)*sizeof(fftwf_complex *));
		for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
	}
	// initialize this atom, if it has not been done yet:
	if (atPot[Znum] == NULL) {
		iKind = Znum;
		// setup cubic spline interpolation:
		splinh(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF);

		atPot[Znum] = (fftwf_complex*) fftwf_malloc(nx*ny*sizeof(fftwf_complex));
		// memset(temp,0,nx*nz*sizeof(fftwf_complex));
		memset(atPot[Znum],0,nx*ny*sizeof(fftwf_complex));
		for (ix=0;ix<nx;ix++) {
			kx = dkx*(ix<nx/2 ? ix : nx-ix);
			for (iy=0;iy<ny;iy++) {
				ky = dky*(iy<ny/2 ? iy : ny-iy);
				s2 = (kx*kx+ky*ky);
				// if this is within the allowed circle:
				if (s2<kmax2) {
					ind = iy+ix*ny;
					// f = fe3D(Znum,k2,muls->tds,1.0,muls->scatFactor);
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					f = seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s2))*exp(-s2*B*0.25);
					phase = PI*(kx*muls->resolutionX*nx+ky*muls->resolutionY*ny);
					atPot[Znum][ind][0] = f*cos(phase);
					atPot[Znum][ind][1] = f*sin(phase);
				}
			}
		}
#if SHOW_SINGLE_POTENTIAL == 1
		imageio = ImageIOPtr(new CImageIO(ny, nx, 0, dkx, dky, std::vector<double>(),
		"potential"));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		imageio->SetThickness(muls->sliceThickness);
		imageio->WriteComplexImage((void**)atPot[Znum], fileName);
#endif
		plan = fftwf_plan_dft_2d(nx,ny,atPot[Znum],atPot[Znum],FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
		for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
				atPot[Znum][iy+ix*ny][0] *= dkx*dky*(OVERSAMP_X*OVERSAMP_X);
		}
		// make sure we don't produce negative potential:
		// if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny][0] -= min;
#if SHOW_SINGLE_POTENTIAL == 1
		imageio = ImageIOPtr(new CImageIO(nx, ny, 0, muls->resolutionX/OVERSAMP_X,
			muls->resolutionY/OVERSAMP_X, std::vector<double>(), "potential"));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		//imageio->SetThickness(nz*muls->sliceThickness/nzPerSlice);
		sprintf(fileName,"potential_%d.img",Znum);
		imageio->WriteComplexImage((void**)atPot[Znum], fileName);
#endif
		printf("Created 2D %d x %d potential array for Z=%d (%d, B=%g A^2)\n",nx,ny,Znum,iKind,B);
	}
	return atPot[Znum];
}
#undef PHI_SCALE
#undef SHOW_SINGLE_POTENTIAL

void writePix(char *outFile,fftw_complex **pict,MULS *muls,int iz) {
	real *sparam;
	real rmin,rmax;
	int i,j, result;

	rmin  = pict[0][0][0];
	rmax  = rmin;

	sparam = float1D( NPARAM, "sparam" );

	for( i=0; i<(*muls).nx; i++)
		for(j=0; j<(*muls).ny; j++)
		{
			if(pict[i][j][0] < rmin ) rmin = pict[i][j][0];
			if(pict[i][j][0] > rmax ) rmax = pict[i][j][0];
		}
		printf("min: %g  max: %g\n",rmin,rmax);

		if (rmin==rmax)
			rmax = rmin+1.0;
		sparam[pRMAX]  = rmax;
		sparam[pIMAX]  = rmax;
		sparam[pRMIN]  = rmin;
		sparam[pIMIN]  = rmin;
		sparam[pXCTILT]  = 0.0f;
		sparam[pYCTILT] = 0.0f;
		sparam[pENERGY] = (*muls).v0;
		sparam[pDX] = (*muls).ax/(real)(*muls).nx;
		sparam[pDY] = (*muls).by/(real)(*muls).ny;
		sparam[pWAVEL] = (real)wavelength((*muls).v0);
		sparam[pNSLICES] = 0.0F;  /* ??? */
		sparam[pDEFOCUS] = 0.0;
		sparam[pOAPERT] = 0.0;
		sparam[pCS] = 0.0;
		sparam[pCAPERT] = 0.0;
		sparam[pDDF] = 0.0;
		sparam[pC]= (*muls).cz[iz];
		result = 1;
		/*
		result = tcreateFloatPixFile(outFile,pict,(long)(*muls).nx,
		(long)(*muls).ny,1,sparam);
		*/

		if (result != 1)
			printf("\ncould not write output file %s\n",outFile);
}


/**********************************************
* This function creates a incident STEM probe
* at position (dx,dy)
* with parameters given in muls
*
* The following Abberation functions are being used:
* 1) ddf = Cc*dE/E + Cc2*(dE/E)^2,
*    Cc, Cc2 = chrom. Abber. (1st, 2nd order) [1]
* 2) chi(qx,qy) = (2*pi/lambda)*{0.5*C1*(qx^2+qy^2)+
*                 0.5*C12a*(qx^2-qy^2)+
*                 C12b*qx*qy+
*                 C21a/3*qx*(qx^2+qy^2)+
*                 ...
*                 +0.5*C3*(qx^2+qy^2)^2
*                 +0.125*C5*(qx^2+qy^2)^3
*                 ... (need to finish)
*
*
*    qx = acos(kx/K), qy = acos(ky/K)
*
* References:
* [1] J. Zach, M. Haider,
*    "Correction of spherical and Chromatic Abberation
*     in a low Voltage SEM", Optik 98 (3), 112-118 (1995)
* [2] O.L. Krivanek, N. Delby, A.R. Lupini,
*    "Towards sub-Angstroem Electron Beams",
*    Ultramicroscopy 78, 1-11 (1999)
*
*********************************************/

#define SMOOTH_EDGE 5 // make a smooth edge on AIS aperture over +/-SMOOTH_EDGE pixels

void probeShiftAndCrop(MULS *muls, WavePtr wave, double dx, double dy, double cnx, double cny)
{
	// Robert A. McLeod
	// 09 April 2014
	// This function takes an input wave and shifts it by [dx,dy] and crops it to [cnx,cny]
	// In general the input wave should have wave.nx > nx + 2*abs(dx) and wave.ny > ny + 2*abs(dy) for a scanned system

	//
	printf("Debug: probeShiftAndCrop does nothing at present\n");
}

void probe(MULS *muls, WavePtr wave, double dx, double dy)
{
	// static char *plotFile = "probePlot.dat",systStr[32];
	int ix, iy, nx, ny, ixmid, iymid;
	int CsDefAstOnly = 0;
	float rmin, rmax, aimin, aimax;
	// float **pixr, **pixi;
	double  kx, ky, ky2,k2, ktheta2, ktheta, k2max, v0, wavlen,ax,by,x,y,
		rx2, ry2,rx,ry, pi, scale, pixel,alpha,
		df, df_eff, chi1, chi2,chi3, sum, chi, time,r,phi;
	double gaussScale = 0.05;
	double envelope,delta,avgRes,edge;

	// FILE *fp=NULL;

	/* temporary fix, necessary, because fftw has rec. space zero
	in center of image:
	*/
	nx = (*muls).nx;
	ny = (*muls).ny;
	ax = nx*(*muls).resolutionX;
	by = ny*(*muls).resolutionY;
	dx = ax-dx;
	dy = by-dy;
	gaussScale = (*muls).gaussScale;
	// average resolution:
	avgRes = sqrt(0.5*(muls->resolutionX*muls->resolutionX+muls->resolutionY*muls->resolutionY));
	edge = SMOOTH_EDGE*avgRes;

	/********************************************************
	* formulas from:
	* http://cimesg1.epfl.ch/CIOL/asu94/ICT_8.html
	*
	* dE_E = dE/E = energy spread of emitted electrons
	* dV_V = dV/V = acc. voltage fluctuations
	* dI_I = dI/I = lens current fluctuations
	* delta defocus in Angstroem (Cc in A)
	*******************************************************/
	delta = muls->Cc*muls->dE_E;
	if (muls->printLevel > 2) printf("defocus offset: %g nm (Cc = %g)\n",delta,muls->Cc);

	if (wave->wave == NULL) {
		printf("Error in probe(): Wave not allocated!\n");
		exit(0);
	}

	/**********************************************************
	*  Calculate misc constants
	*********************************************************/
	time = cputim( );
	pi = 4.0 * atan( 1.0 );

	rx = 1.0/ax;
	rx2 = rx * rx;
	ry = 1.0/by;
	ry2 = ry * ry;

	ixmid = nx/2;
	iymid = ny/2;

	// df = muls->df0;
	v0 = muls->v0;
	wavlen = 12.26/ sqrt( v0*1.e3 + v0*v0*0.9788 );

	/*  printf("Wavelength: %g A\n",wavlen);
	*/


	// chi2 = (*muls).Cs*0.5*wavlen*wavlen;
	// chi3 = (*muls).C5*0.25*wavlen*wavlen*wavlen*wavlen;
	/* delta *= 0.5*delta*pi*pi*wavlen*wavlen; */

	/* convert convergence angle from mrad to rad */
	alpha = 0.001*muls->alpha;
	k2max = sin(alpha)/wavlen;  /* = K0*sin(alpha) */
	k2max = k2max * k2max;

	/*   Calculate MTF
	NOTE zero freg is in the bottom left corner and
	expandes into all other corners - not in the center
	this is required for FFT

	PIXEL = diagonal width of pixel squared
	if a pixel is on the apertur boundary give it a weight
	of 1/2 otherwise 1 or 0
	*/
	pixel = ( rx2 + ry2 );
	scale = 1.0/sqrt((double)nx*(double)ny);

	/*
	if ((muls.a33 == 0) && (muls.a31 == 0) && (muls.a44 == 0) && (muls.a42 == 0) &&
	(muls.a55 == 0) && (muls.a53 == 0) && (muls.a51 == 0) &&
	(muls.a66 == 0) && (muls.a64 == 0) && (muls.a62 == 0) && (muls.C5 == 0)) {
	CsDefAstOnly = 1;
	}
	*/

	for( iy=0; iy<ny; iy++) {
		ky = (double) iy;
		if( iy > iymid ) ky = (double) (iy-ny);
		ky2 = ky*ky*ry2;
		for( ix=0; ix<nx; ix++) {
			kx = (double) ix;
			if( ix > ixmid ) kx = (double) (ix-nx);
			k2 = kx*kx*rx2 + ky2;
			ktheta2 = k2*(wavlen*wavlen);
			ktheta = sqrt(ktheta2);
			phi = atan2(ry*ky,rx*kx);
			// compute the effective defocus from the actual defocus and the astigmatism:
			// df_eff = df + muls->astigMag*cos(muls->astigAngle+phi);

			// chi = chi1*k2*(df_eff +chi2*k2)-2.0*pi*( (dx*kx/ax) + (dy*ky/by) );
			// defocus, astigmatism, and shift:
			chi = ktheta2*(muls->df0+delta + muls->astigMag*cos(2.0*(phi-muls->astigAngle)))/2.0;
			ktheta2 *= ktheta;  // ktheta^3
			if ((muls->a33 > 0) || (muls->a31 > 0)) {
				chi += ktheta2*(muls->a33*cos(3.0*(phi-muls->phi33))+muls->a31*cos(phi-muls->phi31))/3.0;
			}
			ktheta2 *= ktheta;   // ktheta^4
			if ((muls->a44 > 0) || (muls->a42 > 0) || (muls->Cs != 0)) {
				// chi += ktheta2*(muls->a33*cos(3*(phi-muls->phi33))+muls->a31*cos(phi-muls->phi31))/3.0;
				chi += ktheta2*(muls->a44*cos(4.0*(phi-muls->phi44))+muls->a42*cos(2.0*(phi-muls->phi42))+muls->Cs)/4.0;
				//                     1/4*(a(4,4).*cos(4*(kphi-phi(4,4)))+a(4,2).*cos(2*(kphi-phi(4,2)))+c(4)).*ktheta.^4+...
			}
			ktheta2 *= ktheta;    // ktheta^5
			if ((muls->a55 > 0) || (muls->a53 > 0) || (muls->a51 > 0)) {
				chi += ktheta2*(muls->a55*cos(5.0*(phi-muls->phi55))+muls->a53*cos(3.0*(phi-muls->phi53))+muls->a51*cos(phi-muls->phi51))/5.0;
				//                     1/5*(a(5,5).*cos(5*(kphi-phi(5,5)))+a(5,3).*cos(3*(kphi-phi(5,3)))+a(5,1).*cos(1*(kphi-phi(5,1)))).*ktheta.^5+...
			}
			ktheta2 *= ktheta;    // ktheta^6
			if ((muls->a66 > 0) || (muls->a64 > 0) || (muls->a62 = 0) || (muls->C5 != 0)) {
				chi += ktheta2*(muls->a66*cos(6.0*(phi-muls->phi66))+muls->a64*cos(4.0*(phi-muls->phi64))+muls->a62*cos(2.0*(phi-muls->phi62))+muls->C5)/6.0;
				//                     1/6*(a(6,6).*cos(6*(kphi-phi(6,6)))+a(6,4).*cos(4*(kphi-phi(6,4)))+a(6,2).*cos(2*(kphi-phi(6,2)))+c(6)).*ktheta.^6);
			}

			chi *= 2*pi/wavlen;
			chi -= 2.0*pi*( (dx*kx/ax) + (dy*ky/by) );
			// include higher order aberrations


			if ( ( (*muls).ismoth != 0) &&
				( fabs(k2-k2max) <= pixel)) {
					wave->wave[ix][iy][0]= (float) ( 0.5*scale * cos(chi));
					wave->wave[ix][iy][1]= (float) (-0.5*scale* sin(chi));
			}
			else if ( k2 <= k2max ) {
				wave->wave[ix][iy][0]= (float)  scale * cos(chi);
				wave->wave[ix][iy][1]= (float) -scale * sin(chi);
			}
			else {
				wave->wave[ix][iy][0] = wave->wave[ix][iy][1] = 0.0f;
			}
		}
	}
	/* Fourier transform into real space */
	// fftwnd_one(muls->fftPlanInv, &(muls->wave[0][0]), NULL);
#if FLOAT_PRECISION == 1
	fftwf_execute(wave->fftPlanWaveInv);
#else
	fftw_execute(wave->fftPlanWaveInv);
#endif
	/**********************************************************
	* display cross section of probe intensity
	*/

	/* multiply with gaussian in Real Space in order to avoid artifacts */
	if (muls->gaussFlag) {
		for( ix=0; ix<nx; ix++) {
			for( iy=0; iy<ny; iy++) {
				r = exp(-((ix-nx/2)*(ix-nx/2)+(iy-ny/2)*(iy-ny/2))/(nx*nx*gaussScale));
				wave->wave[ix][iy][0] *= (float)r;
				wave->wave[ix][iy][1] *= (float)r;
			}
		}
	}

	/* Apply AIS aperture in Real Space */
	// printf("center: %g,%g\n",dx,dy);
	if (muls->aAIS > 0) {
		for( ix=0; ix<nx; ix++) {
			for( iy=0; iy<ny; iy++) {
				x = ix*muls->resolutionX-dx;
				y = iy*muls->resolutionY-dy;
				r = sqrt(x*x+y*y);
				delta = r-0.5*muls->aAIS+edge;
				if (delta > 0) {
					wave->wave[ix][iy][0] = 0;
					wave->wave[ix][iy][1] = 0;
				}
				else if (delta >= -edge) {
					scale = 0.5*(1-cos(pi*delta/edge));
					wave->wave[ix][iy][0] = scale*wave->wave[ix][iy][0];
					wave->wave[ix][iy][1] = scale*wave->wave[ix][iy][1];
				}
			}
		}
	}

	/*  Normalize probe intensity to unity  */

	sum = 0.0;
	for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
		sum +=  wave->wave[ix][iy][0]*wave->wave[ix][iy][0]
	+ wave->wave[ix][iy][1]*wave->wave[ix][iy][1];

	scale = 1.0 / sum;
	scale = scale * ((double)nx) * ((double)ny);
	scale = (double) sqrt( scale );

	for( ix=0; ix<nx; ix++)
		for( iy=0; iy<ny; iy++) {
			wave->wave[ix][iy][0] *= (float) scale;
			wave->wave[ix][iy][1] *= (float) scale;
		}

		/*  Output results and find min and max to echo
		remember that complex pix are stored in the file in FORTRAN
		order for compatability
		*/

		rmin = wave->wave[0][0][0];
		rmax = rmin;
		aimin = wave->wave[0][0][1];
		aimax = aimin;
		for( iy=0; iy<ny; iy++) {
			for( ix=0; ix<nx; ix++) {
				if( wave->wave[ix][iy][0] < rmin ) rmin = wave->wave[ix][iy][0];
				if( wave->wave[ix][iy][0] > rmax ) rmax = wave->wave[ix][iy][0];
				if( wave->wave[ix][iy][1] < aimin ) aimin = wave->wave[ix][iy][1];
				if( wave->wave[ix][iy][1] > aimax ) aimax = wave->wave[ix][iy][1];
			}
		}
		(*muls).rmin = rmin;
		(*muls).rmax = rmax;
		(*muls).aimin = aimin;
		(*muls).aimax = aimax;

		/**********************************************************/

}  /* end probe() */

/**************************************************************
* The imaginary part of the trans arrays is already allocated
* The projected potential is already located in trans[][][][0]
*
**************************************************************/
#define PHI_SCALE 47.87658
void initSTEMSlices(MULS *muls, int nlayer) {
	int printFlag = 0;
	int ilayer;
	int nbeams;
	double scale,vz,vzscale,mm0,wavlen;
	int nx,ny,ix,iy; // iz;
	real temp,k2max,k2,kx,ky;
	static real *kx2= NULL,*ky2 = NULL; /* *kx= NULL,*ky= NULL, */
	real pi;
	double fftScale;
	double timer1,timer2,time2=0,time1=0;
	// char filename[32];

	pi = (float)PI;

	nx = (*muls).potNx;
	ny = (*muls).potNy;

	/**************************************************/
	/* Setup all the reciprocal lattice vector arrays */
	if ((kx2 == NULL)||(ky2 == NULL)) {
		kx2    = float1D(nx, "kx2" );
		ky2    = float1D(ny, "ky2" );
		/*
		kx     = float1D(nx, "kx" );
		ky     = float1D(ny, "ky" );
		*/
		for(ix=0; ix<nx; ix++) {
			kx = (ix>nx/2) ? (real)(ix-nx)/(*muls).potSizeX :
				(real)ix/(*muls).potSizeX;
		kx2[ix] = kx*kx;
		/*
		kx[ix] = (ix>nx/2) ? (real)(ix-nx)/(*muls).potSizeX :
		(real)ix/(*muls).potSizeX;
		kx2[ix] = kx[ix]*kx[ix];
		*/
		}
		for( iy=0; iy<ny; iy++) {
			ky = (iy>ny/2) ?
				(real)(iy-ny)/(*muls).by : (real)iy/(*muls).potSizeY;
			ky2[iy] = ky*ky;
		}
		if (muls->printLevel > 2) printf("Reciprocal lattice vector arrays initialized ... \n");
	}
	/**************************************************/
	/* setup of all the necessary parameters: */
	wavlen = (real) wavelength((*muls).v0);

	nbeams = 0;

	/* setup of k-vector arrays, etc... */
	k2max = nx/(2.0F*(*muls).potSizeX);
	temp = ny/(2.0F*(*muls).potSizeY);
	if( temp < k2max ) k2max = temp;
	k2max = (real)(BW * k2max);
	if (printFlag) {
		printf("Bandwidth limited to a real space resolution of %f Angstroms\n",
			1.0F/k2max);
		printf("   (= %.2f mrad)  for symmetrical anti-aliasing.",
			wavlen*k2max*1000.0F);
		printf(" (a: %g, b: %g)\n",(*muls).potSizeX,(*muls).potSizeY);
	}
	k2max = k2max*k2max;
	(*muls).k2max = k2max;

	if(muls->trans==NULL) {
		printf("Memory for trans has not been allocated\n");
		exit(0);
	}


	/**************************************************************
	* read in the potential layers and exponentiate them
	* to form the transmission function exp(-i*phi(x)) = transr/i
	*
	* This is independent of kx, ky - so we can load these arrays once,
	* and keep them for all inc. beam angles
	*
	**************************************************************/


	mm0 = 1.0F + muls->v0/511.0F;  // relativistic corr. factor gamma
	// scale = mm0*sigma(muls->v0) / 1000.0;  /* in 1/(volt-Angstroms) */
	// we will need to increase lambda and decrease gamm by the mean inner crystal potential effect.
	scale = mm0*wavelength(muls->v0);
	if (muls->printLevel > 1) printf("Making phase gratings for %d layers (scale=%g rad/VA, gamma=%g, sigma=%g) ... \n",nlayer,scale,mm0,sigma(muls->v0));


	fftScale = 1.0/(nx*ny);
	vzscale= 1.0;
	timer1 = cputim();
	for( ilayer=0;  ilayer<nlayer; ilayer++ ) {
		timer2 = cputim();
		for( iy=0; iy<ny; iy++) for( ix=0; ix<nx; ix++) {
			vz= muls->trans[ilayer][ix][iy][0]*scale;  // scale = lambda*gamma
			// include absorption:
			// vzscale= exp(-(*muls).trans[ilayer][ix][iy][1]*scale);
			/* printf("vz(%d %d) = %g\n",ix,iy,vz); */
			muls->trans[ilayer][ix][iy][0] =  cos(vz);
			muls->trans[ilayer][ix][iy][1] =  sin(vz);
		}
	}

	/*******************************************************************
	* FFT/IFFT the transmit functions in order to bandwidth limit them
	*******************************************************************/
	if (muls->bandlimittrans) {
		timer2 = cputim();
#if FLOAT_PRECISION == 1
		fftwf_execute(muls->fftPlanPotForw);
#else
		fftw_execute(muls->fftPlanPotForw);
#endif
		time2 = cputim()-timer2;
		//     printf("%g sec used for 1st set of FFTs\n",time2);
		for( ilayer=0;  ilayer<nlayer; ilayer++ ) {
			for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
				k2= ky2[iy] + kx2[ix];
				if (k2 < k2max) {
					nbeams++;
					(*muls).trans[ilayer][ix][iy][0] *= fftScale;
					(*muls).trans[ilayer][ix][iy][1] *= fftScale;
				}
				else {
					(*muls).trans[ilayer][ix][iy][0] = 0.0F;
					(*muls).trans[ilayer][ix][iy][1] = 0.0F;
				}
			}
		}  /* end for(ilayer=... */
		timer2 = cputim();
		// old code: fftwnd_one((*muls).fftPlanPotInv, (*muls).trans[ilayer][0], NULL);
#if FLOAT_PRECISION == 1
		fftwf_execute(muls->fftPlanPotInv);
#else
		fftw_execute(muls->fftPlanPotInv);
#endif
		time2 += cputim()-timer2;
	}  /* end of ... if bandlimittrans */
	time1 = cputim()-timer1;

	if (muls->printLevel > 1) {
		if ((*muls).bandlimittrans) {
			printf("%g sec used for making phase grating, %g sec for %d FFTs.\n",
				time1,time2,2*nlayer);
		}
		else
			printf("%g sec used for making phase grating, no bandwidth limiting\n",
			time1);

		if (printFlag)
			printf("Number of symmetrical non-aliasing beams = %d\n", nbeams);
	}

	/*
	printf("Size in pixels Nx x Ny= %d x %d = %d beams\n",
	nx,ny, nx*ny);
	printf("Lattice constant a = %.4f, b = %.4f\n", (*muls).ax,(*muls).by);
	*/
}  // initSTEMSlices

#undef PHI_SCALE




/******************************************************************
* runMulsSTEM() - do the multislice propagation in STEM/CBED mode
*
*    Each probe position is running this function.  Each CPU is thus
*      running a separate instance of the function.  It is nested in
*      the main OpenMP parallel region - specifying critical, single, and
*      barrier OpenMP pragmas should be OK.
*
* waver, wavei are expected to contain incident wave function
* they will be updated at return
*****************************************************************/
int runMulsSTEM(MULS *muls, WavePtr wave) {
	int printFlag = 0;
	int showEverySlice=1;
	int islice,i,ix,iy,mRepeat;
	real cztot=0.0;
	real wavlen,scale,sum=0.0; //,zsum=0.0
	// static int *layer=NULL;
	real x,y;
	int absolute_slice;

	char outStr[64];
	double fftScale;

	printFlag = (muls->printLevel > 3);
	fftScale = 1.0/(muls->nx*muls->ny);

	wavlen = (real)wavelength((*muls).v0);

	/*  calculate the total specimen thickness and echo */
	cztot=0.0;
	for( islice=0; islice<(*muls).slices; islice++) {
		cztot += (*muls).cz[islice];
	}
	if (printFlag)
		printf("Specimen thickness: %g Angstroms\n", cztot);

	scale = 1.0F / (((real)muls->nx) * ((real)muls->ny));

	for (mRepeat = 0; mRepeat < muls->mulsRepeat1; mRepeat++)
	{
		for( islice=0; islice < muls->slices; islice++ )
		{
			absolute_slice = (muls->totalSliceCount+islice);

			// if ((muls->cubez > 0) && (muls->thickness >= muls->cubez)) break;
			//  else if ((muls->cubez == 0) && (muls->thickness >= muls->c)) break;

			/***********************************************************************
			* Transmit is a simple multiplication of wave with trans in real space
			**********************************************************************/
			transmit((void **)wave->wave, (void **)(muls->trans[islice]), muls->nx,muls->ny, wave->iPosX, wave->iPosY);

			//    writeImage_old(wave,(*muls).nx,(*muls).ny,(*muls).thickness,"wavet.img");
			/*****************************************************
			* remember: prop must be here to anti-alias
			* propagate is a simple multiplication of wave with prop
			* but it also takes care of the bandwidth limiting
			*******************************************************/
#if FLOAT_PRECISION == 1
			fftwf_execute(wave->fftPlanWaveForw);
#else
			fftw_execute(wave->fftPlanWaveForw);
#endif
			propagate_slow((void **)wave->wave, muls->nx, muls->ny, muls);

			collectIntensity(muls, wave, muls->totalSliceCount+islice*(1+mRepeat));

			//if (muls->mode != STEM) {
				/* write pendelloesung plots, if this is not STEM */
			//	writeBeams(muls,wave,islice, absolute_slice);
			//}

			// go back to real space:
#if FLOAT_PRECISION == 1
			fftwf_execute(wave->fftPlanWaveInv);
#else
			fftw_execute(wave->fftPlanWaveInv);
#endif
			// old code: fftwnd_one((*muls).fftPlanInv,(fftw_complex *)wave[0][0], NULL);
			fft_normalize((void **)wave->wave,muls->nx,muls->ny);

			/*
			sprintf(outStr,"wave%d.img",islice);
			writeImage_old(wave,(*muls).nx,(*muls).ny,(*muls).thickness,"wavep.img");
			*/

			// write the intermediate TEM wave function:

			/********************************************************************
			* show progress:
			********************************************************************/
			wave->thickness = (absolute_slice+1)*muls->sliceThickness;
			if ((printFlag)) {
				sum = 0.0;
				for( ix=0; ix<(*muls).nx; ix++)  for( iy=0; iy<(*muls).ny; iy++) {
					sum +=  wave->wave[ix][iy][0]* wave->wave[ix][iy][0] +
						wave->wave[ix][iy][1]* wave->wave[ix][iy][1];
				}
				sum *= scale;

				sprintf(outStr,"position (%3d, %3d), slice %4d (%.2f), int. = %f",
					wave->detPosX, wave->detPosY,
					muls->totalSliceCount+islice,wave->thickness,sum );
				if (showEverySlice)
					printf("%s\n",outStr);
				else {
					printf("%s",outStr);
					for (i=0;i<(int)strlen(outStr);i++) printf("\b");
				}
			}

			if ( (muls->mode == TEM) || ((muls->mode == CBED) && (muls->saveLevel > 1)) || ((muls->mode == NBED) && (muls->saveLevel > 1)) )
			{
				// TODO (MCS 2013/04): this restructure probably broke this file saving -
				//   need to rewrite a function to save things for TEM/CBED?
				// This used to call interimWave(muls,wave,muls->totalSliceCount+islice*(1+mRepeat));
				// interimWave(muls,wave,absolute_slice*(1+mRepeat));
				collectIntensity(muls,wave,absolute_slice*(1+mRepeat));
			}
		} /* end for(islice...) */
		// collect intensity at the final slice
		//collectIntensity(muls, wave, muls->totalSliceCount+muls->slices*(1+mRepeat));
	} /* end of mRepeat = 0 ... */
	if (printFlag) printf("\n***************************************\n");

	/****************************************************
	****************************************************
	**           this -WAS- the big loop              **
	****************************************************
	***************************************************/

	// TODO: modifying shared value from multiple threads?
	//#pragma omp single
	muls->rmin  = wave->wave[0][0][0];
	//#pragma omp single
	muls->rmax  = (*muls).rmin;
	//#pragma omp single
	muls->aimin = wave->wave[0][0][1];
	//#pragma omp single
	muls->aimax = (*muls).aimin;

	sum = 0.0;
	for( ix=0; ix<muls->nx; ix++)  for( iy=0; iy<muls->ny; iy++) {
		x =  wave->wave[ix][iy][0];
		y =  wave->wave[ix][iy][1];
		if( x < (*muls).rmin ) (*muls).rmin = x;
		if( x > (*muls).rmax ) (*muls).rmax = x;
		if( y < (*muls).aimin ) (*muls).aimin = y;
		if( y > (*muls).aimax ) (*muls).aimax = y;
		sum += x*x+y*y;
	}
	// TODO: modifying shared value from multiple threads?
	//  Is this sum supposed to be across multiple pixels?
	//#pragma omp critical
	wave->intIntensity = sum*scale;

	if (printFlag) {
		printf( "pix range %g to %g real,\n"
			"          %g to %g imag\n",
			(*muls).rmin,(*muls).rmax,(*muls).aimin,(*muls).aimax);

	}
	if (muls->saveFlag) {
		if ((muls->saveLevel > 1) || (muls->cellDiv > 1)) {
			wave->WriteWave(wave->fileout);
			if (printFlag)
				printf("Created complex image file %s\n",(*wave).fileout);
		}
	}
	return 0;
}  // end of runMulsSTEM


////////////////////////////////////////////////////////////////
// save the current wave function at this intermediate thickness:
void interimWave(MULS *muls,WavePtr wave,int slice) {
	int t;
	char fileName[256];
	std::vector<double> params(9);

	if ((slice < muls->slices*muls->cellDiv-1) && ((slice+1) % muls->outputInterval != 0)) return;

	t = (int)((slice)/muls->outputInterval);
	// write the high tension, too:

	params[0] = muls->v0;  				// high voltage
	params[1] = muls->Cs;				// spherical aberration
	params[2] = muls->df0;				// defocus
	params[3] = muls->astigMag;			// astigmatism
	params[4] = muls->astigAngle;
	params[5] = muls->Cc * sqrt(muls->dE_E*muls->dE_E+muls->dV_V*muls->dV_V+muls->dI_I*muls->dI_I);	// focal spread
	params[6] = muls->alpha;			// illumination convergence angle
	params[7] = muls->btiltx;			// beam tilt in mrad
	params[8] = muls->btilty;			// beam tilt in mrad
	// printf("###  Cc = %f, dE_E = %f, Delta = %f ###\n",muls->Cc,muls->dV_V,muls->Cc * muls->dV_V);

	// produce the following filename:
	// wave_avgCount_thicknessIndex.img or
	// wave_thicknessIndex.img if tds is turned off
	if (muls->tds) sprintf(fileName,"%s/wave_%d_%d.img",muls->folder,muls->avgCount,t);
	else sprintf(fileName,"%s/wave_%d.img",muls->folder,t);
	wave->WriteWave(fileName, "Wave Function", params);
}

/********************************************************************
* collectIntensity(muls, wave, slice)
* collect the STEM signal on the annular detector(s) defined in muls
* and write the appropriate pixel in the image for each detector and thickness
* The number of images is determined by the following formula:
* muls->slices*muls->cellDiv/muls->outputInterval
* There are muls->detectorNum different detectors
*******************************************************************/
void collectIntensity(MULS *muls, WavePtr wave, int slice)
{
	int i,ix,iy,ixs,iys,t;
	real k2;
	double intensity,scale,scaleCBED,scaleDiff,intensity_save;
	char fileName[256],avgName[256];
	float_tt **diffpatAvg = NULL;
	int tCount = 0;

	std::vector<std::vector<DetectorPtr> > detectors;

	scale = muls->electronScale/((double)(muls->nx*muls->ny)*(muls->nx*muls->ny));
	// scaleCBED = 1.0/(scale*sqrt((double)(muls->nx*muls->ny)));
	scaleDiff = 1.0/sqrt((double)(muls->nx*muls->ny));

	tCount = (int)(ceil((double)((muls->slices * muls->cellDiv) / muls->outputInterval)));

	// we write directly to the shared muls object.  This is safe only because
	//    each thread is accessing different pixels in the output images.
	detectors = muls->detectors;

	if (muls->outputInterval == 0) t = 0;
	else if (slice < ((muls->slices*muls->cellDiv)-1))
	{
		t = (int)((slice) / muls->outputInterval);
		//if (t > tCount)
		//{
			// printf("t = %d, which is greater than tCount (%d)\n",t,tCount);
			//t = tCount-1;
		//}
	}
	else
	{
		t = tCount;
	}

	int position_offset = wave->detPosY * muls->scanXN + wave->detPosX;

	// Multiply each image by its number of averages and divide by it later again:
	for (i=0;i<muls->detectorNum;i++)
	{
		detectors[t][i]->image[wave->detPosX][wave->detPosY]  *= detectors[t][i]->Navg;
		detectors[t][i]->image2[wave->detPosX][wave->detPosY] *= detectors[t][i]->Navg;
		detectors[t][i]->error = 0;
	}
	/* add the intensities in the already
	fourier transformed wave function */
	for (ix = 0; ix < muls->nx; ix++)
	{
		for (iy = 0; iy < muls->ny; iy++)
		{
			k2 = muls->kx2[ix]+muls->ky2[iy];
			intensity = (wave->wave[ix][iy][0]*wave->wave[ix][iy][0]+
				wave->wave[ix][iy][1]*wave->wave[ix][iy][1]);
			wave->diffpat[(ix+muls->nx/2)%muls->nx][(iy+muls->ny/2)%muls->ny] = intensity*scaleDiff;
			intensity *= scale;
			for (i=0;i<muls->detectorNum;i++) {
				if ((k2 >= detectors[t][i]->k2Inside) && (k2 <= detectors[t][i]->k2Outside))
				{
					// detector in center of diffraction pattern:
					if ((detectors[t][i]->shiftX == 0) && (detectors[t][i]->shiftY == 0))
					{
						detectors[t][i]->image[wave->detPosX][wave->detPosY] += intensity;
						// misuse the error number for collecting this pixels raw intensity
						detectors[t][i]->error += intensity;
					}
					/* special case for shifted detectors: */
					else
					{
						intensity_save = intensity;
						ixs = (ix+(int)detectors[t][i]->shiftX+muls->nx) % muls->nx;
						iys = (iy+(int)detectors[t][i]->shiftY+muls->ny) % muls->ny;
						intensity = scale * (wave->wave[ixs][iys][0]*wave->wave[ixs][iys][0]+
							wave->wave[ixs][iys][1]*wave->wave[ixs][iys][1]);
						detectors[t][i]->image[wave->detPosX][wave->detPosY] += intensity;
						// repurpose the error number for collecting this pixels raw intensity
						detectors[t][i]->error += intensity;
						/* restore intensity, so that it will not be shifted for the other detectors */
						intensity = intensity_save;
					}
				} /* end of if k2 ... */
			} /* end of for i=0 ... detectorNum */
		} /* end of for iy=0... */
	} /* end of for ix = ... */

	////////////////////////////////////////////////////////////////////////////
	// write the diffraction pattern to disc in case we are working in CBED mode
	if ((muls->mode == CBED) && (muls->saveLevel > 0)) {
		sprintf(avgName,"%s/diff_%d.img",muls->folder,t);
		// for (ix=0;ix<muls->nx*muls->ny;ix++) wave->diffpat[0][ix] *= scaleCBED;
		if (muls->avgCount == 0) {
			wave->WriteDiffPat(avgName);
		}
		else {
			wave->ReadAvgArray(avgName);
			for (ix=0;ix<muls->nx*muls->ny;ix++) {
				wave->avgArray[0][ix] = (muls->avgCount*wave->avgArray[0][ix]+wave->diffpat[0][ix])/(muls->avgCount+1);
			}
			wave->WriteAvgArray(avgName);
		}
	}

	// Divide each image by its number of averages again:
	for (i=0;i<muls->detectorNum;i++) {
		// add intensity squared to image2 for this detector and pixel, then rescale:
		detectors[t][i]->image2[wave->detPosX][wave->detPosY] += detectors[t][i]->error*detectors[t][i]->error;
		detectors[t][i]->image2[wave->detPosX][wave->detPosY] /= detectors[t][i]->Navg+1;

		// do the rescaling for the average image:
		detectors[t][i]->image[wave->detPosX][wave->detPosY] /= detectors[t][i]->Navg+1;
	}
}

/*****  saveSTEMImages *******/
// Saves all detector images (STEM images) that are defined in muls.
//   When saving intermediate STEM images is enabled, this also saves
//   the intermediate STEM images for each detector.
void saveSTEMImages(MULS *muls)
{
	int i, ix, islice;
	double intensity;
	static char fileName[256];
	//imageStruct *header = NULL;
	std::vector<DetectorPtr> detectors;
	float t;

	int tCount = (int)(ceil((double)((muls->slices * muls->cellDiv) / muls->outputInterval)));

	// Loop over slices (intermediates)
	for (islice=0; islice <= tCount; islice++)
	{
		if (islice<tCount)
		{
			t = ((islice+1) * muls->outputInterval ) * muls->sliceThickness;
		}
		else
		{
			t = muls->slices*muls->cellDiv*muls->sliceThickness;
		}
		detectors = muls->detectors[islice];
		// write the output STEM images:
		// This is done only after all pixels have completed, so that image is complete.
		for (i=0; i<muls->detectorNum; i++)
		{
			// calculate the standard error for this image:
			detectors[i]->error = 0;
			intensity             = 0;
			for (ix=0; ix<muls->scanXN * muls->scanYN; ix++)
			{
				detectors[i]->error += (detectors[i]->image2[0][ix]-
					detectors[i]->image[0][ix] * detectors[i]->image[0][ix]);
				intensity += detectors[i]->image[0][ix] * detectors[i]->image[0][ix];
			}
			detectors[i]->error /= intensity;
			if (islice <tCount)
				sprintf(fileName,"%s/%s_%d.img", muls->folder, detectors[i]->name, islice);
			else
				sprintf(fileName,"%s/%s.img", muls->folder, detectors[i]->name, islice);
			//detectors[i]->SetComment(detectors[i]->name);
			// NOTE: the comment for STEM images must be this, or else the MATLAB GUI doesn't recognize it as a STEM image!
			//     That means the quantification and source size dialogs will be disabled.
			detectors[i]->SetComment("STEM image");
			detectors[i]->SetThickness(t);
			detectors[i]->SetParameter(0, (double)muls->avgCount+1);
			detectors[i]->SetParameter(1, (double)detectors[i]->error);

			for (ix=0; ix<muls->scanXN * muls->scanYN; ix++)
			{
				detectors[i]->SetParameter(2+ix, (double)detectors[i]->image2[0][ix]);
			}
			detectors[i]->WriteImage(fileName);
		}
	}
}

void readStartWave(WavePtr wave) {
	wave->ReadWave(wave->fileStart);
}


/******************************************************************
* propagate_slow()
* replicates the original way, mulslice did it:
*****************************************************************/
void propagate_slow(void **w,int nx, int ny,MULS *muls)
{
	int ixa, iya;
	real wr, wi, tr, ti,ax,by;
	real scale,t,dz;
	real dzs=0;
	real *propxr=NULL,*propyr=NULL;
	real *propxi=NULL,*propyi=NULL;
	real *kx2,*ky2;
	real *kx,*ky;
	real k2max=0,wavlen;
#if FLOAT_PRECISION == 1
	fftwf_complex **wave;
	wave = (fftwf_complex **)w;
#else
	fftw_complex **wave;
	wave = (fftw_complex **)w;
#endif

	ax = (*muls).resolutionX*nx;
	by = (*muls).resolutionY*ny;
	dz = (*muls).cz[0];

	if (dz != dzs) {
		if (propxr == NULL) {
			propxr = float1D(nx, "propxr" );
			propxi = float1D(nx, "propxi" );
			propyr = float1D(ny, "propyr" );
			propyi = float1D(ny, "propyi" );
			kx2    = float1D(nx, "kx2" );
			kx     = float1D(nx, "kx" );
			ky2    = float1D(ny, "ky2" );
			ky     = float1D(ny, "ky" );
		}
		dzs = dz;
		scale = dz*PI;
		wavlen = wavelength((*muls).v0);

		for( ixa=0; ixa<nx; ixa++) {
			kx[ixa] = (ixa>nx/2) ? (real)(ixa-nx)/ax :
				(real)ixa/ax;
			kx2[ixa] = kx[ixa]*kx[ixa];
			t = scale * (kx2[ixa]*wavlen);
			propxr[ixa] = (real)  cos(t);
			propxi[ixa] = (real) -sin(t);
		}
		for( iya=0; iya<ny; iya++) {
			ky[iya] = (iya>ny/2) ?
				(real)(iya-ny)/by :
			(real)iya/by;
			ky2[iya] = ky[iya]*ky[iya];
			t = scale * (ky2[iya]*wavlen);
			propyr[iya] = (real)  cos(t);
			propyi[iya] = (real) -sin(t);
		}
		k2max = nx/(2.0F*ax);
		if (ny/(2.0F*by) < k2max ) k2max = ny/(2.0F*by);
		k2max = 2.0/3.0 * k2max;
		// TODO: modifying shared value from multiple threads?
		k2max = k2max*k2max;
		(*muls).kx2=kx2;
		(*muls).ky2=ky2;
		(*muls).kx=kx;
		(*muls).ky=ky;
	}
	/* end of: if dz != dzs */
	/*************************************************************/

	/*************************************************************
	* Propagation
	************************************************************/
	for( ixa=0; ixa<nx; ixa++) {
		if( kx2[ixa] < k2max ) {
			for( iya=0; iya<ny; iya++) {
				if( (kx2[ixa] + ky2[iya]) < k2max ) {

					wr = wave[ixa][iya][0];
					wi = wave[ixa][iya][1];
					tr = wr*propyr[iya] - wi*propyi[iya];
					ti = wr*propyi[iya] + wi*propyr[iya];
					wave[ixa][iya][0] = tr*propxr[ixa] - ti*propxi[ixa];
					wave[ixa][iya][1] = tr*propxi[ixa] + ti*propxr[ixa];

				} else
					wave[ixa][iya][0] = wave[ixa][iya][1] = 0.0F;
			} /* end for(iy..) */

		} else for( iya=0; iya<ny; iya++)
			wave[ixa][iya][0] = wave[ixa][iya][1] = 0.0F;
	} /* end for(ix..) */
} /* end propagate_slow() */


/*------------------------ transmit() ------------------------*/
/*
transmit the wavefunction thru one layer
(simply multiply wave by transmission function)

waver,i[ix][iy]  = real and imaginary parts of wavefunction
transr,i[ix][iy] = real and imag parts of transmission functions

nx, ny = size of array

on entrance waver,i and transr,i are in real space

only waver,i will be changed by this routine
*/
void transmit(void **wave, void **trans,int nx, int ny,int posx,int posy) {
	int ix, iy;
	double wr, wi, tr, ti;

#if FLOAT_PRECISION == 1
	fftwf_complex **w, **t;
	w = (fftwf_complex **)wave;
	t = (fftwf_complex **)trans;
#else
	fftw_complex **w,**t;
	w = (fftw_complex **)wave;
	t = (fftw_complex **)trans;
#endif
	/*  trans += posx; */
	for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
		wr = w[ix][iy][0];
		wi = w[ix][iy][1];
		tr = t[ix+posx][iy+posy][0];
		ti = t[ix+posx][iy+posy][1];

		w[ix][iy][0] = wr*tr - wi*ti;
		w[ix][iy][1] = wr*ti + wi*tr;
	} /* end for(iy.. ix .) */
} /* end transmit() */

void fft_normalize(void **array,int nx, int ny) {
	int ix,iy;
	double fftScale;
#if FLOAT_PRECISION == 1
	fftwf_complex **carray;
	carray = (fftwf_complex **)array;
#else
	fftw_complex **carray;
	carray = (fftw_complex **)array;
#endif

	fftScale = 1.0/(double)(nx*ny);
	for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
		carray[ix][iy][0] *= fftScale;
		carray[ix][iy][1] *= fftScale;
	}
}

void showPotential(fftw_complex ***pot,int nz,int nx,int ny,double dx,double dy,double dz) {
	char *fileName = "potential.dat";
	char systStr[256];
	FILE *fpPot;
	int ix,iz;
	static fftw_complex *data = NULL;
	int length;
	float r;


	/*
	make data array:
	*/
	length = (nx < ny) ? nx : ny ;
	length +=2;
	if (data == NULL)
		data = (fftw_complex *)malloc(length * sizeof(fftw_complex));

	/*
	copy data to array
	*/

	for (ix=0;ix<nx;ix++) {
		data[ix][0] = pot[0][ix][ix][0];
		data[ix][1] = pot[0][ix][ix][1];
		for (iz=1;iz<nz;iz++) {
			data[ix][0] += 2*pot[iz][ix][ix][0];
			data[ix][1] += 2*pot[iz][ix][ix][1];
		}
		/*    printf("ix: %d, pot: %g\n",ix,data[ix]); */
	}

	if ( (fpPot = fopen( fileName, "w" )) == NULL ) {
		printf("Could not open %s for writing!\n",fileName);
		return;
	}
	for (ix=0;ix<nx;ix++) {
		r = sqrt(ix*ix*(dx*dx+dy*dy));
		fprintf( fpPot, "%g", r );
		fprintf( fpPot, "\t%g\t%g", data[ix][0], data[ix][1] );
		/*    for (iz = 0;iz < ((nz>10) ? 10 : nz);iz++) {
		r = sqrt(ix*ix*(dx*dx+dy*dy)+iz*iz*dz*dz);
		fprintf(fpPot,"\t%g",pot[iz][ix][ix][0]*r);
		}
		*/
		fprintf( fpPot, "\n" );
	}
	fclose( fpPot );

	sprintf(systStr,"xmgr -nxy %s &",fileName);
	system(systStr);
}


/*****************************************************************
* This function will write a data file with the pendeloesungPlot
* for selected beams
****************************************************************/
void writeBeams(MULS *muls, WavePtr wave, int ilayer, int absolute_slice) {
	static char fileAmpl[32];
	static char filePhase[32];
	static char fileBeam[32];
	static FILE *fp1 = NULL,*fpAmpl = NULL,*fpPhase=NULL;
	int ib;
	static int *hbeam=NULL,*kbeam=NULL;
	static real zsum = 0.0f,scale;
	real rPart,iPart,ampl,phase;
	static char systStr[64];
	// static int counter=0;

	if (!muls->lbeams)
		return;

	if ((muls->mode != REFINE) && ((*muls).mode != CBED)) {
		if (ilayer < 0) {
			if (fp1 != NULL) fclose(fp1);
			if (fpAmpl != NULL) fclose(fpAmpl);
			if (fpPhase != NULL) fclose(fpPhase);
			fp1 = fpAmpl = fpPhase = NULL;
			sprintf(systStr,"xmgr -nxy %s &",fileAmpl);
			system(systStr);
			return;
		}

		if ((fp1 == NULL) || (fpAmpl == NULL) || (fpPhase == NULL)) {
			scale = 1.0F / ( ((real)muls->nx) * ((real)muls->ny) );
			hbeam = (*muls).hbeam;
			kbeam = (*muls).kbeam;
			if ((hbeam == NULL) || (kbeam == NULL)) {
				printf("ERROR: hbeam or kbeam == NULL!\n");
				exit(0);
			}

			sprintf(fileAmpl,"%s/beams_amp.dat",(*muls).folder);
			sprintf(filePhase,"%s/beams_phase.dat",(*muls).folder);
			sprintf(fileBeam,"%s/beams_all.dat",(*muls).folder);
			fp1 = fopen(fileBeam, "w" );
			fpAmpl = fopen( fileAmpl, "w" );
			fpPhase = fopen( filePhase, "w" );
			if(fp1==NULL) {
				printf("can't open file %s\n", fileBeam);
				exit(0);
			}
			if(fpAmpl==NULL) {
				printf("can't open amplitude file %s\n",fileAmpl);
				exit(0);
			}
			if(fpPhase==NULL) {
				printf("can't open phase file %s\n", filePhase);
				exit(0);
			}
			fprintf(fp1, " (h,k) = ");
			for(ib=0; ib<(*muls).nbout; ib++) {
				fprintf(fp1," (%d,%d)", muls->hbeam[ib],  muls->kbeam[ib]);
			}
			fprintf( fp1, "\n" );
			fprintf( fp1, "nslice, (real,imag) (real,imag) ...\n\n");
			for( ib=0; ib<muls->nbout; ib++)
			{
				// printf("beam: %d [%d,%d]",ib,hbeam[ib],kbeam[ib]);
				if(hbeam[ib] < 0 ) hbeam[ib] = muls->nx + hbeam[ib];
				if(kbeam[ib] < 0 ) kbeam[ib] = muls->ny + kbeam[ib];
				if(hbeam[ib] < 0 ) hbeam[ib] = 0;
				if(kbeam[ib] < 0 ) kbeam[ib] = 0;
				if(hbeam[ib] > muls->nx-1 ) hbeam[ib] = muls->nx-1;
				if(kbeam[ib] > muls->ny-1 ) kbeam[ib] = muls->ny-1;
				// printf(" => [%d,%d] %d %d\n",hbeam[ib],kbeam[ib],muls->nx,muls->ny);
			}
			/****************************************************/
			/* setup of beam files, include the t=0 information */
			fprintf( fpAmpl, "%g",0.0);
			fprintf( fpPhase, "%g",0.0);
			for( ib=0; ib<muls->nbout; ib++) {
				ampl = 0.0;
				if ((hbeam[ib] == 0) && (kbeam[ib]==0))
					ampl = 1.0;
				fprintf(fpAmpl,"\t%g",ampl);
				fprintf(fpPhase,"\t%g",0.0);
			}
			fprintf( fpAmpl, "\n");
			fprintf( fpPhase, "\n");
		} /* end of if fp1 == NULL ... i.e. setup */


		zsum += (*muls).cz[ilayer];

		fprintf( fp1, "%g", zsum);
		fprintf( fpAmpl, "%g",zsum);
		fprintf( fpPhase, "%g",zsum);
		for( ib=0; ib<(*muls).nbout; ib++) {
			fprintf(fp1, "\t%g\t%g",
				rPart = scale*(*wave).wave[hbeam[ib]][kbeam[ib]][0],
				iPart = scale*(*wave).wave[hbeam[ib]][kbeam[ib]][1]);
			ampl = (real)sqrt(rPart*rPart+iPart*iPart);
			phase = (real)atan2(iPart,rPart);
			fprintf(fpAmpl,"\t%g",ampl);
			fprintf(fpPhase,"\t%g",phase);
		}
		fprintf( fp1, "\n");
		fprintf( fpAmpl, "\n");
		fprintf( fpPhase, "\n");
	} /* end of if muls.mode != REFINE */

	if (muls->mode == TEM) {
		if (muls->pendelloesung == NULL) {
			muls->pendelloesung =
				float2D((*muls).nbout,
				(*muls).slices*(*muls).mulsRepeat1*(*muls).mulsRepeat2*(*muls).cellDiv,
				"pendelloesung");
			scale = 1.0/(muls->nx*muls->ny);
			printf("Allocated memory for pendelloesung plot (%d x %d)\n",
				(*muls).nbout,(*muls).slices*(*muls).mulsRepeat1*(*muls).mulsRepeat2);
		}
		for( ib=0; ib<muls->nbout; ib++) {
			rPart = (*wave).wave[muls->hbeam[ib]][muls->kbeam[ib]][0];
			iPart = (*wave).wave[muls->hbeam[ib]][muls->kbeam[ib]][1];
			muls->pendelloesung[ib][absolute_slice] = scale*(real)(rPart*rPart+iPart*iPart);
			// printf("slice: %d beam: %d [%d,%d], intensity: %g\n",muls->nslic0,ib,muls->hbeam[ib],muls->kbeam[ib],muls->pendelloesung[ib][muls->nslic0]);
		} // end of ib=0 ...
	}

}







//////////////////////////////////////////////////////////////////////////////////
// DO NOT use this function, gives wrong results.
//
// This method of sampling z at the slice interval gives very wrong HOLZ results,
// since this means a linear interpolation of a highly non-linear function.  Atoms
// are quasi split up into two layers.  This causes a cutting in half of the z-periodicity
// for some unit cells (e.g. STO sampled 2 slices per unit cell)
//
#define S_SCALE 0.5
#define PHI_SCALE 47.87658
#define SHOW_SINGLE_POTENTIAL 1
fftwf_complex *getAtomPotential3D_3DFFT(int Znum, MULS *muls,double B) {
	int ix,iy,iz,iiz,ind3d,ind3dd,iKind;
	double zScale,kzmax,kzborder;
	fftwf_plan plan;
	static double f,phase,s2,kmax2,kx,ky,kz,dkx,dky,dkz,dx2,dy2,dz2;
	static int nx,ny,nz;
	static fftwf_complex **atPot = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
	//static imageStruct *header = NULL;

	fftwf_complex *ptr = NULL;
	char fileName[256];
	ImageIOPtr imageIO = ImageIOPtr(new CImageIO(muls->potNx,muls->potNy,
				muls->sliceThickness,muls->resolutionX/OVERSAMP_X,
				muls->resolutionY/OVERSAMP_X));;
#endif
	static double *splinb=NULL;
	static double *splinc=NULL;
	static double *splind=NULL;


	// scattering factors in:
	// float scatPar[4][30]
	if (atPot == NULL) {
		splinb = double1D(30, "splinb" );
		splinc = double1D(30, "splinc" );
		splind = double1D(30, "splind" );


		// Why do I use nx+2 and make dkx=1/nx? ... don't know anymore. :(
		nx = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX)+2;
		ny = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY)+2;
		// make the atom sphere have an even number of slices in z-direction
		// the 2 center slices will be equal.
		// The FFT-resolution in the z-direction must be high enough to avoid
		// artifacts due to premature cutoff of the rec. space scattering factor
		nz = 2*OVERSAMP_Z*(int)ceil(muls->atomRadius/muls->sliceThickness);
		dkx = OVERSAMP_X/((nx-2)*muls->resolutionX);  // 0.5-factor s->k
		dky = OVERSAMP_X/((ny-2)*muls->resolutionY);
		dkz = OVERSAMP_Z/(double)(nz*muls->sliceThickness);
		// if (nz<=2) dkz=0;
		// else dkz = 1/((nz-2)*muls->sliceThickness);
		kmax2 = 0.5*nx*dkx/(double)OVERSAMP_X;  // largest k that we'll admit

		printf("Cutoff scattering angle:kmax=%g, smax=%g (1/A)\n",kmax2,S_SCALE*kmax2);
		scatPar[0][29] = 1.2*S_SCALE*kmax2;
		scatPar[0][28] = 1.1*S_SCALE*kmax2;
		scatPar[0][27] = S_SCALE*kmax2;
		if (scatPar[0][26] > scatPar[0][27]) {
			// set additional scattering parameters to zero:
			for (ix = 0;ix < 20;ix++) {
				if (scatPar[0][26-ix] < scatPar[0][27]-0.1*(ix+1)) break;
				scatPar[0][26-ix] = scatPar[0][27]-0.1*(ix+1);
				for (iy=1; iy<=8;iy++) scatPar[iy][26-ix] = 0;
			}
			if (muls->printLevel > 1)
				printf("getAtomPotential3D: reduced angular range of scattering factor to %g/A!\n",scatPar[0][26-ix]);
		}
		kmax2 *= kmax2;

		atPot = (fftwf_complex **)malloc((NZMAX+1)*sizeof(fftwf_complex *));
		for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
	}
	// initialize this atom, if it has not been done yet:
	if (atPot[Znum] == NULL) {
		switch (Znum) {
	case 38: iKind = 1; break;  // Sr
	case 22: iKind = 2; break;  // Ti
	case  8: iKind = 3; break;  // O
	case 49: iKind = 4; break;  // In
	case 15: iKind = 5; break;  // P
	case  2: iKind = 6; break;  // He
	case 17: iKind = 7; break;  // Cl
	case 14: iKind = 8; break;  // Si
	case 20: iKind = 9; break;  // Ca
	case 56: iKind = 10; break;  // Ba
	case 26: iKind = 11; break;  // Fe
	default:
		printf("This atom kind (%d) is not supported yet - sorry!\n",Znum);
		exit(0);
		}


		// setup cubic spline interpolation:
		splinh(scatPar[0],scatPar[iKind],splinb,splinc,splind,30);

		atPot[Znum] = (fftwf_complex*)fftwf_malloc(nx*ny*nz*sizeof(fftwf_complex));
		memset(atPot[Znum],0,nx*ny*nz*sizeof(fftwf_complex));
		kzmax    = dkz*nz/2.0;
		kzborder = dkz*(nz/(2*OVERSAMP_Z) -1);
		for (iz=0;iz<nz;iz++) {
			kz = dkz*(iz<nz/2 ? iz : iz-nz);
			// We also need to taper off the potential in z-direction
			// in order to avoid cutoff artifacts.
			zScale = fabs(kz) <= kzborder ? 1.0 :
				0.5+0.5*cos(PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
			// printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
			for (ix=0;ix<nx;ix++) {
				kx = dkx*(ix<nx/2 ? ix : ix-nx);
				for (iy=0;iy<ny;iy++) {
					ky = dky*(iy<ny/2 ? iy : iy-ny);
					s2 = S_SCALE*S_SCALE*(kx*kx+ky*ky+kz*kz);
					// if this is within the allowed circle:
					if (s2<S_SCALE*S_SCALE*kmax2) {
						ind3d = iy+ix*ny+iz*nx*ny;
						// f = fe3D(Znum,k2,muls->tds,1.0,muls->scatFactor);
						// multiply scattering factor with Debye-Waller factor:
						// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
						f = seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,30,sqrt(s2))*exp(-s2*B);
						// note that the factor 2 is missing in the phase (2pi k*r)
						// this places the atoms in the center of the box.
						phase = PI*(kx*muls->resolutionX*nx/(OVERSAMP_X)+ky*muls->resolutionY*ny/(OVERSAMP_X));
						phase += PI*kz/OVERSAMP_Z*(muls->sliceThickness*(nz+1));
						atPot[Znum][ind3d][0] = zScale*f*cos(phase);
						atPot[Znum][ind3d][1] = zScale*f*sin(phase);
						// if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
					}
				}
			}
		}

		plan = fftwf_plan_dft_3d(nz,nx,ny,atPot[Znum],atPot[Znum],FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
		dx2 = muls->resolutionX*muls->resolutionX/(OVERSAMP_X*OVERSAMP_X);
		dy2 = muls->resolutionY*muls->resolutionY/(OVERSAMP_X*OVERSAMP_X);
		dz2 = muls->sliceThickness*muls->sliceThickness/(OVERSAMP_Z*OVERSAMP_Z);
		// Here we make sure that our atom is not bigger than the desired radius, i.e. that
		// we really have round blobs.
		// We also make sure that the potential touches zero at least somewhere.  This will avoid
		// sharp edges that could produce ringing artifacts.
		// It is certainly debatable whether this is a good apprach, or not.
		// printf("Setting up %d x %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,ny,nz,Znum,iKind,dkx*dky*dkz);
		// min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
		for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) for (iz=0;iz<nz/OVERSAMP_Z;iz++) {
			ind3d = iy+ix*ny+iz*nx*ny;
			// Integrate over 3 neighboring layers here:::::::::
			for (zScale=0,iiz=0;iiz<OVERSAMP_Z;iiz++) {
				ind3dd = iy+ix*ny+(OVERSAMP_Z*iz+iiz)*nx*ny;
				// if (sqrt((ix-nx/2)*(ix-nx/2)*dx2+(iy-ny/2)*(iy-ny/2)*dy2+(iz-0.5-nz/2)*(iz-0.5-nz/2)*dz2) > muls->atomRadius) {
				if ((sqrt((ix-nx/2)*(ix-nx/2)*dx2+(iy-ny/2)*(iy-ny/2)*dy2+
					(OVERSAMP_Z*iz+iiz+0.5-nz/2)*(OVERSAMP_Z*iz+iiz+0.5-nz/2)*dz2)) > muls->atomRadius) {
						// atPot[Znum][ind3d][0] *= dkx*dky*dkz;
						atPot[Znum][ind3dd][0] = 0;
				}
				zScale += atPot[Znum][ind3dd][0];
			}
			// assign the iz-th slice the sum of the 3 other slices:
			// and divide by unit cell area (volume, if in 3D):
			atPot[Znum][ind3d][0] = zScale*dkx*dky/nz;
			if (atPot[Znum][ind3d][0] < 0) atPot[Znum][ind3d][0] = 0;
			// if (atPot[Znum][ind3d][0] < min) min = atPot[Znum][ind3d][0];

			atPot[Znum][ind3d][1]= 0;
		}
		// printf("Found minimum potential value of %g ... subtracting it from 3D potential.\n",min);
		for (zScale = 0,ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++)  for (iz=0;iz<nz/OVERSAMP_Z;iz++){
			if ((sqrt((ix-nx/2)*(ix-nx/2)*dx2+(iy-ny/2)*(iy-ny/2)*dy2+
				(iz+0.5-nz/(2*OVERSAMP_Z))*(iz+0.5-nz/(2*OVERSAMP_Z))*dz2)) <= muls->atomRadius) {
					//      if (sqrt((ix-nx/2)*(ix-nx/2)*dx2+(iy-ny/2)*(iy-ny/2)*dy2+(iz-nz/2)*(iz-nz/2)*dz2) <= muls->atomRadius) {
					ind3d = iy+ix*ny+iz*nx*ny;
					if ((ix == nx/2) && (iy == ny/2)) zScale += atPot[Znum][ind3d][0];
					// atPot[Znum][ind3d][0] -= min;
			}
		}
		// make sure we don't produce negative potential:
		// if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny][0] -= min;
#if SHOW_SINGLE_POTENTIAL == 1
		for (iz=0;iz<nz/OVERSAMP_Z;iz++) {
			sprintf(fileName,"potential_%d_%d.img",Znum,iz);
			imageIO->SetThickness(iz);
			ptr = &(atPot[Znum][iz*nx*ny]);
			imageIO->WriteRealImage((void **)ptr,fileName);
		}
#endif
		if (muls->printLevel > 0)
			printf("Created 3D %d x %d x %d potential array for Z=%d (%d, B=%g, sum=%g)\n",nx,ny,nz/OVERSAMP_Z,Znum,iKind,B,zScale);
	}
	return atPot[Znum];
}
#undef SHOW_SINGLE_POTENTIAL
