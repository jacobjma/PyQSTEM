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

#define VERSION 2.30
#define VIB_IMAGE_TEST

#ifndef _WIN32
#define UNIX
#endif
/* #define USE_FFT_POT */
// for memory leak checking in windows.  Should not affect speed of release builds.
#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef _WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif
#include <string.h>
#ifndef _WIN32
#ifdef __cplusplus
#include <cmath>
#else
#include <math.h>
#endif
#else
#include <math.h>
#endif

#include <time.h>
#include <ctype.h>
#include <sys/stat.h>
// #include <stat.h>

#include <omp.h>
#include <iostream>

#include "memory_fftw3.h"	/* memory allocation routines */
#include "readparams.h"
#include "imagelib_fftw3.h"
#include "fileio_fftw3.h"
#include "matrixlib.h"
#include "stemlib.h"
#include "stemutil.h"
// #include "weblib.h"
#include "customslice.h"
#include "data_containers.h"

#define NCINMAX 1024
#define NPARAM	64    /* number of parameters */
#define MAX_SCANS 1   /* maximum number of linescans per graph window */
#define PHASE_GRATING 0
#define BUF_LEN 256

#define DELTA_T 1     /* number of unit cells between pictures */
#define PICTS 5      /* number of different thicknesses */
#define NBITS 8	       /* number of bits for writeIntPix */
#define RAD2DEG 57.2958
#define SQRT_2 1.4142135

const char *resultPage = "result.html";
/* global variable: */
MULS muls;
// int fftMeasureFlag = FFTW_MEASURE;
int fftMeasureFlag = FFTW_ESTIMATE;
extern char *elTable;

void makeAnotation(real **pict,int nx,int ny,char *text);
void initMuls();
void writeIntPix(char *outFile,real **pict,int nx,int ny);
void runMuls(int lstart);
void saveLineScan(int run);
void readBeams(FILE *fpBeams);
void doCBED();
void doNBED();
void doSTEM();
void doTEM();
void doMSCBED();
void doTOMO();
void readFile();
void displayParams();

void usage() {
	printf("usage: stem [input file='stem.dat']\n\n");
}


/***************************************************************
***************************************************************
* MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN *
***************************************************************
**************************************************************/


int main(int argc, char *argv[]) {
	int i;
	double timerTot;
	char fileName[512];
	char cinTemp[BUF_LEN];

	timerTot = cputim();
	for (i=0;i<BUF_LEN;i++)
		cinTemp[i] = 0;
	muls.nCellX = 1; muls.nCellY = 1; muls.nCellZ = 1;

#ifdef UNIX
	system("date");
#endif

	/*************************************************************
	* read in the parameters
	************************************************************/
	if (argc < 2)
		sprintf(fileName,"stem.dat");
	else
		strcpy(fileName,argv[1]);
	if (parOpen(fileName) == 0)
	{
		printf("could not open input file %s!\n",fileName);
		usage();
		exit(0);
	}
	readFile();

	displayParams();
#ifdef _OPENMP
	omp_set_dynamic(1);
#endif
	if (muls.mode == STEM) {
		// sprintf(systStr,"mkdir %s",muls.folder);
		// system(systStr);
		for (muls.avgCount=0;muls.avgCount<muls.avgRuns;muls.avgCount++) {
			muls.dE_E = muls.dE_EArray[muls.avgCount];
			// printf("dE/E: %g\n",muls.dE_E);
		}

		// drawStructure();

		muls.showProbe = 0;
		muls.avgCount = -1;
	}

	switch (muls.mode) {
	  case CBED:   doCBED();   break;
	  case STEM:   doSTEM();   break;
	  case TEM:    doTEM();    break;
	  case MSCBED: doMSCBED(); break;
	  case TOMO:   doTOMO();   break;
	  case NBED:   doNBED();   break;
	  // case REFINE: doREFINE(); break;
	  default:
		  printf("Mode not supported\n");
	}

	parClose();
#if _DEBUG
	_CrtDumpMemoryLeaks();
#endif

	// Console.Read(); //  <--- Right here
	printf( "DEBUG Main: Press any key to exit()" );
	std::cin.ignore();
	return 0;
}

/***************************************************************
***************************************************************
** End of MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN  **
***************************************************************
**************************************************************/

void initMuls() {
	int sCount,i,slices;

	slices = muls.slices;

	/* general setup: */
	muls.lpartl = 0;

	muls.atomRadius = 5.0;  /* radius in A for making the potential boxes */

	for (sCount =0;sCount<slices;sCount++)
		muls.cin2[sCount] = 'a'+sCount;
	for (sCount = slices;sCount < NCINMAX;sCount++)
		muls.cin2[sCount] = 0;
	muls.nlayer = slices;
	muls.saveFlag = 0;

	muls.sigmaf = 0;
	muls.dfdelt = 0;
	muls.acmax = 0;
	muls.acmin = 0;
	muls.aobj = 0;
	muls.Cs = 0;
	muls.aAIS = 0;
	// muls.areaAIS = 1.0;

	// Tomography parameters:
	muls.tomoTilt = 0;
	muls.tomoStart = 0;
	muls.tomoStep = 0;
	muls.tomoCount = 0;  // indicate: NO Tomography simulation.

	/* make multislice read the inout files and assign transr and transi: */
	muls.trans = NULL;
	muls.cz = NULL;  // (real *)malloc(muls.slices*sizeof(real));

	muls.onlyFresnel = 0;
	muls.showPhaseplate = 0;
	muls.czOffset = 0;  /* defines the offset for the first slice in
						fractional coordinates        */
	muls.normHolog = 0;
	muls.gaussianProp = 0;


	muls.sparam = (float *)malloc(NPARAM*sizeof(float));
	for (i=0;i<NPARAM;i++)
		muls.sparam[i] = 0.0;
	muls.kx = NULL;
	muls.kx2= NULL;
	muls.ky = NULL;
	muls.ky2= NULL;

	/****************************************************/
	/* copied from slicecell.c                          */
	muls.pendelloesung = NULL;
}

int DirExists(char *filename) {
  struct stat status;
  status.st_mode = 0;

  // Attempt to get the file attributes
  if (stat(filename,&status) == 0) return 0;
  if (status.st_mode & S_IFDIR )  return 1;
  return 0;
}


/************************************************************************
*
***********************************************************************/
void displayProgress(int flag) {
	// static double timer;
	static double timeAvg = 0;
	static double intensityAvg = 0;
	static time_t time0,time1;
	double curTime;
	int jz;

	if (flag < 0) {
		time(&time0);
		// timer = cputim();
		return;
	}
	time(&time1);
	curTime = difftime(time1,time0);
	/*   curTime = cputim()-timer;
	if (curTime < 0) {
	printf("timer: %g, curr. time: %g, diff: %g\n",timer,cputim(),curTime);
	}
	*/
	if (muls.printLevel > 0) {

		if (muls.tds) {
			timeAvg = ((muls.avgCount)*timeAvg+curTime)/(muls.avgCount+1);
			intensityAvg = ((muls.avgCount)*intensityAvg+muls.intIntensity)/(muls.avgCount+1);
			printf("\n********************** run %3d ************************\n",muls.avgCount+1);
			// if (muls.avgCount < 1) {
			printf("* <u>: %3d |",muls.Znums[0]);
			for (jz=1;jz<muls.atomKinds;jz++) printf(" %8d |",muls.Znums[jz]);
			printf(" intensity | time(sec) |    chi^2  |\n");
			// }
			/*
			printf("* %9g | %9g | %9g \n",muls.u2,muls.intIntensity,curTime);
			}
			else {
			*/
			printf("*");
			for (jz=0;jz<muls.atomKinds;jz++) printf(" %8f |",(float)(muls.u2[jz]));
			printf(" %9f | %9f | %9f |\n",muls.intIntensity,curTime,muls.avgCount > 0 ? muls.chisq[muls.avgCount-1] : 0);
			printf("*");
			for (jz=0;jz<muls.atomKinds;jz++) printf(" %8f |",(float)(muls.u2avg[jz]));
			printf(" %9f | %9f \n",intensityAvg,timeAvg);
		}
		else {
			printf("\n**************** finished after %.1f sec ******************\n",curTime);
		}
	}  // end of printLevel check.

	time(&time0);
	//  timer = cputim();

}

void displayParams() {
	FILE *fpDir;
	char systStr[64];
	double k2max,temp;
	int i,j;
	static char Date[16],Time[16];
	time_t caltime;
	struct tm *mytime;
	const double pi=3.1415926535897;

	if (muls.printLevel < 1) {
		if ((fpDir = fopen(muls.folder,"r"))) {
			fclose(fpDir);
			// printf(" (already exists)\n");
		}
		else {
			sprintf(systStr,"mkdir %s",muls.folder);
			system(systStr);
			// printf(" (created)\n");
		}
		return;
	}
	caltime = time( NULL );
	mytime = localtime( &caltime );
	strftime( Date, 12, "%Y:%m:%d", mytime );
	strftime( Time, 9, "%H:%M:%S", mytime );

	printf("\n*****************************************************\n");
	printf("* Running program STEM3 (version %.2f) in %s mode\n",VERSION,
		(muls.mode == STEM) ? "STEM" : (muls.mode==TEM) ? "TEM" :
		(muls.mode == CBED) ? "CBED" : (muls.mode==TOMO)? "TOMO" :
		(muls.mode == NBED) ? "NBED" :
		"???");
	printf("* Date: %s, Time: %s\n",Date,Time);
	printf("*****************************************************\n");
	printf("* Print level:          %d\n",muls.printLevel);
	printf("* Save level:           %d\n",muls.saveLevel);
	printf("* Input file:           %s\n",muls.atomPosFile);
	if (muls.savePotential)
		printf("* Potential file name:  %s\n",muls.fileBase);
	/* create the data folder ... */
	printf("* Data folder:          ./%s/ ",muls.folder);
	if (DirExists(muls.folder)) {
	// if ((fpDir = fopen(muls.folder,"r"))!= NULL) {
	// 	fclose(fpDir);
		printf(" (already exists)\n");
	}
	else {
		sprintf(systStr,"mkdir %s",muls.folder);
		system(systStr);
		printf(" (created)\n");
	}

	if ((muls.cubex == 0) || (muls.cubey == 0) || (muls.cubez == 0))
		printf("* Unit cell:            ax=%g by=%g cz=%g\n",
		muls.ax,muls.by,muls.c);
	else {
		printf("* Size of Cube:         ax=%g by=%g cz=%g\n",
			muls.cubex,muls.cubey,muls.cubez);
		printf("* Cube size adjusted:   %s\n",muls.adjustCubeSize ? "yes" : "no");
	}

	printf("* Super cell:           %d x %d x %d unit cells\n",muls.nCellX,muls.nCellY,muls.nCellZ);
	printf("* Number of atoms:      %d (super cell)\n",muls.natom);
	printf("* Crystal tilt:         x=%g deg, y=%g deg, z=%g deg\n",
		muls.ctiltx*RAD2DEG,muls.ctilty*RAD2DEG,muls.ctiltz*RAD2DEG);
	printf("* Beam tilt:            x=%g deg, y=%g deg (tilt back == %s)\n",muls.btiltx*RAD2DEG,muls.btilty*RAD2DEG,
		(muls.tiltBack == 1 ? "on" : "off"));
	printf("* Model dimensions:     ax=%gA, by=%gA, cz=%gA (after tilt)\n"
		"*                       sampled every %g x %g x %g A\n",
		muls.ax,muls.by,muls.c,muls.resolutionX,muls.resolutionY,muls.sliceThickness);
	printf("* Atom species:         %d (Z=%d",muls.atomKinds,muls.Znums[0]);
	for (i=1;i<muls.atomKinds;i++) printf(", %d",muls.Znums[i]); printf(")\n");
	printf("* Super cell divisions: %d (in z direction) %s\n",muls.cellDiv,muls.equalDivs ? "equal" : "non-equal");
	printf("* Slices per division:  %d (%gA thick slices [%scentered])\n",
		muls.slices,muls.sliceThickness,(muls.centerSlices) ? "" : "not ");
	printf("* Output every:         %d slices\n",muls.outputInterval);

	printf("* Potential:            ");
	if (muls.potential3D) printf("3D"); else printf("2D");
	if (muls.fftpotential) printf(" (fast method)\n"); else printf(" (slow method)\n");
	printf("* Pot. array offset:    (%g,%g,%g)A\n",muls.potOffsetX,muls.potOffsetY,muls.czOffset);
	printf("* Potential periodic:   (x,y): %s, z: %s\n",
		(muls.nonPeriod) ? "no" : "yes",(muls.nonPeriodZ) ? "no" : "yes");

	printf("* Beams:                %d x %d \n",muls.nx,muls.ny);
	printf("* Acc. voltage:         %g (lambda=%gA)\n",muls.v0,wavelength(muls.v0));
	printf("* C_3 (C_s):            %g mm\n",muls.Cs*1e-7);
	printf("* C_1 (Defocus):        %g nm%s\n",0.1*muls.df0,
		(muls.Scherzer == 1) ? " (Scherzer)" : (muls.Scherzer==2) ? " (opt.)":"");
	printf("* Astigmatism:          %g nm, %g deg\n",0.1*muls.astigMag,RAD2DEG*muls.astigAngle);

	// more aberrations:
	if (muls.a33 > 0)
		printf("* a_3,3:                %g nm, phi=%g deg\n",muls.a33*1e-1,muls.phi33*RAD2DEG);
	if (muls.a31 > 0)
		printf("* a_3,1:                %g nm, phi=%g deg\n",muls.a31*1e-1,muls.phi31*RAD2DEG);

	if (muls.a44 > 0)
		printf("* a_4,4:                %g um, phi=%g deg\n",muls.a44*1e-4,muls.phi44*RAD2DEG);
	if (muls.a42 > 0)
		printf("* a_4,2:                %g um, phi=%g deg\n",muls.a42*1e-4,muls.phi42*RAD2DEG);

	if (muls.a55 > 0)
		printf("* a_5,5:                %g um, phi=%g deg\n",muls.a55*1e-4,muls.phi55*RAD2DEG);
	if (muls.a53 > 0)
		printf("* a_5,3:                %g um, phi=%g deg\n",muls.a53*1e-4,muls.phi53*RAD2DEG);
	if (muls.a51 > 0)
		printf("* a_5,1:                %g um, phi=%g deg\n",muls.a51*1e-4,muls.phi51*RAD2DEG);

	if (muls.a66 > 0)
		printf("* a_6,6:                %g um, phi=%g deg\n",muls.a66*1e-7,muls.phi66*RAD2DEG);
	if (muls.a64 > 0)
		printf("* a_6,4:                %g um, phi=%g deg\n",muls.a64*1e-7,muls.phi64*RAD2DEG);
	if (muls.a62 > 0)
		printf("* a_6,2:                %g um, phi=%g deg\n",muls.a62*1e-7,muls.phi62*RAD2DEG);
	if (muls.C5 != 0)
		printf("* C_5:                  %g mm\n",muls.C5*1e-7);

	printf("* C_c:                  %g mm\n",muls.Cc*1e-7);
	printf("* Aperture half angle:  %g mrad\n",muls.alpha);
	printf("* AIS aperture:         ");
	if (muls.aAIS > 0) printf("%g A\n",muls.aAIS);
	else printf("none\n");
	printf("* beam current:         %g pA\n",muls.beamCurrent);
	printf("* dwell time:           %g msec (%g electrons)\n",
		muls.dwellTime,muls.electronScale);

	printf("* Damping dE/E: %g / %g \n",sqrt(muls.dE_E*muls.dE_E+muls.dV_V*muls.dV_V+muls.dI_I*muls.dI_I)*muls.v0*1e3,muls.v0*1e3);

	/*
	if (muls.ismoth) printf("Type 1 (=smooth aperture), ");
	if (muls.gaussFlag) printf("will apply gaussian smoothing");
	printf("\n");
	*/
	printf("* Temperature:          %gK\n",muls.tds_temp);
	if (muls.tds)
		printf("* TDS:                  yes (%d runs)\n",muls.avgRuns);
	else
		printf("* TDS:                  no\n");
	if (muls.imageGamma == 0)
		printf("* Gamma for diff. patt: logarithmic\n");
	else
		printf("* Gamma for diff. patt: %g\n",muls.imageGamma);

	/**********************************************************
	* FFTW specific data structures (stores in row major order)
	*/
	if (fftMeasureFlag == FFTW_MEASURE)
		printf("* Probe array:          %d x %d pixels (optimized)\n",muls.nx,muls.ny);
	else
		printf("* Probe array:          %d x %d pixels (estimated)\n",muls.nx,muls.ny);
	printf("*                       %g x %gA\n",
		muls.nx*muls.resolutionX,muls.ny*muls.resolutionY);

	if (fftMeasureFlag == FFTW_MEASURE)
		printf("* Potential array:      %d x %d (optimized)\n",muls.potNx,muls.potNy);
	else
		printf("* Potential array:      %d x %d (estimated)\n",muls.potNx,muls.potNy);
	printf("*                       %g x %gA\n",muls.potSizeX,muls.potSizeY);
	printf("* Scattering factors:   %d\n",muls.scatFactor);
	/***************************************************/
	/*  printf("Optimizing fftw plans according to probe array (%d x %dpixels = %g x %gA) ...\n",
	muls.nx,muls.ny,muls.nx*muls.resolutionX,muls.ny*muls.resolutionY);
	*/
	k2max = muls.nx/(2.0*muls.potSizeX);
	temp = muls.ny/(2.0*muls.potSizeY);
	if( temp < k2max ) k2max = temp;
	k2max = (BW * k2max);

	printf("* Real space res.:      %gA (=%gmrad)\n",
		1.0/k2max,wavelength(muls.v0)*k2max*1000.0);
	printf("* Reciprocal space res: dkx=%g, dky=%g\n",
		1.0/(muls.nx*muls.resolutionX),1.0/(muls.ny*muls.resolutionY));
	if (muls.mode == STEM) {
		printf("*\n"
			"* STEM parameters:\n");
		printf("* Maximum scattering angle:  %.0f mrad\n",
			0.5*2.0/3.0*wavelength(muls.v0)/muls.resolutionX*1000);
		printf("* Number of detectors:  %d\n",muls.detectorNum);
		for (i=0;i<muls.detectorNum;i++) {
			printf("* %d (\"%s\"):",i+1,muls.detectors[0][i]->name);
			for (j=0;j<14-strlen(muls.detectors[0][i]->name);j++) printf(" ");
			printf(" %g .. %g mrad = (%.2g .. %.2g 1/A)\n",
				muls.detectors[0][i]->rInside,
				muls.detectors[0][i]->rOutside,
				muls.detectors[0][i]->k2Inside,
				muls.detectors[0][i]->k2Outside);
			if ((muls.detectors[0][i]->shiftX != 0) ||(muls.detectors[0][i]->shiftY != 0))
				printf("*   center shifted:     dkx=%g, dky=%g\n",
				muls.detectors[0][i]->shiftX,muls.detectors[0][i]->shiftY);
		}
		printf("* Scan window:          (%g,%g) to (%g,%g)A, %d x %d = %d pixels\n",
			muls.scanXStart,muls.scanYStart,muls.scanXStop,muls.scanYStop,
			muls.scanXN,muls.scanYN,muls.scanXN*muls.scanYN);
	} /* end of if mode == STEM */

	/***********************************************************************
	* TOMOGRAPHY Mode
	**********************************************************************/
	if (muls.mode == TOMO) {
		printf("*\n"
			"* TOMO parameters:\n");
		printf("* Starting angle:       %g mrad (%g deg)\n",
			muls.tomoStart,muls.tomoStart*0.18/pi);
		printf("* Angular increase:     %g mrad (%g deg)\n",
			muls.tomoStep,muls.tomoStep*0.180/pi);
		printf("* Number of dp's:       %d\n",muls.tomoCount);
		printf("* Zoom factor:          %g\n",muls.zoomFactor);


	}

	printf("*\n*****************************************************\n");

	/* k2max = muls.alpha * 0.001;  k2max = aperture in rad */
}


void readArray(char *title,double *array,int N) {
	int i;
	char buf[512],*str;

	if (!readparam(title,buf,1)) printf("%s array not found - exit\n",title), exit(0);
	i=0;
	str = buf;
	if (strchr(" \t\n",*str) != NULL) str = strnext(str," \t");
	while (i<N) {
		array[i++] = atof(str);
		str = strnext(str," \t\n");
		while (str == NULL) {
			if (!readNextLine(buf,511))
				printf("Incomplete reading of %s array - exit\n",title), exit(0);
			str = buf;
			if (strchr(" \t\n",*str) != NULL) str = strnext(str," \t\n");
		}
	}
}

/***********************************************************************
* readSFactLUT() reads the scattering factor lookup table from the
* input file
**********************************************************************/
void readSFactLUT() {
	int Nk,i,j;
	double **sfTable=NULL;
	double *kArray = NULL;
	char buf[256], elem[8];

	if (readparam("Nk:",buf,1))
		Nk = atoi(buf);
	else {
		printf("Could not find number of k-points for custom scattering factors (Nk)\n");
		exit(0);
	}

	// allocate memory for sfTable and kArray:
	sfTable = double2D(muls.atomKinds,Nk+1,"sfTable");
	kArray  = double1D(Nk+1,"kArray");

	// read the k-values:
	readArray("k:",kArray,Nk);
	kArray[Nk] = 2.0*kArray[Nk-1];

	for (j=0;j<muls.atomKinds;j++) {
		elem[3] = '\0';
		elem[2] = ':';
		elem[0] = elTable[2*muls.Znums[j]-2];
		elem[1] = elTable[2*muls.Znums[j]-1];
		if (elem[1] == ' ') {
			elem[1] = ':';
			elem[2] = '\0';
		}
		// printf("%s\n",elem);
		readArray(elem,sfTable[j],Nk);
		sfTable[j][Nk] = 0.0;
	}

	if (0) {
		printf("k: ");
		for (i=0;i<=Nk;i++) printf("%.3f ",kArray[i]);
		for (j=0;j<muls.atomKinds;j++) {
			printf("\n%2d: ",muls.Znums[j]);
			for (i=0;i<=Nk;i++) printf("%.3f ",sfTable[j][i]);
		}
		printf("\n");
	}
	muls.sfTable = sfTable;
	muls.sfkArray = kArray;
	muls.sfNk = Nk+1;
}


/************************************************************************
* readFile()
*
* reads the parameters from the input file and does some
* further setup accordingly
*
***********************************************************************/
void readFile() {
	char answer[256];
	FILE *fpTemp;
	float ax,by,c;
	char buf[BUF_LEN],*strPtr;
	int i,ix;
	int potDimensions[2];
	long ltime;
	unsigned long iseed;
	double dE_E0,x,y,dx,dy;
	const double pi=3.1415926535897;


	ltime = (long) time(NULL);
	iseed = (unsigned) ltime;
	muls.cubex = 0.0;
	muls.cubey = 0.0;
	muls.cubez = 0.0;

	muls.mode = STEM;
	if (readparam("mode:",buf,1)) {
		if (strstr(buf,"STEM")) muls.mode = STEM;
		else if (strstr(buf,"TEM")) muls.mode = TEM;
		else if (strstr(buf,"CBED")) muls.mode = CBED;
		else if (strstr(buf, "NBED")) muls.mode = NBED;
		else if (strstr(buf, "TOMO")) muls.mode = TOMO;
		else if (strstr(buf, "REFINE")) muls.mode = REFINE;
	}

	muls.printLevel = 2;
	if (readparam("print level:",buf,1)) sscanf(buf,"%d",&(muls.printLevel));
	muls.saveLevel = 0;
	if (readparam("save level:",buf,1)) sscanf(buf,"%d",&(muls.saveLevel));


	/************************************************************************
	* Basic microscope/simulation parameters,
	*/
	//muls.fileBase = (char *)malloc(512);
	//muls.atomPosFile = (char *)malloc(512);
	if ( !readparam( "filename:", buf, 1 ) )
	{
		perror( "Error:" );
		printf( "readFile did not find crystal .cfg file for parsing." );
		exit( 0 );
	}
	sscanf(buf,"%s",muls.fileBase);
	// printf("buf: %s\n",buf);
	// printf("fileBase: %s\n",muls.fileBase);

	// search for second '"', in case the filename is in quotation marks:
	if (muls.fileBase[0] == '"') {
		strPtr = strchr(buf,'"');
		strcpy(muls.fileBase,strPtr+1);
		strPtr = strchr(muls.fileBase,'"');
		*strPtr = '\0';
	}

	// RAM: add support to read-in a single wavefunction
	if ( readparam( "wavename:", buf, 1) )
	{
		sscanf(buf, "%s", muls.fileWaveIn);
	}
	if ( muls.fileWaveIn[0] == '"' )
	{ // search for second '"', in case the filename is in quotation marks:
		strPtr = strchr( buf, '"' );
		strcpy( muls.fileWaveIn, strPtr + 1 );
		strPtr = strchr( muls.fileWaveIn, '"' );
		*strPtr = '\0';
	}

	if (readparam("NCELLX:",buf,1)) sscanf(buf,"%d",&(muls.nCellX));
	if (readparam("NCELLY:",buf,1)) sscanf(buf,"%d",&(muls.nCellY));

	muls.cellDiv = 1;
	if (readparam("NCELLZ:",buf,1)) {
		sscanf(buf,"%s",answer);
		if ((strPtr = strchr(answer,'/')) != NULL) {
			strPtr[0] = '\0';
			muls.cellDiv = atoi(strPtr+1);
		}
		muls.nCellZ = atoi(answer);
	}

	/*************************************************
	* Read the beam tilt parameters
	*/
	muls.btiltx = 0.0;
	muls.btilty = 0.0;
	muls.tiltBack = 1;
	answer[0] = '\0';
	if (readparam("Beam tilt X:",buf,1)) {
		sscanf(buf,"%g %s",&(muls.btiltx),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.btiltx *= pi/180.0;
	}
	answer[0] = '\0';
	if (readparam("Beam tilt Y:",buf,1)) {
		sscanf(buf,"%g %s",&(muls.btilty),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.btilty *= pi/180.0;
	}
	if (readparam("Tilt back:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.tiltBack  = (tolower(answer[0]) == (int)'y');
	}


	/*************************************************
	* Read the crystal tilt parameters
	*/
	muls.ctiltx = 0.0;  /* tilt around X-axis in mrad */
	muls.ctilty = 0.0;  /* tilt around y-axis in mrad */
	muls.ctiltz = 0.0;  /* tilt around z-axis in mrad */
	answer[0] = '\0';
	if (readparam("Crystal tilt X:",buf,1)) {
		sscanf(buf,"%g %s",&(muls.ctiltx),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.ctiltx *= pi/180.0;
	}
	answer[0] = '\0';
	if (readparam("Crystal tilt Y:",buf,1)) {
		sscanf(buf,"%g %s",&(muls.ctilty),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.ctilty *= pi/180.0;
	}
	answer[0] = '\0';
	if (readparam("Crystal tilt Z:",buf,1)) {
		sscanf(buf,"%g %s",&(muls.ctiltz),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.ctiltz *= pi/180.0;
	}
	muls.cubex = 0; muls.cubey = 0; muls.cubez = 0;
	if (readparam("Cube:",buf,1)) {
		sscanf(buf,"%g %g %g",&(muls.cubex),&(muls.cubey),&(muls.cubez)); /* in A */
	}

	muls.adjustCubeSize = 0;
	if (readparam("Adjust cube size with tilt:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.adjustCubeSize  = (tolower(answer[0]) == (int)'y');
	}

	/***************************************************************************
	* temperature related data must be read before reading the atomic positions:
	***************************************************************************/
	if (readparam("tds:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.tds = (tolower(answer[0]) == (int)'y');
	}
	else muls.tds = 0;
	if (readparam("temperature:",buf,1)) sscanf(buf,"%g",&(muls.tds_temp));
	else muls.tds_temp = 300.0;
	muls.Einstein = 1;
	//muls.phononFile = NULL;
	if (readparam("phonon-File:",buf,1)) {
		sscanf(buf,"%s",muls.phononFile);
		muls.Einstein = 0;
	}

	/**********************************************************************
	* Read the atomic model positions !!!
	*********************************************************************/
	sprintf(muls.atomPosFile,muls.fileBase);
	/* remove directory in front of file base: */
	while ((strPtr = strchr(muls.fileBase,'\\')) != NULL) strcpy(muls.fileBase,strPtr+1);

	/* add a '_' to fileBase, if not existent */
	// RAM: ERROR THIS CODE SEEMS TO BREAK IF ONE USES UNDERSCORES IN FILENAME
	// CHANGE TO CHECK FOR UNDERSCORE AT END OF STRING
	if (strrchr(muls.fileBase,'_') != muls.fileBase+strlen(muls.fileBase)-1) {
		if ((strPtr = strchr(muls.fileBase,'.')) != NULL) sprintf(strPtr,"_");
		else strcat(muls.fileBase,"_");
	}

	if (strchr(muls.atomPosFile,'.') == NULL) {
		/*
		strPtr = strrchr(muls.atomPosFile,'_');
		if (strPtr != NULL)
		*(strPtr) = 0;
		*/
		// take atomPosFile as is, or add an ending to it, if it has none yet
		// RAM: This code can break if there's more than one period in a filename, TO DO: do strcmp's on last four/five characters instead
		if (strrchr(muls.atomPosFile,'.') == NULL)
		{
			sprintf(buf,"%s.cssr",muls.atomPosFile);
			if ((fpTemp=fopen(buf,"r")) == NULL)
			{
				sprintf(buf,"%s.cfg",muls.atomPosFile);
				if ((fpTemp=fopen(buf,"r")) == NULL)
				{
					printf("Could not find input file %s.cssr or %s.cfg\n",
						muls.atomPosFile,muls.atomPosFile);
					exit(0);
				}
				strcat(muls.atomPosFile,".cfg");
				fclose(fpTemp);
			}
			else
			{
				strcat(muls.atomPosFile,".cssr");
				fclose(fpTemp);
			}
		}
	}
	// We need to initialize a few variables, before reading the atomic
	// positions for the first time.
	muls.natom = 0;
	muls.atoms = NULL;
	muls.Znums = NULL;
	muls.atomKinds = 0;
	muls.u2 = NULL;
	muls.u2avg = NULL;

	muls.xOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam("xOffset:",buf,1)) sscanf(buf,"%g",&(muls.xOffset));
	muls.yOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam("yOffset:",buf,1)) sscanf(buf,"%g",&(muls.yOffset));
	// printf("Reading Offset: %f, %f\n",muls.xOffset,muls.yOffset);

	// the last parameter is handleVacancies.  If it is set to 1 vacancies
	// and multiple occupancies will be handled.
	// _CrtSetDbgFlag  _CRTDBG_CHECK_ALWAYS_DF();
	// printf("memory check: %d, ptr= %d\n",_CrtCheckMemory(),(int)malloc(32*sizeof(char)));

	muls.atoms = readUnitCell(&(muls.natom),muls.atomPosFile,&muls,1);

	// printf("memory check: %d, ptr= %d\n",_CrtCheckMemory(),(int)malloc(32*sizeof(char)));


	if (muls.atoms == NULL) {
		printf("Error reading atomic positions!\n");
		exit(0);
	}
	if (muls.natom == 0) {
		printf("No atom within simulation boundaries!\n");
		exit(0);
	}
	// printf("hello!\n");
	ax = muls.ax/muls.nCellX;
	by = muls.by/muls.nCellY;;
	c =  muls.c/muls.nCellZ;


	/*****************************************************************
	* Done reading atomic positions
	****************************************************************/

	if (!readparam("nx:",buf,1)) exit(0); sscanf(buf,"%d",&(muls.nx));
	if (readparam("ny:",buf,1)) sscanf(buf,"%d",&(muls.ny));
	else muls.ny = muls.nx;

	muls.resolutionX = 0.0;
	muls.resolutionY = 0.0;
	if (readparam("resolutionX:",buf,1)) sscanf(buf,"%g",&(muls.resolutionX));
	if (readparam("resolutionY:",buf,1)) sscanf(buf,"%g",&(muls.resolutionY));
	if (!readparam("v0:",buf,1)) exit(0); sscanf(buf,"%g",&(muls.v0));

	muls.centerSlices = 0;
	if (readparam("center slices:",buf,1)) {
		// answer[0] =0;
		sscanf(buf,"%s",answer);
		// printf("center: %s (%s)\n",answer,buf);
		muls.centerSlices = (tolower(answer[0]) == (int)'y');
	}
	// just in case the answer was not exactly 1 or 0:
	// muls.centerSlices = (muls.centerSlices) ? 1 : 0;

	muls.sliceThickness = 0.0;
	if (readparam("slice-thickness:",buf,1)) {
		sscanf(buf,"%g",&(muls.sliceThickness));
		if (readparam("slices:",buf,1)) {
			sscanf(buf,"%d",&(muls.slices));
		}
		else {
			if (muls.cubez >0)
				muls.slices = (int)(muls.cubez/(muls.cellDiv*muls.sliceThickness)+0.99);
			else
				muls.slices = (int)(muls.c/(muls.cellDiv*muls.sliceThickness)+0.99);
		}
		muls.slices += muls.centerSlices;
	}
	else {
		muls.slices = 0;
		if (readparam("slices:",buf,1)) {
			sscanf(buf,"%d",&(muls.slices));
			// muls.slices = (int)(muls.slices*muls.nCellZ/muls.cellDiv);
			if (muls.sliceThickness == 0.0) {
				if ((muls.slices == 1) && (muls.cellDiv == 1)) {
					if (muls.cubez >0)
						muls.sliceThickness = (muls.centerSlices) ? 2.0*muls.cubez/muls.cellDiv : muls.cubez/muls.cellDiv;
					else
						// muls.sliceThickness = (muls.centerSlices) ? 2.0*muls.c/(muls.cellDiv) : muls.c/(muls.cellDiv);
						muls.sliceThickness = (muls.centerSlices) ? 1.0*muls.c/(muls.cellDiv) : muls.c/(muls.cellDiv);
				}
				else {
					if (muls.cubez >0) {
						muls.sliceThickness = muls.cubez/(muls.cellDiv*muls.slices-muls.centerSlices);
					}
					else {
						muls.sliceThickness = muls.c/(muls.cellDiv*muls.slices);
					}
				}
			}
			else {
				muls.cellDiv = (muls.cubez >0) ? (int)ceil(muls.cubez/(muls.slices*muls.sliceThickness)) :
					(int)ceil(muls.c/(muls.slices*muls.sliceThickness));
			if (muls.cellDiv < 1) muls.cellDiv = 1;
			}
		}
	}
	if (muls.slices == 0) {
		if (muls.printLevel > 0) printf("Error: Number of slices = 0\n");
		exit(0);
	}
	/* Find out whether we need to recalculate the potential every time, or not
	*/

	muls.equalDivs = ((!muls.tds)  && (muls.nCellZ % muls.cellDiv == 0) &&
		(fabs(muls.slices*muls.sliceThickness-muls.c/muls.cellDiv) < 1e-5));

	// read the output interval:
	muls.outputInterval = muls.slices;
	if (readparam("slices between outputs:",buf,1)) sscanf(buf,"%d",&(muls.outputInterval));
	if (muls.outputInterval < 1) muls.outputInterval= muls.slices;



	initMuls();
	muls.czOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam("zOffset:",buf,1)) sscanf(buf,"%g",&(muls.czOffset));



	/***********************************************************************
	* Fit the resolution to the wave function array, if not specified different
	*/
	if (muls.resolutionX == 0.0)
		muls.resolutionX = muls.ax / (double)muls.nx;
	if (muls.resolutionY == 0.0)
		muls.resolutionY = muls.by / (double)muls.ny;



	/************************************************************************
	* Optional parameters:
	* determine whether potential periodic or not, etc.:
	*/
	muls.nonPeriodZ = 1;
	muls.nonPeriod = 1;
	if (readparam("periodicXY:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.nonPeriod = (tolower(answer[0]) != (int)'y');
	}
	if (readparam("periodicZ:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.nonPeriodZ = (tolower(answer[0]) != (int)'y'); /* if 'y' -> nonPeriodZ=0 */
	}
	if ((muls.nonPeriodZ == 0) && (muls.cellDiv > 1)) {
		printf("****************************************************************\n"
			"* Warning: cannot use cell divisions >1 and Z-periodic potential\n"
			"* periodicZ = NO\n"
			"****************************************************************\n");
		muls.nonPeriodZ = 1;
	}

	muls.bandlimittrans = 1;
	if (readparam("bandlimit f_trans:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.bandlimittrans = (tolower(answer[0]) == (int)'y');
	}
	muls.readPotential = 0;
	if (readparam("read potential:",buf,1)) {
		sscanf(buf," %s",answer);
		muls.readPotential = (tolower(answer[0]) == (int)'y');
	}
	muls.savePotential = 0;
	if (readparam("save potential:",buf,1)) {
		sscanf(buf," %s",answer);
		muls.savePotential = (tolower(answer[0]) == (int)'y');
	}
	muls.saveTotalPotential = 0;
	if (readparam("save projected potential:",buf,1)) {
		sscanf(buf," %s",answer);
		muls.saveTotalPotential = (tolower(answer[0]) == (int)'y');
	}
	muls.plotPotential = 0;
	if (readparam("plot V(r)*r:",buf,1)) {
		sscanf(buf," %s",answer);
		muls.plotPotential = (tolower(answer[0]) == (int)'y');
	}
	muls.fftpotential = 1;
	if (readparam("one time integration:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.fftpotential = (tolower(answer[0]) == (int)'y');
	}
	muls.potential3D = 1;
	if (readparam("potential3D:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.potential3D = (tolower(answer[0]) == (int)'y');
	}
	muls.avgRuns = 10;
	if (readparam("Runs for averaging:",buf,1))
		sscanf(buf,"%d",&(muls.avgRuns));

	muls.storeSeries = 0;
	if (readparam("Store TDS diffr. patt. series:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.storeSeries = (tolower(answer[0]) == (int)'y');
	}

	if (!muls.tds) muls.avgRuns = 1;

	muls.scanXStart = muls.ax/2.0;
	muls.scanYStart = muls.by/2.0;
	muls.scanXN = 1;
	muls.scanYN = 1;
	muls.scanXStop = muls.scanXStart;
	muls.scanYStop = muls.scanYStart;

  printf("%f \n\n",muls.scanXStart);
	switch (muls.mode) {
		/////////////////////////////////////////////////////////
		// read the position for doing CBED:
	case CBED:
		if (readparam("scan_x_start:",buf,1)) sscanf(buf,"%g",&(muls.scanXStart));
		if (readparam("scan_y_start:",buf,1)) sscanf(buf,"%g",&(muls.scanYStart));
		muls.scanXStop = muls.scanXStart;
		muls.scanYStop = muls.scanYStart;
		break;
		/////////////////////////////////////////////////////////
		// read the position for doing NBED:
	case NBED:
		if (readparam("scan_x_start:", buf, 1)) sscanf(buf, "%g", &(muls.scanXStart));
		if (readparam("scan_y_start:", buf, 1)) sscanf(buf, "%g", &(muls.scanYStart));
		muls.scanXStop = muls.scanXStart;
		muls.scanYStop = muls.scanYStart;
		break;

		/////////////////////////////////////////////////////////
		// Read STEM scanning parameters

	case STEM:
		/* Read in scan coordinates: */
		if (!readparam("scan_x_start:",buf,1)) exit(0);
		sscanf(buf,"%g",&(muls.scanXStart));
		if (!readparam("scan_x_stop:",buf,1)) exit(0);
		sscanf(buf,"%g",&(muls.scanXStop));
		if (!readparam("scan_x_pixels:",buf,1)) exit(0);
		sscanf(buf,"%d",&(muls.scanXN));
		if (!readparam("scan_y_start:",buf,1)) exit(0);
		sscanf(buf,"%g",&(muls.scanYStart));
		if (!readparam("scan_y_stop:",buf,1)) exit(0);
		sscanf(buf,"%g",&(muls.scanYStop));
		if (!readparam("scan_y_pixels:",buf,1)) exit(0);
		sscanf(buf,"%d",&(muls.scanYN));

		if (muls.scanXN < 1) muls.scanXN = 1;
		if (muls.scanYN < 1) muls.scanYN = 1;
		// if (muls.scanXStart > muls.scanXStop) muls.scanXN = 1;

		muls.displayProgInterval = muls.scanYN*muls.scanYN;
		if (readparam("propagation progress interval:",buf,1))
			sscanf(buf,"%d",&(muls.displayProgInterval));
	}
	muls.displayPotCalcInterval = 100000; // RAM: default, but normally read-in by .CFG file in next code fragment
	if ( readparam( "potential progress interval:", buf, 1 ) )
	{
		sscanf( buf, "%d", &(muls.displayPotCalcInterval) );
	}

	/**********************************************************************
	* Read STEM/CBED probe parameters
	*/
	muls.dE_E = 0.0;
	muls.dI_I = 0.0;
	muls.dV_V = 0.0;
	muls.Cc = 0.0;
	if (readparam("dE/E:",buf,1))
		muls.dE_E = atof(buf);
	if (readparam("dI/I:",buf,1))
		muls.dI_I = atof(buf);
	if (readparam("dV/V:",buf,1))
		muls.dV_V = atof(buf);
	if (readparam("Cc:",buf,1))
		muls.Cc = 1e7*atof(buf);


	/* memorize dE_E0, and fill the array of well defined energy deviations */
	dE_E0 = sqrt(muls.dE_E*muls.dE_E+muls.dI_I*muls.dI_I+muls.dV_V*muls.dV_V);
	muls.dE_EArray = (double *)malloc((muls.avgRuns+1)*sizeof(double));
	muls.dE_EArray[0] = 0.0;
	/***********************************************************
	* Statistical gaussian energy spread
	* (takes too long for the statistics to become gaussian)
	*/

	/* for (i = 1;i <= muls.avgRuns*muls.tds; i++) {
	muls.dE_EArray[i] = rangauss(&iseed)*dE_E0;
	printf("dE/E[%d]: %g\n",i,muls.dE_EArray[i]);
	}
	*/

	/**********************************************************
	* quick little fix to calculate gaussian energy distribution
	* without using statistics (better for only few runs)
	*/
	if (muls.printLevel > 0) printf("avgRuns: %d\n",muls.avgRuns);
	// serious bug in Visual C - dy comes out enormous.
	//dy = sqrt((double)pi)/((double)2.0*(double)(muls.avgRuns));
	// using precalculated sqrt(pi):
	dy = 1.772453850905/((double)2.0*(double)(muls.avgRuns));
	dx = pi/((double)(muls.avgRuns+1)*20);
	for (ix=1,x=0,y=0;ix<muls.avgRuns;x+=dx) {
		y += exp(-x*x)*dx;
		if (y>=ix*dy) {
			muls.dE_EArray[ix++] = x*2*dE_E0/pi;
			if (muls.printLevel > 2) printf("dE[%d]: %g eV\n",ix,muls.dE_EArray[ix-1]*muls.v0*1e3);
			if (ix < muls.avgRuns) {
				muls.dE_EArray[ix] = -muls.dE_EArray[ix-1];
				ix ++;
				if (muls.printLevel > 2) printf("dE[%d]: %g eV\n",ix,muls.dE_EArray[ix-1]*muls.v0*1e3);
			}
		}
	}


	if (!readparam("Cs:",buf,1))  exit(0);
	sscanf(buf,"%g",&(muls.Cs)); /* in mm */
	muls.Cs *= 1.0e7; /* convert Cs from mm to Angstroem */

	muls.C5 = 0;
	if (readparam("C5:",buf,1)) {
		sscanf(buf,"%g",&(muls.C5)); /* in mm */
		muls.C5 *= 1.0e7; /* convert C5 from mm to Angstroem */
	}

	/* assume Scherzer defocus as default */
	muls.df0 = -(float)sqrt(1.5*muls.Cs*(wavelength(muls.v0))); /* in A */
	muls.Scherzer = 1;
	if (readparam("defocus:",buf,1)) {
		sscanf(buf,"%s",answer);
		/* if Scherzer defocus */
		if (tolower(answer[0]) == 's') {
			muls.df0 = -(float)sqrt(1.5*muls.Cs*(wavelength(muls.v0)));
			muls.Scherzer = 1;
		}
		else if (tolower(answer[0]) == 'o') {
			muls.df0 = -(float)sqrt(muls.Cs*(wavelength(muls.v0)));
			muls.Scherzer = 2;
		}
		else {
			sscanf(buf,"%g",&(muls.df0)); /* in nm */
			muls.df0 = 10.0*muls.df0;       /* convert defocus to A */
			muls.Scherzer = (-(float)sqrt(1.5*muls.Cs*(wavelength(muls.v0)))==muls.df0);
		}
	}
	// Astigmatism:
	muls.astigMag = 0;
	if (readparam("astigmatism:",buf,1)) sscanf(buf,"%g",&(muls.astigMag));
	// convert to A from nm:
	muls.astigMag = 10.0*muls.astigMag;
	muls.astigAngle = 0;
	if (readparam("astigmatism angle:",buf,1)) sscanf(buf,"%g",&(muls.astigAngle));
	// convert astigAngle from deg to rad:
	muls.astigAngle *= pi/180.0;

	////////////////////////////////////////////////////////
	// read in more aberrations:
	muls.a33 = 0;
	muls.a31 = 0;
	muls.a44 = 0;
	muls.a42 = 0;
	muls.a55 = 0;
	muls.a53 = 0;
	muls.a51 = 0;
	muls.a66 = 0;
	muls.a64 = 0;
	muls.a62 = 0;

	muls.phi33 = 0;
	muls.phi31 = 0;
	muls.phi44 = 0;
	muls.phi42 = 0;
	muls.phi55 = 0;
	muls.phi53 = 0;
	muls.phi51 = 0;
	muls.phi66 = 0;
	muls.phi64 = 0;
	muls.phi62 = 0;

	if (readparam("a_33:",buf,1)) {sscanf(buf,"%g",&(muls.a33)); }
	if (readparam("a_31:",buf,1)) {sscanf(buf,"%g",&(muls.a31)); }
	if (readparam("a_44:",buf,1)) {sscanf(buf,"%g",&(muls.a44)); }
	if (readparam("a_42:",buf,1)) {sscanf(buf,"%g",&(muls.a42)); }
	if (readparam("a_55:",buf,1)) {sscanf(buf,"%g",&(muls.a55)); }
	if (readparam("a_53:",buf,1)) {sscanf(buf,"%g",&(muls.a53)); }
	if (readparam("a_51:",buf,1)) {sscanf(buf,"%g",&(muls.a51)); }
	if (readparam("a_66:",buf,1)) {sscanf(buf,"%g",&(muls.a66)); }
	if (readparam("a_64:",buf,1)) {sscanf(buf,"%g",&(muls.a64)); }
	if (readparam("a_62:",buf,1)) {sscanf(buf,"%g",&(muls.a62)); }

	if (readparam("phi_33:",buf,1)) {sscanf(buf,"%g",&(muls.phi33)); }
	if (readparam("phi_31:",buf,1)) {sscanf(buf,"%g",&(muls.phi31)); }
	if (readparam("phi_44:",buf,1)) {sscanf(buf,"%g",&(muls.phi44)); }
	if (readparam("phi_42:",buf,1)) {sscanf(buf,"%g",&(muls.phi42)); }
	if (readparam("phi_55:",buf,1)) {sscanf(buf,"%g",&(muls.phi55)); }
	if (readparam("phi_53:",buf,1)) {sscanf(buf,"%g",&(muls.phi53)); }
	if (readparam("phi_51:",buf,1)) {sscanf(buf,"%g",&(muls.phi51)); }
	if (readparam("phi_66:",buf,1)) {sscanf(buf,"%g",&(muls.phi66)); }
	if (readparam("phi_64:",buf,1)) {sscanf(buf,"%g",&(muls.phi64)); }
	if (readparam("phi_62:",buf,1)) {sscanf(buf,"%g",&(muls.phi62)); }

	muls.phi33 /= (float)RAD2DEG;
	muls.phi31 /= (float)RAD2DEG;
	muls.phi44 /= (float)RAD2DEG;
	muls.phi42 /= (float)RAD2DEG;
	muls.phi55 /= (float)RAD2DEG;
	muls.phi53 /= (float)RAD2DEG;
	muls.phi51 /= (float)RAD2DEG;
	muls.phi66 /= (float)RAD2DEG;
	muls.phi64 /= (float)RAD2DEG;
	muls.phi62 /= (float)RAD2DEG;


	if (!readparam("alpha:",buf,1)) exit(0);
	sscanf(buf,"%g",&(muls.alpha)); /* in mrad */

	muls.aAIS = 0;  // initialize AIS aperture to 0 A
	if (readparam("AIS aperture:",buf,1))
		sscanf(buf,"%g",&(muls.aAIS)); /* in A */

	///// read beam current and dwell time ///////////////////////////////
	muls.beamCurrent = 1;  // pico Ampere
	muls.dwellTime = 1;    // msec
	if (readparam("beam current:",buf,1)) {
		sscanf(buf,"%g",&(muls.beamCurrent)); /* in pA */
	}
	if (readparam("dwell time:",buf,1)) {
		sscanf(buf,"%g",&(muls.dwellTime)); /* in msec */
	}
	muls.electronScale = muls.beamCurrent*muls.dwellTime*MILLISEC_PICOAMP;
	//////////////////////////////////////////////////////////////////////

	muls.sourceRadius = 0;
	if (readparam("Source Size (diameter):",buf,1))
		muls.sourceRadius = atof(buf)/2.0;

	if (readparam("smooth:",buf,1)) sscanf(buf,"%s",answer);
	muls.ismoth = (tolower(answer[0]) == (int)'y');
	muls.gaussScale = 0.05f;
	muls.gaussFlag = 0;
	if (readparam("gaussian:",buf,1)) {
		sscanf(buf,"%s %g",answer,&(muls.gaussScale));
		muls.gaussFlag = (tolower(answer[0]) == (int)'y');
	}

	/**********************************************************************
	* Parameters for image display and directories, etc.
	*/
	muls.imageGamma = 1.0;
	if (readparam("Display Gamma:",buf,1)) {
		muls.imageGamma = atof(buf);
	}
	muls.showProbe = 0;
	if (readparam("show Probe:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.showProbe = (tolower(answer[0]) == (int)'y');
	}

	// Setup folder for saved data
	sprintf(muls.folder,"data");
	if (readparam("Folder:",buf,1))
		sscanf(buf," %s",muls.folder);

	// search for second '"', in case the folder name is in quotation marks:
	printf( "DEBUG stem3.readFile: folder is %s \n", muls.folder );
	if ( muls.folder[0] == '"' ) {
		strPtr = strchr( buf, '"' );
		strcpy( muls.folder, strPtr + 1 );
		strPtr = strchr( muls.folder, '"' );
		*strPtr = '\0';
	}

	if ( muls.folder[strlen( muls.folder ) - 1] == '/' )
	{
		muls.folder[strlen( muls.folder ) - 1] = '\0';
	}

	// Web update
	muls.webUpdate = 0;
	if (readparam("update Web:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.webUpdate = (tolower(answer[0]) == (int)'y');
	}



	/*  readBeams(parFp); */
	/************************************************************************/
	/* read the different detector configurations                           */
	resetParamFile();
	muls.detectorNum = 0;

	if (muls.mode == STEM)
	{
		int tCount = (int)(ceil((double)((muls.slices * muls.cellDiv) / muls.outputInterval)));

		/* first determine number of detectors */
		while (readparam("detector:",buf,0)) muls.detectorNum++;
		/* now read in the list of detectors: */
		resetParamFile();

		// loop over thickness planes where we're going to record intermediates
		// TODO: is this too costly in terms of memory?  It simplifies the parallelization to
		//       save each of the thicknesses in memory, then save to disk afterwards.
		for (int islice=0; islice<=tCount; islice++)
		{
			std::vector<DetectorPtr> detectors;
			resetParamFile();
			while (readparam("detector:",buf,0)) {
				DetectorPtr det = DetectorPtr(new Detector(muls.scanXN, muls.scanYN,
					(muls.scanXStop-muls.scanXStart)/(float)muls.scanXN,
					(muls.scanYStop-muls.scanYStart)/(float)muls.scanYN));

				sscanf(buf,"%g %g %s %g %g",&(det->rInside),
					&(det->rOutside), det->name, &(det->shiftX),&(det->shiftY));

				/* determine v0 specific k^2 values corresponding to the angles */
				det->k2Inside =
					(float)(sin(det->rInside*0.001)/(wavelength(muls.v0)));
				det->k2Outside =
					(float)(sin(det->rOutside*0.001)/(wavelength(muls.v0)));
				// printf("Detector %d: %f .. %f, lambda = %f (%f)\n",i,muls.detectors[i]->k2Inside,muls.detectors[i]->k2Outside,wavelength(muls.v0),muls.v0);
				/* calculate the squares of the ks */
				det->k2Inside *= det->k2Inside;
				det->k2Outside *= det->k2Outside;
				detectors.push_back(det);
			}
			muls.detectors.push_back(detectors);
		}
	}
	/************************************************************************/

	// in case this file has been written by the tomography function, read the current tilt:
	if (readparam("tomo tilt:",buf,1)) {
		sscanf(buf,"%lf %s",&(muls.tomoTilt),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.tomoTilt *= 1000*pi/180.0;
	}
	/************************************************************************
	* Tomography Parameters:
	***********************************************************************/
	if (muls.mode == TOMO) {
		if (readparam("tomo start:",buf,1)) {
			sscanf(buf,"%lf %s",&(muls.tomoStart),answer); /* in mrad */
			if (tolower(answer[0]) == 'd')
				muls.tomoStart *= 1000*pi/180.0;
		}
		if (readparam("tomo step:",buf,1)) {
			sscanf(buf,"%lf %s",&(muls.tomoStep),answer); /* in mrad */
			if (tolower(answer[0]) == 'd')
				muls.tomoStep *= 1000*pi/180.0;
		}

		if (readparam("tomo count:",buf,1))
			muls.tomoCount = atoi(buf);
		if (readparam("zoom factor:",buf,1))
			sscanf(buf,"%lf",&(muls.zoomFactor));
		if ((muls.tomoStep == 0) && (muls.tomoStep > 1))
			muls.tomoStep = -2.0*muls.tomoStart/(double)(muls.tomoCount - 1);
	}
	/***********************************************************************/


	/*******************************************************************
	* Read in parameters related to the calculation of the projected
	* Potential
	*******************************************************************/
	muls.atomRadius = 5.0;
	if (readparam("atom radius:",buf,1))
		sscanf(buf,"%g",&(muls.atomRadius)); /* in A */
	// why ??????  so that number of subdivisions per slice >= number of fitted points!!!
	/*
	if (muls.atomRadius < muls.sliceThickness)
	muls.atomRadius = muls.sliceThickness;
	*/
	muls.scatFactor = DOYLE_TURNER;
	if (readparam("Structure Factors:",buf,1)) {
		sscanf(buf," %s",answer);
		switch (tolower(answer[0])) {
	case 'w':
		if (tolower(answer[1])=='k') muls.scatFactor = WEICK_KOHL;
		break;
	case 'd':  // DOYLE_TURNER
		muls.scatFactor = DOYLE_TURNER;
		break;
	case 'c':  // CUSTOM - specify k-lookup table and values for all atoms used
		muls.scatFactor = CUSTOM;
		// we already have the kinds of atoms stored in
		// int *muls.Znums and int muls.atomKinds
		readSFactLUT();
		break;
	default:
		muls.scatFactor = DOYLE_TURNER;
		}
	}




	/***************************************************************
	* We now need to determine the size of the potential array,
	* and the offset from its edges in A.  We only need to calculate
	* as much potential as we'll be illuminating later with the
	* electron beam.
	**************************************************************/
	muls.potOffsetX = 0;
	muls.potOffsetY = 0;

	if ((muls.mode == STEM) || (muls.mode == CBED)) {
		/* we are assuming that there is enough atomic position data: */
		muls.potOffsetX = muls.scanXStart - 0.5*muls.nx*muls.resolutionX;
		muls.potOffsetY = muls.scanYStart - 0.5*muls.ny*muls.resolutionY;
		/* adjust scanStop so that it coincides with a full pixel: */
		muls.potNx = (int)((muls.scanXStop-muls.scanXStart)/muls.resolutionX);
		muls.potNy = (int)((muls.scanYStop-muls.scanYStart)/muls.resolutionY);
		muls.scanXStop = muls.scanXStart+muls.resolutionX*muls.potNx;
		muls.scanYStop = muls.scanYStart+muls.resolutionY*muls.potNy;
		muls.potNx+=muls.nx;
		muls.potNy+=muls.ny;
		muls.potSizeX = muls.potNx*muls.resolutionX;
		muls.potSizeY = muls.potNy*muls.resolutionY;
	}
	else {
		muls.potNx = muls.nx;
		muls.potNy = muls.ny;
		muls.potSizeX = muls.potNx*muls.resolutionX;
		muls.potSizeY = muls.potNy*muls.resolutionY;
		muls.potOffsetX = muls.scanXStart - 0.5*muls.potSizeX;
		muls.potOffsetY = muls.scanYStart - 0.5*muls.potSizeY;
	}
	/**************************************************************
	* Check to see if the given scan parameters really fit in cell
	* dimensions:
	*************************************************************/
	if ((muls.scanXN <=0) ||(muls.scanYN <=0)) {
		printf("The number of scan pixels must be >=1\n");
		exit(0);
	}
	if ((muls.scanXStart<0) || (muls.scanYStart<0) ||
		(muls.scanXStop<0) || (muls.scanYStop<0) ||
		(muls.scanXStart>muls.ax) || (muls.scanYStart>muls.by) ||
		(muls.scanXStop>muls.ax) || (muls.scanYStop>muls.by)) {
			printf("Scanning window is outside model dimensions (%g,%g .. %g,%g) [ax = %g, by = %g]!\n",muls.scanXStart,muls.scanYStart,muls.scanXStop,muls.scanYStop,muls.ax,muls.by);
			exit(0);
	}
	/*************************************************************
	* read in the beams we want to plot in the pendeloesung plot
	* Only possible if not in STEM or CBED mode
	*************************************************************/
	muls.lbeams = 0;   /* flag for beam output */
	muls.nbout = 0;    /* number of beams */
	resetParamFile();
	if ((muls.mode != STEM) && (muls.mode != CBED)) {
		if (readparam("Pendelloesung plot:",buf,1)) {
			sscanf(buf,"%s",answer);
			muls.lbeams = (tolower(answer[0]) == (int)'y');
		}
		if (muls.lbeams) {
			while (readparam("beam:",buf,0)) muls.nbout++;
			printf("will record %d beams\n",muls.nbout);
			muls.hbeam = (int*)malloc(muls.nbout*sizeof(int));
			muls.kbeam = (int*)malloc(muls.nbout*sizeof(int));
			/* now read in the list of detectors: */
			resetParamFile();
			for (i=0;i<muls.nbout;i++) {
				if (!readparam("beam:",buf,0)) break;
				muls.hbeam[i] = 0;
				muls.kbeam[i] = 0;
				sscanf(buf,"%d %d",muls.hbeam+i,muls.kbeam+i);
				muls.hbeam[i] *= muls.nCellX;
				muls.kbeam[i] *= muls.nCellY;

				muls.hbeam[i] = (muls.hbeam[i]+muls.nx) % muls.nx;
				muls.kbeam[i] = (muls.kbeam[i]+muls.ny) % muls.ny;
				printf("beam %d [%d %d]\n",i,muls.hbeam[i],muls.kbeam[i]); 			}
		}
	}

	// muls.btiltx = 0;
	// muls.btilty = 0;
	//wave->thickness = 0.0;

	/* TODO: possible breakage here - MCS 2013/04 - made muls.cfgFile be allocated on the struct
	       at runtim - thus this null check doesn't make sense anymore.  Change cfgFile set
	   Old comment:
		if cfgFile != NULL, the program will later write a the atomic config to this file */
	//muls.cfgFile = NULL;
	// RAM Apr14: Fix this broken code by adding an else statement to generate a default cfg-file extension for TDS/tilted specimens
	if ( readparam( "CFG-file:", buf, 1 ) )
	{
		printf( "Debug, MSC breakge for filename: %s", muls.cfgFile );
		sscanf( buf, "%s", muls.cfgFile );
	}
	else
	{
		// RAM: Build a default filename with added _t to indicated tilted or tds crystal.
		strcpy( muls.cfgFile, muls.fileBase );
		strcat( muls.cfgFile, "t" );
		printf( "DEBUG: tilt/tds default filename made = %s \n", muls.cfgFile );
	}

	/* allocate memory for wave function */

	potDimensions[0] = muls.potNx;
	potDimensions[1] = muls.potNy;
#if FLOAT_PRECISION == 1
	muls.trans = complex3Df(muls.slices,muls.potNx,muls.potNy,"trans");
	// printf("allocated trans %d %d %d\n",muls.slices,muls.potNx,muls.potNy);
	muls.fftPlanPotForw = fftwf_plan_many_dft(2,potDimensions, muls.slices,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy, FFTW_FORWARD, fftMeasureFlag);
	muls.fftPlanPotInv = fftwf_plan_many_dft(2,potDimensions, muls.slices,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy, FFTW_BACKWARD, fftMeasureFlag);
#else

	muls.trans = complex3D(muls.slices,muls.potNx,muls.potNy,"trans");
	muls.fftPlanPotForw = fftw_plan_many_dft(2,potDimensions, muls.slices,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy, FFTW_FORWARD, fftMeasureFlag);
	muls.fftPlanPotInv = fftw_plan_many_dft(2,potDimensions, muls.slices,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy, FFTW_BACKWARD, fftMeasureFlag);
#endif

	////////////////////////////////////
	if (muls.printLevel >= 4)
		printf("Memory for transmission function (%d x %d x %d) allocated and plans initiated\n",muls.slices,muls.potNx,muls.potNy);


	// printf("%d %d %d %d\n",muls.nx,muls.ny,sizeof(fftw_complex),(int)(&muls.wave[2][2])-(int)(&muls.wave[2][1]));


} /* end of readFile() */



/************************************************************************
* doTOMO performs a Diffraction Tomography simulation
*
* This routine creates a script file which can then be run by a separate
* command.
* For now this routine will only allow tomography about the y-axis.
* To do anything else one can start with a previously rotated super-cell.
*
* Important parameters: tomoStart, tomoStep, tomoCount, zoomFactor
***********************************************************************/
void doTOMO() {
	double boxXmin=0,boxXmax=0,boxYmin=0,boxYmax=0,boxZmin=0,boxZmax=0;
	double mAx,mBy,mCz;
	int ix,iy,iz,iTheta,i;
	double u[3],**Mm = NULL;
	double theta = 0;
	atom *atoms = NULL;
	char cfgFile[512],stemFile[512],scriptFile[64],diffAnimFile[64];
	FILE *fpScript,*fpDiffAnim;


	Mm = muls.Mm;
	atoms = (atom *)malloc(muls.natom*sizeof(atom));

	boxXmin = boxXmax = muls.ax/2.0;
	boxYmin = boxYmax = muls.by/2.0;
	boxZmin = boxZmax = muls.c/2.0;

	// For all tomography tilt angles:
	for (iTheta = 0;iTheta < muls.tomoCount;iTheta++) {
		theta = muls.tomoStart + iTheta*muls.tomoStep; // theta in mrad
		// Try different corners of the box, and see, how far they poke out.
		for (ix=-1;ix<=1;ix++) for (iy=-1;iy<=1;iy++) for (iz=-1;iz<=1;iz++) {
			// Make center of unit cell rotation center
			u[0]=ix*muls.ax/2; u[1]=iy*muls.by/2.0; u[2]=iz*muls.c/2.0;

			// rotate about y-axis
			rotateVect(u,u,0,theta*1e-3,0);

			// shift origin back to old (0,0,0):
			u[0]+=muls.ax/2; u[1]+=muls.by/2.0; u[2]+=muls.c/2.0;

			boxXmin = boxXmin>u[0] ? u[0] : boxXmin; boxXmax = boxXmax<u[0] ? u[0] : boxXmax;
			boxYmin = boxYmin>u[1] ? u[1] : boxYmin; boxYmax = boxYmax<u[1] ? u[1] : boxYmax;
			boxZmin = boxZmin>u[2] ? u[2] : boxZmin; boxZmax = boxZmax<u[2] ? u[2] : boxZmax;

		}
	} /* for iTheta ... */

	// find max. box size:
	boxXmax -= boxXmin;
	boxYmax -= boxYmin;
	boxZmax -= boxZmin;
	printf("Minimum box size for tomography tilt series: %g x %g x %gA, zoom Factor: %g\n",
		boxXmax,boxYmax,boxZmax,muls.zoomFactor);
	boxXmax /= muls.zoomFactor;
	boxYmax = boxXmax*muls.by/muls.ax;

	// boxMin will now be boxCenter:
	boxXmin = 0.5*boxXmax;
	boxYmin = 0.5*boxYmax;
	boxZmin = 0.5*boxZmax;

	// We have to save the original unit cell dimensions
	mAx = muls.ax; mBy = muls.by; mCz = muls.c;
	muls.ax=boxXmax; muls.by=boxYmax; muls.c=boxZmax;

	printf("Will use box sizes: %g x %g x %gA (kept original aspect ratio). \n"
		"Writing structure files now, please wait ...\n",
		boxXmax,boxYmax,boxZmax);

	// open the script file and write the stem instructions in there
	sprintf(scriptFile,"%s/run_tomo",muls.folder);
	if ((fpScript=fopen(scriptFile,"w")) == NULL) {
		printf("doTOMO: unable to open scriptFile %s for writing\n",scriptFile);
		exit(0);
	}
	fprintf(fpScript,"#!/bin/bash\n\n");

	// open the diffraction animation file
	sprintf(diffAnimFile,"%s/diff_anim",muls.folder);
	if ((fpDiffAnim=fopen(diffAnimFile,"w")) == NULL) {
		printf("doTOMO: unable to open diffraction animation %s for writing\n",diffAnimFile);
		exit(0);
	}
	fprintf(fpDiffAnim,"#!/bin/bash\n\n");


	// We will now rotate all the unit cell and write the answer to
	// separate structure output files.
	for (iTheta = 0;iTheta < muls.tomoCount;iTheta++) {
		theta = muls.tomoStart + iTheta*muls.tomoStep; // theta in mrad
		muls.tomoTilt = theta;

		// rotate the structure and write result to local atom array
		for(i=0;i<(muls.natom);i++) {
			u[0] = muls.atoms[i].x - mAx/2.0;
			u[1] = muls.atoms[i].y - mBy/2.0;
			u[2] = muls.atoms[i].z - mCz/2.0;
			rotateVect(u,u,0,theta*1e-3,0);
			atoms[i].x = u[0]+boxXmin;
			atoms[i].y = u[1]+boxYmin;
			atoms[i].z = u[2]+boxZmin;
			atoms[i].Znum = muls.atoms[i].Znum;
			atoms[i].occ = muls.atoms[i].occ;
			atoms[i].dw = muls.atoms[i].dw;
		}

		sprintf(cfgFile,"%s/tomo_%dmrad.cfg",muls.folder,(int)theta);
		sprintf(stemFile,"%s/tomo_%dmrad.dat",muls.folder,(int)theta);
		printf("Writing file %s | ",cfgFile);
		writeCFG(atoms,muls.natom,cfgFile,&muls);
		sprintf(cfgFile,"tomo_%dmrad.cfg",(int)theta);
		writeSTEMinput(stemFile,cfgFile,&muls);

		// add to script files:
		fprintf(fpScript,"stem tomo_%dmrad.dat\n",(int)theta);
	}
	muls.ax = mAx; muls.by = mBy; muls.c = mCz;
	sprintf(stemFile,"copy fparams.dat %s/",muls.folder);
	system(stemFile);

	// close script files again
	fprintf(fpDiffAnim,"convert -delay 20 diff*.jpg diff.gif\n");
	fclose(fpScript);
	fclose(fpDiffAnim);
	sprintf(stemFile,"chmod +x %s",scriptFile);
	system(stemFile);
	sprintf(stemFile,"chmod +x %s",diffAnimFile);
	system(stemFile);

	exit(0);
}

/************************************************************************
* doMSCBED performs a super-CBED calculation based on the MS algorithm
* including phonons.
*
***********************************************************************/
void doMSCBED() {


}

/************************************************************************
* doNBED performs a NBED calculation
*  -- added by Robert A. McLeod 04 April 2014
*  -- designed to read an .IMG file as input and use that instead of a STEM probe
*    otherwise code should be very similar to CBED.
*	From a programming perspective it would be better to update CBED to be more general, but since I share this code and don't own it, not so great.
*
***********************************************************************/

void doNBED()
{
	// RAM: image reading is found in imagelib_fftw3
	// Basic idea is to use the writeIMG function in Matlab and use it to pass in a complex-value probe to stem3.exe which can be rastered like a CBED probe
	// RAM: All we really need to do is initialize the WavePtr wave with a filename found by readFile() above ^^
	// FIXME: With phonon calcs the original WavePtr is not being reloaded for each series, which produces some funky results.


	int ix, iy, i, pCount, result;
	FILE *avgFp, *fpWave, *fpPos, *fpNBED, *fpTest = 0;
	double timer, timerTot;
	double probeCenterX, probeCenterY, probeOffsetX, probeOffsetY;
	char buf[BUF_LEN], avgName[32], systStr[64];
	real t = 0;
	real **avgPendelloesung = NULL;
	int oldMulsRepeat1 = 1;
	int oldMulsRepeat2 = 1;
	long iseed = 0;
	WavePtr wave = WavePtr(new WAVEFUNC(muls.nx, muls.ny, muls.resolutionX, muls.resolutionY));
	ImageIOPtr imageIO = ImageIOPtr(new CImageIO(muls.nx, muls.ny, t, muls.resolutionX, muls.resolutionY));
	std::vector<double> params(2);

	//printf("Debug doNBED: wavefile: %s\n",muls.fileWaveIn);

	// Try and test
	//fpWave = fopen(muls.fileWaveIn, "r+");
	//printf("Debug doNBED1 \n");
	//if( fpWave != NULL ) fclose(fpWave);
	//printf("Debug doNBED2 \n");
	// RAM: read-in the entry wave .IMG
	wave->ReadWave(muls.fileWaveIn);

	//printf("Debug doNBED3 \n");

	muls.chisq = std::vector<double>(muls.avgRuns);

	if (iseed == 0) iseed = -(long)time(NULL);

	if (muls.lbeams) {
		muls.pendelloesung = NULL;
		if (avgPendelloesung == NULL) {
			avgPendelloesung = float2D(muls.nbout,
				muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
				"pendelloesung");
		}
	}
	probeCenterX = muls.scanXStart;
	probeCenterY = muls.scanYStart;

	timerTot = 0; /* cputim();*/
	displayProgress(-1);

	for (muls.avgCount = 0; muls.avgCount < muls.avgRuns; muls.avgCount++) {
		muls.totalSliceCount = 0;
		pCount = 0;
		/* make sure we start at the beginning of the file
		so we won't miss any line that contains a sequence,
		because we will not do any EOF wrapping
		*/
		resetParamFile();

		/* probe(&muls,xpos,ypos); */
		/* make incident probe wave function with probe exactly in the center */
		/* if the potential array is not big enough, the probe can
		* then also be adjusted, so that it is off-center
		*/

		probeOffsetX = muls.sourceRadius*gasdev(&iseed)*SQRT_2;
		probeOffsetY = muls.sourceRadius*gasdev(&iseed)*SQRT_2;
		muls.scanXStart = probeCenterX + probeOffsetX;
		muls.scanYStart = probeCenterY + probeOffsetY;


		// RAM: This is where I need to remove the existing STEM probe and in the future add some circular shift and crop?
		// RAM: made a new function in stemlib, probeShiftAndCrop(&muls, wave, muls.scanXStart - muls.potOffsetX, muls.scanYStart - muls.potOffsetY)
		// probe(&muls, wave, muls.scanXStart - muls.potOffsetX, muls.scanYStart - muls.potOffsetY);
		probeShiftAndCrop(&muls, wave, muls.scanXStart - muls.potOffsetX, muls.scanYStart - muls.potOffsetY, muls.nx, muls.ny);
		// RAM: function is empty for now (returns wave)

		if (muls.saveLevel > 2)
		{

			sprintf(systStr, "%s/wave_probe.img", muls.folder);
			wave->WriteWave(systStr);
		}
		// printf("Probe: (%g, %g)\n",muls.scanXStart,muls.scanYStart);


		if (muls.sourceRadius > 0) {
			if (muls.avgCount == 0) fpPos = fopen("probepos.dat", "w");
			else fpPos = fopen("probepos.dat", "a");
			if (fpPos == NULL) {
				printf("Was unable to open file probepos.dat for writing\n");
			}
			else {
				fprintf(fpPos, "%g %g\n", muls.scanXStart, muls.scanYStart);
				fclose(fpPos);
			}
		}

		if ((muls.showProbe) && (muls.avgCount == 0)) {
#ifndef WIN32
			//probePlot(&muls);
			sprintf(buf, "ee %s/probePlot_0.jpg &", muls.folder);
			system(buf);
#endif
		}
		//muls.nslic0 = 0;

		result = readparam("sequence: ", buf, 0);
		while (result) {
			if (((buf[0] < 'a') || (buf[0] > 'z')) &&
				((buf[0] < '1') || (buf[0] > '9')) &&
				((buf[0] < 'A') || (buf[0] > 'Z'))) {
				// printf("Stacking sequence: %s\n",buf);
				printf("Can only work with old stacking sequence\n");
				break;
			}
			muls.mulsRepeat1 = 1;
			muls.mulsRepeat2 = 1;
			sscanf(buf, "%d %d", &muls.mulsRepeat1, &muls.mulsRepeat2);
			for (i = 0; i<(int)strlen(buf); i++) buf[i] = 0;
			if (muls.mulsRepeat2 < 1) muls.mulsRepeat2 = 1;
			sprintf(muls.cin2, "%d", muls.mulsRepeat1);


			/***********************************************************
			* make sure we have enough memory for the pendelloesung plot
			*/
			if ((muls.lbeams) &&
				((oldMulsRepeat1 != muls.mulsRepeat1) ||
				(oldMulsRepeat2 != muls.mulsRepeat2))) {
				oldMulsRepeat1 = muls.mulsRepeat1;
				oldMulsRepeat2 = muls.mulsRepeat2;
				if (muls.pendelloesung != NULL)
					free(muls.pendelloesung[0]);
				free(avgPendelloesung[0]);
				muls.pendelloesung = NULL;
				avgPendelloesung = float2D(muls.nbout,
					muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
					"pendelloesung");
			}
			/*********************************************************/

			// printf("Stacking sequence: %s\n",buf);


			/************************************************
			*make3DSlicesFFT(&muls,muls.slices,atomPosFile,NULL);
			*exit(0);
			************************************************/
			if (muls.equalDivs) {
				if (muls.scatFactor == CUSTOM)
					make3DSlicesFT(&muls);
				else
					make3DSlices(&muls, muls.slices, muls.atomPosFile, NULL);
				initSTEMSlices(&muls, muls.slices);
			}

			muls.saveFlag = 0;
			/****************************************
			* do the (small) loop
			*****************************************/
			for (pCount = 0; pCount<muls.mulsRepeat2*muls.cellDiv; pCount++) {

				// printf("pCount = %d\n",pCount);
				/*******************************************************
				* build the potential slices from atomic configuration
				******************************************************/
				if (!muls.equalDivs) {
					if (muls.scatFactor == CUSTOM)
						make3DSlicesFT(&muls);
					else
						make3DSlices(&muls, muls.slices, muls.atomPosFile, NULL);
					initSTEMSlices(&muls, muls.slices);
				}

				timer = cputim();
				// what probe should runMulsSTEM use here?
				runMulsSTEM(&muls, wave);

				printf("Thickness: %gA, int.=%g, time: %gsec\n",
					wave->thickness, wave->intIntensity, cputim() - timer);

				/***************** Only if Save level > 2: ****************/
				if ((muls.avgCount == 0) && (muls.saveLevel > 2)) {
					sprintf(systStr, "%s/wave_final.img", muls.folder);
					wave->WriteWave(systStr);
				}
#ifdef VIB_IMAGE_TEST_CBED
				wave->WriteWave(sysStr)
#endif
					muls.totalSliceCount += muls.slices;

			} // end of for pCount = 0...
			result = readparam("sequence: ", buf, 0);
		}
		/*    printf("Total CPU time = %f sec.\n", cputim()-timerTot ); */

		sprintf(avgName, "%s/diff.img", muls.folder);
		// TODO: Why are we reading in a DP at this point?  Do we have one yet?
		//     What happens if it isn't there?
		// RAM: Assume these are Michael's comments above, this crashes because the diffraction pattern isn't there yet.  Provide if/else fix
		fpTest = fopen( avgName, "rb" );
		if ( fpTest != NULL )
		{
			printf( "DEBUG doNBED: Found filename %s \n", avgName );
			fclose( fpTest ); // close file handle before opening it again in data_container.cpp
			wave->ReadDiffPat( avgName );
		}
		else
		{
			printf( "DEBUG doNBED: Did not find filename %s \n", avgName );
			wave->WriteDiffPat( avgName, "DEBUG Wave intensity" );
		}
			//fpWave = fopen(muls.fileWaveIn, "r+");
			//printf("Debug doNBED1 \n");
			//if( fpWave != NULL ) fclose(fpWave);



		if (muls.avgCount == 0) {
			memcpy((void *)wave->avgArray[0], (void *)wave->diffpat[0],
				(size_t)(muls.nx*muls.ny*sizeof(float_tt)));
			/* move the averaged (raw data) file to the target directory as well */
			sprintf(avgName, "%s/diffAvg_%d.img", muls.folder, muls.avgCount + 1);
			sprintf(systStr, "mv %s/diff.img %s", muls.folder, avgName);
			system(systStr);
			if (muls.lbeams) {
				for (iy = 0; iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv; iy++) {
					for (ix = 0; ix<muls.nbout; ix++) {
						avgPendelloesung[ix][iy] = muls.pendelloesung[ix][iy];
					}
				}
			}
		} // of if muls.avgCount == 0 ...
		else {
			muls.chisq[muls.avgCount - 1] = 0.0;
			for (ix = 0; ix<muls.nx; ix++) for (iy = 0; iy<muls.ny; iy++) {
				t = ((real)muls.avgCount*wave->avgArray[ix][iy] +
					wave->diffpat[ix][iy]) / ((real)(muls.avgCount + 1));
				muls.chisq[muls.avgCount - 1] += (wave->avgArray[ix][iy] - t)*(wave->avgArray[ix][iy] - t);
				wave->avgArray[ix][iy] = t;

			}
			muls.chisq[muls.avgCount - 1] = muls.chisq[muls.avgCount - 1] / (double)(muls.nx*muls.ny);
			sprintf(avgName, "%s/diffAvg_%d.img", muls.folder, muls.avgCount + 1);
			params[0] = muls.tomoTilt;
			params[1] = 1.0 / wavelength(muls.v0);
			wave->WriteAvgArray(avgName, "Averaged Diffraction pattern, unit: 1/A", params);

			muls.storeSeries = 1;
			if (muls.saveLevel == 0)	muls.storeSeries = 0;
			else if (muls.avgCount % muls.saveLevel != 0) muls.storeSeries = 0;

			if (muls.storeSeries == 0) {
				// printf("Removing old file \n");
				sprintf(avgName, "%s/diffAvg_%d.img", muls.folder, muls.avgCount);
				sprintf(systStr, "rm %s", avgName);
				system(systStr);
			}

			/* write the data to a file */
			if (muls.saveFlag >-1) {
				sprintf(systStr, "%s/avgresults.dat", muls.folder);
				if ((avgFp = fopen(systStr, "w")) == NULL)
					printf("Sorry, could not open data file for averaging\n");
				else {
					for (ix = 0; ix<muls.avgCount; ix++) {
						fprintf(avgFp, "%d %g\n", ix + 1, muls.chisq[ix]);
					}
					fclose(avgFp);
				}
			}
			/*************************************************************/

			/***********************************************************
			* Average over the pendelloesung plot as well
			*/
			if (muls.lbeams) {
				for (iy = 0; iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv; iy++) {
					for (ix = 0; ix<muls.nbout; ix++) {
						avgPendelloesung[ix][iy] =
							((real)muls.avgCount*avgPendelloesung[ix][iy] +
							muls.pendelloesung[ix][iy]) / (real)(muls.avgCount + 1);
					}
				}
			}
		} /* else ... if avgCount was greater than 0 */

		if (muls.lbeams) {
			/**************************************************************
			* The diffraction spot intensities of the selected
			* diffraction spots are now stored in the 2 dimensional array
			* muls.pendelloesung[beam][slice].
			* We can write the array to a file and display it, just for
			* demonstration purposes
			*************************************************************/
			sprintf(systStr, "%s/pendelloesung.dat", muls.folder);
			if ( (fpNBED = fopen( systStr, "w" )) != NULL ) {
				printf("Writing Pendelloesung data\n");
				for (iy = 0; iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv; iy++) {
					/* write the thicknes in the first column of the file */
					fprintf( fpNBED, "%g", iy*muls.c / ((float)(muls.slices*muls.cellDiv)) );
					/* write the beam intensities in the following columns */
					for (ix = 0; ix<muls.nbout; ix++) {
						fprintf( fpNBED, "\t%g", avgPendelloesung[ix][iy] );
					}
					/* close the line, and start a new one for the next set of
					* intensities
					*/
					fprintf( fpNBED, "\n" );
				}
				fclose( fpNBED );
			}
			else {
				printf("Could not open file for pendelloesung plot\n");
			}
		} /* end of if lbemas ... */
		displayProgress(1);
	} /* end of for muls.avgCount=0.. */
	//delete(wave);
}
/************************************************************************
* End of doNBED()
***********************************************************************/

/************************************************************************
* doCBED performs a CBED calculation
*
***********************************************************************/

void doCBED() {
	int ix,iy,i,pCount,result;
	FILE *avgFp, *fpCBED, *fpPos = 0, *fpTest = 0;
	double timer,timerTot;
	double probeCenterX,probeCenterY,probeOffsetX,probeOffsetY;
	char buf[BUF_LEN],avgName[32],systStr[64];
	real t=0;
	real **avgPendelloesung = NULL;
	int oldMulsRepeat1 = 1;
	int oldMulsRepeat2 = 1;
	long iseed=0;
	WavePtr wave = WavePtr(new WAVEFUNC(muls.nx,muls.ny, muls.resolutionX, muls.resolutionY));
	ImageIOPtr imageIO = ImageIOPtr(new CImageIO(muls.nx, muls.ny, t, muls.resolutionX, muls.resolutionY));
	std::vector<double> params(2);

	muls.chisq = std::vector<double>(muls.avgRuns);

	if (iseed == 0) iseed = -(long) time( NULL );

	if (muls.lbeams) {
		muls.pendelloesung = NULL;
		if (avgPendelloesung == NULL) {
			avgPendelloesung = float2D(muls.nbout,
				muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
				"pendelloesung");
		}
	}
	probeCenterX = muls.scanXStart;
	probeCenterY = muls.scanYStart;

	timerTot = 0; /* cputim();*/
	displayProgress(-1);

	for (muls.avgCount = 0;muls.avgCount < muls.avgRuns;muls.avgCount++) {
		muls.totalSliceCount = 0;
		pCount = 0;
		/* make sure we start at the beginning of the file
		so we won't miss any line that contains a sequence,
		because we will not do any EOF wrapping
		*/
		resetParamFile();

		/* probe(&muls,xpos,ypos); */
		/* make incident probe wave function with probe exactly in the center */
		/* if the potential array is not big enough, the probe can
		* then also be adjusted, so that it is off-center
		*/

		probeOffsetX = muls.sourceRadius*gasdev(&iseed)*SQRT_2;
		probeOffsetY = muls.sourceRadius*gasdev(&iseed)*SQRT_2;
		muls.scanXStart = probeCenterX+probeOffsetX;
		muls.scanYStart = probeCenterY+probeOffsetY;
		probe(&muls, wave,muls.scanXStart-muls.potOffsetX,muls.scanYStart-muls.potOffsetY);
		if (muls.saveLevel > 2) {
			sprintf(systStr,"%s/wave_probe.img",muls.folder);
			wave->WriteWave(systStr);
		}
		// printf("Probe: (%g, %g)\n",muls.scanXStart,muls.scanYStart);
		/*****************************************************************
		* For debugging only!!!
		*
		muls->WriteWave("probe.img")
		*****************************************************************/


		if (muls.sourceRadius > 0) {
			if (muls.avgCount == 0) fpPos = fopen("probepos.dat","w");
			else fpPos = fopen("probepos.dat","a");
			if (fpPos == NULL) {
				printf("Was unable to open file probepos.dat for writing\n");
			}
			else {
				fprintf(fpPos,"%g %g\n",muls.scanXStart,muls.scanYStart);
				fclose(fpPos);
			}
		}

		if ((muls.showProbe) && (muls.avgCount == 0)) {
#ifndef WIN32
			//probePlot(&muls);
			sprintf(buf,"ee %s/probePlot_0.jpg &",muls.folder);
			system(buf);
#endif
		}
		//muls.nslic0 = 0;

		result = readparam("sequence: ",buf,0);
		while (result) {
			if (((buf[0] < 'a') || (buf[0] > 'z')) &&
				((buf[0] < '1') || (buf[0] > '9')) &&
				((buf[0] < 'A') || (buf[0] > 'Z'))) {
					// printf("Stacking sequence: %s\n",buf);
					printf("Can only work with old stacking sequence\n");
					break;
			}
			muls.mulsRepeat1 = 1;
			muls.mulsRepeat2 = 1;
			sscanf(buf,"%d %d",&muls.mulsRepeat1,&muls.mulsRepeat2);
			for (i=0;i<(int)strlen(buf);i++) buf[i] = 0;
			if (muls.mulsRepeat2 < 1) muls.mulsRepeat2 = 1;
			sprintf(muls.cin2,"%d",muls.mulsRepeat1);


			/***********************************************************
			* make sure we have enough memory for the pendelloesung plot
			*/
			if ((muls.lbeams) &&
				((oldMulsRepeat1 !=muls.mulsRepeat1) ||
				(oldMulsRepeat2 !=muls.mulsRepeat2))) {
					oldMulsRepeat1 = muls.mulsRepeat1;
					oldMulsRepeat2 = muls.mulsRepeat2;
					if (muls.pendelloesung != NULL)
						free(muls.pendelloesung[0]);
					free(avgPendelloesung[0]);
					muls.pendelloesung = NULL;
					avgPendelloesung = float2D(muls.nbout,
						muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
						"pendelloesung");
			}
			/*********************************************************/

			// printf("Stacking sequence: %s\n",buf);


			/************************************************
			*make3DSlicesFFT(&muls,muls.slices,atomPosFile,NULL);
			*exit(0);
			************************************************/
			if (muls.equalDivs) {
				if (muls.scatFactor == CUSTOM)
					make3DSlicesFT(&muls);
				else
					make3DSlices(&muls,muls.slices,muls.atomPosFile,NULL);
				initSTEMSlices(&muls,muls.slices);
			}

			muls.saveFlag = 0;
			/****************************************
			* do the (small) loop
			*****************************************/
			for (pCount = 0;pCount<muls.mulsRepeat2*muls.cellDiv;pCount++) {

				// printf("pCount = %d\n",pCount);
				/*******************************************************
				* build the potential slices from atomic configuration
				******************************************************/
				if (!muls.equalDivs) {
					if (muls.scatFactor == CUSTOM)
						make3DSlicesFT(&muls);
					else
						make3DSlices(&muls,muls.slices,muls.atomPosFile,NULL);
					initSTEMSlices(&muls,muls.slices);
				}

				timer = cputim();
				// what probe should runMulsSTEM use here?
				runMulsSTEM(&muls,wave);

				printf("Thickness: %gA, int.=%g, time: %gsec\n",
					wave->thickness,wave->intIntensity,cputim()-timer);

				/***************** Only if Save level > 2: ****************/
				if ((muls.avgCount == 0) && (muls.saveLevel > 2)) {
					sprintf(systStr,"%s/wave_final.img",muls.folder);
					wave->WriteWave(systStr);
				}
#ifdef VIB_IMAGE_TEST_CBED
				wave->WriteWave(sysStr)
#endif
				muls.totalSliceCount += muls.slices;

			} // end of for pCount = 0...
			result = readparam("sequence: ",buf,0);
		}
		/*    printf("Total CPU time = %f sec.\n", cputim()-timerTot ); */

		sprintf(avgName,"%s/diff.img",muls.folder);
		// TODO: Why are we reading in a DP at this point?  Do we have one yet?
		//     What happens if it isn't there?
		// RAM: fix to CBED code as well.  In the future doCBED and doNBED should be merged (pending C++ refactor)
		fpTest = fopen(avgName, "rb");
		if (fpTest != NULL)
		{
			printf("DEBUG doCBED: Found filename %s \n", avgName);
			fclose(fpTest); // close file handle before opening it again in data_container.cpp
			wave->ReadDiffPat(avgName);
		}
		else
		{
			printf("DEBUG doCBED: Did not find filename %s \n", avgName);
			wave->WriteDiffPat(avgName, "DEBUG Wave intensity");
		}
		// RAM: old code
		//wave->ReadDiffPat(avgName);

		if (muls.avgCount == 0) {
			memcpy((void *)wave->avgArray[0],(void *)wave->diffpat[0],
				(size_t)(muls.nx*muls.ny*sizeof(float_tt)));
			/* move the averaged (raw data) file to the target directory as well */
			sprintf(avgName,"%s/diffAvg_%d.img",muls.folder,muls.avgCount+1);
			sprintf(systStr,"mv %s/diff.img %s",muls.folder,avgName);
			system(systStr);
			if (muls.lbeams) {
				for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
					for (ix=0;ix<muls.nbout;ix++) {
						avgPendelloesung[ix][iy] = muls.pendelloesung[ix][iy];
					}
				}
			}
		} // of if muls.avgCount == 0 ...
		else {
			muls.chisq[muls.avgCount-1] = 0.0;
			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				t = ((real)muls.avgCount*wave->avgArray[ix][iy]+
					wave->diffpat[ix][iy])/((real)(muls.avgCount+1));
				muls.chisq[muls.avgCount-1] += (wave->avgArray[ix][iy]-t)*(wave->avgArray[ix][iy]-t);
				wave->avgArray[ix][iy] = t;

			}
			muls.chisq[muls.avgCount-1] = muls.chisq[muls.avgCount-1]/(double)(muls.nx*muls.ny);
			sprintf(avgName,"%s/diffAvg_%d.img",muls.folder,muls.avgCount+1);
			params[0] = muls.tomoTilt;
			params[1] = 1.0/wavelength(muls.v0);
			wave->WriteAvgArray(avgName,"Averaged Diffraction pattern, unit: 1/A",params);

			muls.storeSeries = 1;
			if (muls.saveLevel == 0)	muls.storeSeries = 0;
			else if (muls.avgCount % muls.saveLevel != 0) muls.storeSeries = 0;

			if (muls.storeSeries == 0) {
				// printf("Removing old file \n");
				sprintf(avgName,"%s/diffAvg_%d.img",muls.folder,muls.avgCount);
				sprintf(systStr,"rm %s",avgName);
				system(systStr);
			}

			/* write the data to a file */
			if (muls.saveFlag >-1) {
				sprintf(systStr,"%s/avgresults.dat",muls.folder);
				if ((avgFp = fopen(systStr,"w")) == NULL )
					printf("Sorry, could not open data file for averaging\n");
				else {
					for (ix =0;ix<muls.avgCount;ix++) {
						fprintf(avgFp,"%d %g\n",ix+1,muls.chisq[ix]);
					}
					fclose(avgFp);
				}
			}
			/*************************************************************/

			/***********************************************************
			* Average over the pendelloesung plot as well
			*/
			if (muls.lbeams) {
				for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
					for (ix=0;ix<muls.nbout;ix++) {
						avgPendelloesung[ix][iy] =
							((real)muls.avgCount*avgPendelloesung[ix][iy]+
							muls.pendelloesung[ix][iy])/(real)(muls.avgCount+1);
					}
				}
			}
		} /* else ... if avgCount was greater than 0 */

		if (muls.lbeams) {
			/**************************************************************
			* The diffraction spot intensities of the selected
			* diffraction spots are now stored in the 2 dimensional array
			* muls.pendelloesung[beam][slice].
			* We can write the array to a file and display it, just for
			* demonstration purposes
			*************************************************************/
			sprintf(systStr,"%s/pendelloesung.dat",muls.folder);
			if ( (fpCBED = fopen( systStr, "w" )) != NULL ) {
				printf("Writing Pendelloesung data\n");
				for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
					/* write the thicknes in the first column of the file */
					fprintf( fpCBED, "%g", iy*muls.c / ((float)(muls.slices*muls.cellDiv)) );
					/* write the beam intensities in the following columns */
					for (ix=0;ix<muls.nbout;ix++) {
						fprintf( fpCBED, "\t%g", avgPendelloesung[ix][iy] );
					}
					/* close the line, and start a new one for the next set of
					* intensities
					*/
					fprintf( fpCBED, "\n" );
				}
				fclose( fpCBED );
			}
			else {
				printf("Could not open file for pendelloesung plot\n");
			}
		} /* end of if lbemas ... */
		displayProgress(1);
	} /* end of for muls.avgCount=0.. */
	//delete(wave);
}
/************************************************************************
* End of doCBED()
***********************************************************************/

/************************************************************************
* doTEM performs a TEM calculation (with through focus reconstruction)
*
***********************************************************************/

void doTEM() {
	const double pi=3.1415926535897;
	int ix,iy,i,pCount,result;
	FILE *avgFp,*fpTEM; // *fpPos=0;
	double timer,timerTot;
	double x,y,ktx,kty;
	char buf[BUF_LEN],avgName[256],systStr[512];
	char *comment;
	real t;
	real **avgPendelloesung = NULL;
	int oldMulsRepeat1 = 1;
	int oldMulsRepeat2 = 1;
	long iseed=0;
	std::vector<double> params;
	WavePtr wave = WavePtr(new WAVEFUNC(muls.nx,muls.ny,muls.resolutionX,muls.resolutionY));
	fftwf_complex **imageWave = NULL;

	if (iseed == 0) iseed = -(long) time( NULL );

	muls.chisq=std::vector<double>(muls.avgRuns);

	if (muls.lbeams) {
		muls.pendelloesung = NULL;
		if (avgPendelloesung == NULL) {
			avgPendelloesung = float2D(muls.nbout,
				muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
				"pendelloesung");
		}
	}

	timerTot = 0; /* cputim();*/
	displayProgress(-1);
	for (muls.avgCount = 0;muls.avgCount < muls.avgRuns;muls.avgCount++) {
		muls.totalSliceCount = 0;

		pCount = 0;
		resetParamFile();

		/* make incident probe wave function with probe exactly in the center */
		/* if the potential array is not big enough, the probe can
		* then also be adjusted, so that it is off-center
		*/

		// muls.scanXStart = muls.nx/2*muls.resolutionX+muls.sourceRadius*gasdev(&iseed)*sqrt(2);
		// muls.scanYStart = muls.ny/2*muls.resolutionY+muls.sourceRadius*gasdev(&iseed)*sqrt(2);
		// probe(&muls,muls.scanXStart,muls.scanYStart);

		//muls.nslic0 = 0;
		// produce an incident plane wave:
		if ((muls.btiltx == 0) && (muls.btilty == 0)) {
			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				wave->wave[ix][iy][0] = 1;	wave->wave[ix][iy][1] = 0;
			}
		}
		else {
			// produce a tilted wave function (btiltx,btilty):
			ktx = 2.0*pi*sin(muls.btiltx)/wavelength(muls.v0);
			kty = 2.0*pi*sin(muls.btilty)/wavelength(muls.v0);
			for (ix=0;ix<muls.nx;ix++) {
				x = muls.resolutionX*(ix-muls.nx/2);
				for (iy=0;iy<muls.ny;iy++) {
					y = muls.resolutionY*(ix-muls.nx/2);
					wave->wave[ix][iy][0] = (float)cos(ktx*x+kty*y);
					wave->wave[ix][iy][1] = (float)sin(ktx*x+kty*y);
				}
			}
		}

		result = readparam("sequence: ",buf,0);
		while (result) {
			if (((buf[0] < 'a') || (buf[0] > 'z')) &&
				((buf[0] < '1') || (buf[0] > '9')) &&
				((buf[0] < 'A') || (buf[0] > 'Z'))) {
					// printf("Stacking sequence: %s\n",buf);
					printf("Can only work with old stacking sequence\n");
					break;
			}

			muls.mulsRepeat1 = 1;
			muls.mulsRepeat2 = 1;
			sscanf(buf,"%d %d",&muls.mulsRepeat1,&muls.mulsRepeat2);
			for (i=0;i<(int)strlen(buf);i++)
				buf[i] = 0;
			if (muls.mulsRepeat2 < 1) muls.mulsRepeat2 = 1;
			sprintf(muls.cin2,"%d",muls.mulsRepeat1);
			/***********************************************************
			* make sure we have enough memory for the pendelloesung plot
			*/
			if ((muls.lbeams) &&
				((oldMulsRepeat1 !=muls.mulsRepeat1) ||
				(oldMulsRepeat2 !=muls.mulsRepeat2))) {
					oldMulsRepeat1 = muls.mulsRepeat1;
					oldMulsRepeat2 = muls.mulsRepeat2;
					if (muls.pendelloesung != NULL)  free(muls.pendelloesung[0]);
					free(avgPendelloesung[0]);
					muls.pendelloesung = NULL;
					avgPendelloesung = float2D(muls.nbout,
						muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
						"pendelloesung");
			}
			/*********************************************************/

			// printf("Stacking sequence: %s\n",buf);


			/************************************************
			*make3DSlicesFFT(&muls,muls.slices,atomPosFile,NULL);
			*exit(0);
			************************************************/
      printf("slices: %d \n",muls.slices);

			if (muls.equalDivs) {
				if (muls.printLevel > 1) printf("found equal unit cell divisions\n");
				make3DSlices(&muls,muls.slices,muls.atomPosFile,NULL);
				initSTEMSlices(&muls,muls.slices);
			}

			muls.saveFlag = 0;
			/****************************************
			* do the (small) loop through the slabs
			*****************************************/
			for (pCount = 0;pCount<muls.mulsRepeat2*muls.cellDiv;pCount++) {

				/*******************************************************
				* build the potential slices from atomic configuration
				******************************************************/
				// if ((muls.tds) || (muls.nCellZ % muls.cellDiv != 0)) {
				if (!muls.equalDivs) {
					make3DSlices(&muls,muls.slices,muls.atomPosFile,NULL);
					initSTEMSlices(&muls,muls.slices);
				}



				timer = cputim();
				runMulsSTEM(&muls,wave);
				muls.totalSliceCount += muls.slices;

				if (muls.printLevel > 0) {
					printf("t=%gA, int.=%g time: %gsec (avgCount=%d)\n",
						wave->thickness,wave->intIntensity,cputim()-timer,muls.avgCount);
				}

				/***************** FOR DEBUGGING ****************/
				if ((muls.avgCount == 0) && (muls.saveLevel >=0) && (pCount+1==muls.mulsRepeat2*muls.cellDiv)) {
					if (muls.tds) comment = "Test wave function for run 0";
					else comment = "Exit face wave function for no TDS";
					sprintf(systStr,"%s/wave.img",muls.folder);
					if ((muls.tiltBack) && ((muls.btiltx != 0) || (muls.btilty != 0))) {
						ktx = -2.0*pi*sin(muls.btiltx)/wavelength(muls.v0);
						kty = -2.0*pi*sin(muls.btilty)/wavelength(muls.v0);
						for (ix=0;ix<muls.nx;ix++) {
							x = muls.resolutionX*(ix-muls.nx/2);
							for (iy=0;iy<muls.ny;iy++) {
								y = muls.resolutionY*(ix-muls.nx/2);
								wave->wave[ix][iy][0] *= cos(ktx*x+kty*y);
								wave->wave[ix][iy][1] *= sin(ktx*x+kty*y);
							}
						}
						if (muls.printLevel > 1) printf("** Applied beam tilt compensation **\n");
					}

					wave->WriteWave(systStr, comment);
				}
#ifdef VIB_IMAGE_TEST  // doTEM
				if ((muls.tds) && (muls.saveLevel > 2)) {
					sprintf(systStr,"%s/wave_%d.img",muls.folder,muls.avgCount);
					params = std::vector<double>(9);
					params[0] = muls.v0;  				// high voltage
					params[1] = muls.Cs;				// spherical aberration
					params[2] = muls.df0;				// defocus
					params[3] = muls.astigMag;			// astigmatism
					params[4] = muls.astigAngle;
					params[5] = muls.Cc * sqrt(muls.dE_E*muls.dE_E+muls.dV_V*muls.dV_V+muls.dI_I*muls.dI_I);	// focal spread
					// printf("****  Cc = %f, dE_E = %f, Delta = %f ****\n",muls.Cc,muls.dV_V,muls.Cc * muls.dV_V);
					params[6] = muls.alpha;				// illumination convergence angle
					// beam tilt:
					params[7] = muls.btiltx;			// beam tilt in mrad
					params[8] = muls.btilty;			// beam tilt in mrad

					comment = "complex exit face Wave function";

					wave->WriteWave(systStr, comment, params);
				}
#endif

			}
			result = readparam("sequence: ",buf,0);
		}
		/////////////////////////////////////////////////////////////////////////////
		// finished propagating through whole sample, we're at the exit surface now.
		// This means the wave function is used for nothing else than producing image(s)
		// and diffraction patterns.
		//////////////////////////////////////////////////////////////////////////////

		sprintf(avgName,"%s/diff.img",muls.folder);

		wave->ReadDiffPat(avgName);

		if (muls.avgCount == 0) {
			/***********************************************************
			* Save the diffraction pattern
			**********************************************************/
			memcpy((void *)wave->avgArray[0],(void *)wave->diffpat[0],(size_t)(muls.nx*muls.ny*sizeof(real)));
			/* move the averaged (raw data) file to the target directory as well */
#ifndef WIN32
			sprintf(avgName,"diffAvg_%d.img",muls.avgCount+1);
			sprintf(systStr,"mv %s/diff.img %s/%s",muls.folder,muls.folder,avgName);
			system(systStr);
#else
			sprintf(avgName,"diffAvg_%d.img",muls.avgCount+1);
			sprintf(systStr,"move %s\\diff.img %s/%s",muls.folder,muls.folder,avgName);
			// system(systStr);
#endif
			/***********************************************************
			* Save the Pendelloesung Plot
			**********************************************************/
			if (muls.lbeams) {
				for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
					for (ix=0;ix<muls.nbout;ix++) {
						avgPendelloesung[ix][iy] = muls.pendelloesung[ix][iy];
					}
				}
			}
			/***********************************************************
			* Save the defocused image, we can do with the wave what
			* we want, since it is not used after this anymore.
			* We will therefore multiply with the transfer function for
			* all the different defoci, inverse FFT and save each image.
			* diffArray will be overwritten with the image.
			**********************************************************/
			if (imageWave == NULL) imageWave = complex2Df(muls.nx,muls.ny,"imageWave");
			// multiply wave (in rec. space) with transfer function and write result to imagewave
			fftwf_execute(wave->fftPlanWaveForw);
			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				// here, we apply the CTF:
				imageWave[ix][iy][0] = wave->wave[ix][iy][0];
				imageWave[ix][iy][1] = wave->wave[ix][iy][1];
			}
			fftwf_execute_dft(wave->fftPlanWaveInv,imageWave[0],imageWave[0]);
			// get the amplitude squared:
			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				wave->diffpat[ix][iy] = imageWave[ix][iy][0]*imageWave[ix][iy][0]+imageWave[ix][iy][1]*imageWave[ix][iy][1];
			}
			sprintf(avgName,"%s/waveIntensity.img",muls.folder);
			wave->WriteDiffPat(avgName, "Wave intensity");
			// End of Image writing (if avgCount = 0)
			//////////////////////////////////////////////////////////////////////

		} // of if muls.avgCount == 0 ...
		else {
			/* 	 readRealImage_old(avgArray,muls.nx,muls.ny,&t,"diffAvg.img"); */
			muls.chisq[muls.avgCount-1] = 0.0;
			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				t = ((real)muls.avgCount*wave->avgArray[ix][iy]+
					wave->diffpat[ix][iy])/((real)(muls.avgCount+1));
				muls.chisq[muls.avgCount-1] += (wave->avgArray[ix][iy]-t)*(wave->avgArray[ix][iy]-t);
				wave->avgArray[ix][iy] = t;
			}
			muls.chisq[muls.avgCount-1] = muls.chisq[muls.avgCount-1]/(double)(muls.nx*muls.ny);
			sprintf(avgName,"%s/diffAvg_%d.img",muls.folder,muls.avgCount+1);
			wave->WriteAvgArray(avgName, "Diffraction pattern");

			/* write the data to a file */
			if ((avgFp = fopen("avgresults.dat","w")) == NULL )
				printf("Sorry, could not open data file for averaging\n");
			else {
				for (ix =0;ix<muls.avgCount;ix++) {
					fprintf(avgFp,"%d %g\n",ix+1,muls.chisq[ix]);
				}
				fclose(avgFp);
			}
			/*************************************************************/

			/***********************************************************
			* Average over the pendelloesung plot as well
			*/
			if (muls.lbeams) {
				for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
					for (ix=0;ix<muls.nbout;ix++) {
						avgPendelloesung[ix][iy] =
							((real)muls.avgCount*avgPendelloesung[ix][iy]+
							muls.pendelloesung[ix][iy])/(real)(muls.avgCount+1);
					}
				}
			}
			/***********************************************************
			* Save the defocused image, we can do with the wave what
			* we want, since it is not used after this anymore.
			* We will therefore multiply with the transfer function for
			* all the different defoci, inverse FFT and save each image.
			* diffArray will be overwritten with the image.
			**********************************************************/
			if (imageWave == NULL) imageWave = complex2Df(muls.nx,muls.ny,"imageWave");
			// multiply wave (in rec. space) with transfer function and write result to imagewave
#if FLOAT_PRECISION == 1
			fftwf_execute(wave->fftPlanWaveForw);
#elif FLOAT_PRECISION == 2
			fftw_execute(wave->fftPlanWaveForw);
#endif

			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				imageWave[ix][iy][0] = wave->wave[ix][iy][0];
				imageWave[ix][iy][1] = wave->wave[ix][iy][1];
			}
#if FLOAT_PRECISION == 1
			fftwf_execute_dft(wave->fftPlanWaveInv,imageWave[0],imageWave[0]);
#elif FLOAT_PRECISION == 2
			fftw_execute_dft(wave->fftPlanWaveInv,imageWave[0],imageWave[0]);
#endif

			// save the amplitude squared:
			sprintf(avgName,"%s/image.img",muls.folder);
			wave->ReadDiffPat(avgName);
			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				t = ((real)muls.avgCount*wave->diffpat[ix][iy]+
					imageWave[ix][iy][0]*imageWave[ix][iy][0]+imageWave[ix][iy][1]*imageWave[ix][iy][1])/(real)(muls.avgCount+1);
				wave->diffpat[ix][iy] = t;
			}
			wave->WriteDiffPat(avgName, "Image intensity");
			// End of Image writing (if avgCount > 0)
			//////////////////////////////////////////////////////////////////////

		} /* else ... if avgCount was greater than 0 */


		/////////////////////////////////////////////////////
		// Save the Pendelloesung plot:
		if (muls.lbeams) {
			/**************************************************************
			* The diffraction spot intensities of the selected
			* diffraction spots are now stored in the 2 dimensional array
			* muls.pendelloesung[beam][slice].
			* We can write the array to a file and display it, just for
			* demonstration purposes
			*************************************************************/
			sprintf(avgName,"%s/pendelloesung.dat",muls.folder);
			if ( (fpTEM = fopen( avgName, "w" )) != NULL ) {
				printf("Writing Pendelloesung data\n");
				for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
					/* write the thicknes in the first column of the file */
					fprintf( fpTEM, "%g", iy*muls.c / ((float)(muls.slices*muls.cellDiv)) );
					/* write the beam intensities in the following columns */
					for (ix=0;ix<muls.nbout;ix++) {
						// store the AMPLITUDE:
						fprintf( fpTEM, "\t%g", sqrt( avgPendelloesung[ix][iy] / (muls.nx*muls.ny) ) );
					}
					/* close the line, and start a new one for the next set of
					* intensities
					*/
					fprintf( fpTEM, "\n" );
				}
				fclose( fpTEM );
			}
			else {
				printf("Could not open file for pendelloesung plot\n");
			}
		} /* end of if lbemas ... */
		displayProgress(1);
	} /* end of for muls.avgCount=0.. */
}
/************************************************************************
* end of doTEM
***********************************************************************/



/************************************************************************
* doSTEM performs a STEM calculation
*
***********************************************************************/

void doSTEM() {
	int ix=0,iy=0,i,pCount,picts,ixa,iya,totalRuns;
	double timer, total_time=0;
	char buf[BUF_LEN];
	real t;
	static real **avgArray=NULL;
	double collectedIntensity;

	std::vector<WavePtr> waves;
	WavePtr wave;

	//pre-allocate several waves (enough for one row of the scan.
	for (int th=0; th<omp_get_max_threads(); th++)
	{
		waves.push_back(WavePtr(new WAVEFUNC(muls.nx, muls.ny, muls.resolutionX, muls.resolutionY)));
	}

	muls.chisq = std::vector<double>(muls.avgRuns);
	totalRuns = muls.avgRuns;
	timer = cputim();

	/* average over several runs of for TDS */
	displayProgress(-1);

	for (muls.avgCount = 0;muls.avgCount < totalRuns; muls.avgCount++) {
		total_time = 0;
		collectedIntensity = 0;
		muls.totalSliceCount = 0;
		muls.dE_E = muls.dE_EArray[muls.avgCount];


		/****************************************
		* do the (big) loop
		*****************************************/
		pCount = 0;
		/* make sure we start at the beginning of the file
		so we won't miss any line that contains a sequence,
		because we will not do any EOF wrapping
		*/
		resetParamFile();
		while (readparam("sequence: ",buf,0)) {
			if (((buf[0] < 'a') || (buf[0] > 'z')) &&
				((buf[0] < '1') || (buf[0] > '9')) &&
				((buf[0] < 'A') || (buf[0] > 'Z'))) {
					printf("Stacking sequence: %s\n",buf);
					printf("Can only work with old stacking sequence\n");
					break;
			}

			// printf("Stacking sequence: %s\n",buf);

			picts = 0;
			/* for the dislocation models picts will be 1, because the atomcoordinates
			* are expressed explicitly for every atom in the whole specimen
			* For perfect Si samples mulsRepeat1 will remain 1, but picts will give
			* the number of unit cells in Z-direction, where the unit cell is defined by
			* the unit cell in the cssr file multiplied by NCELLZ.
			* cellDiv will usually be 1 in that case.
			*/
			sscanf(buf,"%d %d",&muls.mulsRepeat1,&picts);
			for (i=0;i<(int)strlen(buf);i++) buf[i] = 0;
			if (picts < 1) picts = 1;
			muls.mulsRepeat2 = picts;
			sprintf(muls.cin2,"%d",muls.mulsRepeat1);
			/* if the unit cell is divided into slabs, we need to multiply
			* picts by that number
			*/
			if ((picts > 1)&& (muls.cubex >0) && (muls.cubey >0) && (muls.cubez>0)) {
				printf("Warning: cube size of height %gA has been defined, ignoring sequence\n",muls.cubez);
				picts = 1;
			}
			picts *= muls.cellDiv;

			if (muls.equalDivs) {
				make3DSlices(&muls, muls.slices, muls.atomPosFile, NULL);
				initSTEMSlices(&muls, muls.slices);
				timer = cputim();
			}

			/****************************************
			* do the (small) loop over slabs
			*****************************************/
			for (pCount=0;pCount<picts;pCount++) {
				/*******************************************************
				* build the potential slices from atomic configuration
				******************************************************/
				if (!muls.equalDivs) {
					make3DSlices(&muls,muls.slices,muls.atomPosFile,NULL);
					initSTEMSlices(&muls,muls.slices);
					timer = cputim();
				}

				muls.complete_pixels=0;
				/**************************************************
				* scan through the different probe positions
				*************************************************/
				// default(none) forces us to specify all of the variables that are used in the parallel section.
				//    Otherwise, they are implicitly shared (and this was cause of several bugs.)
#pragma omp parallel \
	private(ix, iy, ixa, iya, wave, t, timer) \
	shared(pCount, picts, muls, collectedIntensity, total_time, waves) \
	default(none)
#pragma omp for
				for (i=0; i < (muls.scanXN * muls.scanYN); i++)
				{
					timer=cputim();
					ix = i / muls.scanYN;
					iy = i % muls.scanYN;

					wave = waves[omp_get_thread_num()];

					//printf("Scanning: %d %d %d %d\n",ix,iy,pCount,muls.nx);

					/* if this is run=0, create the inc. probe wave function */
					if (pCount == 0)
					{
						probe(&muls, wave, muls.nx/2*muls.resolutionX, muls.ny/2*muls.resolutionY);

						// TODO: modifying shared value from multiple threads?
						//muls.nslic0 = 0;
						//wave->thickness = 0.0;
					}

					else
					{
						/* load incident wave function and then propagate it */
						sprintf(wave->fileStart, "%s/mulswav_%d_%d.img", muls.folder, ix, iy);
						readStartWave(wave);  /* this also sets the thickness!!! */
						// TODO: modifying shared value from multiple threads?
						//muls.nslic0 = pCount;
					}
					/* run multislice algorithm
					   and save exit wave function for this position
					   (done by runMulsSTEM),
					   but we need to define the file name */
					sprintf(wave->fileout,"%s/mulswav_%d_%d.img",muls.folder,ix,iy);
					muls.saveFlag = 1;

					wave->iPosX =(int)(ix*(muls.scanXStop-muls.scanXStart)/
									  ((float)muls.scanXN*muls.resolutionX));
					wave->iPosY = (int)(iy*(muls.scanYStop-muls.scanYStart)/
									   ((float)muls.scanYN*muls.resolutionY));
					if (wave->iPosX > muls.potNx-muls.nx)
					{
						wave->iPosX = muls.potNx-muls.nx;
					}
					if (wave->iPosY > muls.potNy-muls.ny)
					{
						wave->iPosY = muls.potNy-muls.ny;
					}

					// MCS - update the probe wavefunction with its position
					wave->detPosX=ix;
					wave->detPosY=iy;

					runMulsSTEM(&muls,wave);


					/***************************************************************
					* In order to save some disk space we will add the diffraction
					* patterns to their averages now.  The diffraction pattern
					* should be stored in wave->diffpat (which each thread has independently),
					* if collectIntensity() has been executed correctly.
					***************************************************************/

					#pragma omp atomic
					collectedIntensity += wave->intIntensity;

					if (pCount == picts-1)  /* if this is the last slice ... */
					{
						sprintf(wave->avgName,"%s/diffAvg_%d_%d.img",muls.folder,ix,iy);
						// printf("Will copy to avgArray %d %d (%d, %d)\n",muls.nx, muls.ny,(int)(muls.diffpat),(int)avgArray);

						if (muls.saveLevel > 0)
						{
							if (muls.avgCount == 0)
							{
								// initialize the avgArray from the diffpat
								for (ixa=0;ixa<muls.nx;ixa++)
								{
									for (iya=0;iya<muls.ny;iya++)
									{
										wave->avgArray[ixa][iya]=wave->diffpat[ixa][iya];
									}
								}
							}
							else
							{
								// printf("Will read image %d %d\n",muls.nx, muls.ny);
								wave->ReadAvgArray(wave->avgName);
								for (ixa=0;ixa<muls.nx;ixa++) for (iya=0;iya<muls.ny;iya++) {
									t = ((real)muls.avgCount * wave->avgArray[ixa][iya] +
										wave->diffpat[ixa][iya]) / ((real)(muls.avgCount + 1));
									if (muls.avgCount>1)
									{
										#pragma omp atomic
										muls.chisq[muls.avgCount-1] += (wave->avgArray[ixa][iya]-t)*
											(wave->avgArray[ixa][iya]-t);
									}
									wave->avgArray[ixa][iya] = t;
								}
							}
							// Write the array to a file, resize and crop it,
							wave->WriteAvgArray(wave->avgName);
							}
							else {
								if (muls.avgCount > 0)	muls.chisq[muls.avgCount-1] = 0.0;
							}
					} /* end of if pCount == picts, i.e. conditional code, if this
						  * was the last slice
						  */

					#pragma omp atomic
					++muls.complete_pixels;

					if (muls.displayProgInterval > 0) if ((muls.complete_pixels) % muls.displayProgInterval == 0)
					{
						#pragma omp atomic
						total_time += cputim()-timer;
						printf("Pixels complete: (%d/%d), int.=%.3f, avg time per pixel: %.2fsec\n",
							muls.complete_pixels, muls.scanXN*muls.scanYN, wave->intIntensity,
							(total_time)/muls.complete_pixels);
						timer=cputim();
					}
				} /* end of looping through STEM image pixels */
				/* save STEM images in img files */
				saveSTEMImages(&muls);
				muls.totalSliceCount += muls.slices;
			} /* end of loop through thickness (pCount) */
		} /* end of  while (readparam("sequence: ",buf,0)) */
		// printf("Total CPU time = %f sec.\n", cputim()-timerTot );

		/*************************************************************/
		if (muls.avgCount>1)
			muls.chisq[muls.avgCount-1] = muls.chisq[muls.avgCount-1]/(double)(muls.nx*muls.ny);
		muls.intIntensity = collectedIntensity/(muls.scanXN*muls.scanYN);
		displayProgress(1);
	} /* end of loop over muls.avgCount */

}
