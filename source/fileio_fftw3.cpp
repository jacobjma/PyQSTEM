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

/* file: fileio_fftw3.c */

#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

// #include "../lib/floatdef.h"
#include "stemtypes_fftw3.h"
#include "memory_fftw3.h"	/* memory allocation routines */
// #include "tiffsubs.h"
#include "matrixlib.h"
#include "readparams.h"
#include "fileio_fftw3.h"
// #include "stemlib.h"

#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

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

#define N_A 6.022e+23
#define K_B 1.38062e-23      /* Boltzman constant */
#define PID 3.14159265358979 /* pi */
#define SQRT_PI 1.77245385090552 /* sqrt(pi) */
#define HPLANCK 6.6262e-34   /* planck's constant */
#define AMU 1.66053e-27      /* atomic mass unit */
#define THZ_AMU_PI_HBAR 0.05012012918415 /*  A°^2*THz*amu/(pi*hbar) */
#define THZ_AMU_HBAR 0.15745702964189    /*   A°^2*THz*amu/(hbar)   */
// 4.46677327584453 /* 1e10/sqrt(THz*amu/(pi*hbar)) */ 
#define THZ_HBAR_2KB  3.81927135604119     /* THz*hbar/(2*kB) */
#define THZ_HBAR_KB   1.90963567802059     /* THz*hbar/kB */
#define AMU_THZ2_A2_KB   1.20274224623720     /* AMU*THz^2*A^2/kB */

char *elTable = {
	"H HeLiBeB C N O F NeNaMgAlSiP S Cl"
	"ArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBr"
	"KrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTe"
	"I XeCsBaLaCePrNdPmSmEuGdTbDyHoErTm"
	"YbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn"
	"FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLr"
};

#define MAX_MASS_INDEX 59
int idArraySize=0;
int *idArray = NULL;
int idArrayPtr = 0;
double massArray[MAX_MASS_INDEX]={1.008,                                         4.0026,
6.939,9.012,      10.81,12.01,14.01,16.00,19.00,20.18,
23.00,24.31,      26.98,28.09,30.97,32.06,35.45,39.95,
39.10,40.08,
44.96,47.90,50.94,52.00,54.94,55.85,58.93,58.71,63.55,65.38,
69.72,72.59,74.92,78.96,79.90,83.80,
85.47,87.62,88.91,91.22,92.91,95.94,98.91,10.07,102.9,106.4,107.9,112.4,
114.8,118.7,121.8,127.6,126.9,131.3,
132.9054,137.33,138.9055,178.49,180.9479};
/* so far this list goes up to Xe (from Gerthsen) .. Ta (webelememts) */
double chargeTable[MAX_MASS_INDEX];


int atomCompareZnum(const void *atPtr1,const void *atPtr2) {
	int comp = 0;
	atom *atom1, *atom2;

	atom1 = (atom *)atPtr1;
	atom2 = (atom *)atPtr2;
	/* Use the fact that z is the first element in the atom struct */
	comp = (atom1->Znum == atom2->Znum) ? 0 : ((atom1->Znum > atom2->Znum) ? -1 : 1); 
	return comp;
}

int atomCompareZYX(const void *atPtr1,const void *atPtr2) {
	int comp = 0;
	atom *atom1, *atom2;

	atom1 = (atom *)atPtr1;
	atom2 = (atom *)atPtr2;
	/* Use the fact that z is the first element in the atom struct */
	comp = (atom1->z == atom2->z) ? 0 : ((atom1->z > atom2->z) ? 1 : -1); 
	if (comp == 0) {
		comp = (atom1->y == atom2->y) ? 0 : ((atom1->y > atom2->y) ? 1 : -1); 
		if (comp == 0) {
			comp = (atom1->x == atom2->x) ? 0 : ((atom1->x > atom2->x) ? 1 : -1); 
		}
	}
	return comp;
}


/********************************************************
* writePDB(atoms,natoms,fileName)
* This function will write the atomic coordinates in the 
* atoms array to the file <fileName> in pdb format, which
* is readable by AtomicEye
* The return value is 1 for success, or 0 for failure.
********************************************************/

int writePDB(atom *atoms,int natoms,char *fileName,MULS *muls) {
	FILE *writefp;
	int j,i;
	static char *elTable = {
		"H HeLiBeB C N O F NeNaMgAlSiP S Cl"
		"ArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBr"
		"KrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTe"
		"I XeCsBaLaCePrNdPmSmEuGdTbDyHoErTm"
		"YbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn"
		"FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLr"};
		char elem[16];
		double ax,by,cz;

		if (natoms < 1) {
			printf("Atom array empty - no file written\n");
			return 1;
		}

		printf( "DEBUG: writePDB: filename is %s \n", fileName );
		writefp = fopen(fileName, "w" );
		if( writefp == NULL ) {
			printf("Cannot open file %s\n",fileName);
			return 0;
		}

		ax = muls->ax;
		by=muls->by;
		cz=muls->c;


		fprintf( writefp,"HEADER    libAtoms:Config_save_as_pdb; %d atoms\n",natoms);
		fprintf( writefp,"CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1\n",
			muls->ax,muls->by,muls->c);
		elem[2] = '\0';
		for (j=0;j<natoms;j++) {
			elem[0] = elTable[2*atoms[j].Znum-2];
			elem[1] = elTable[2*atoms[j].Znum-1];
			if (elem[1] == ' ') elem[1] = '\0';
			fprintf( writefp,"ATOM   %4d %s",j+1,elem);
			for (i=strlen(elem);i<13;i++)
				fprintf( writefp," ");
			fprintf( writefp,"1   ");
			fprintf( writefp," %8.3f%8.3f%8.3f\n",atoms[j].x,atoms[j].y,atoms[j].z);
			/*
			fprintf( writefp," %8.3f%8.3f%8.3f\n",atoms[j].x/muls->ax,
			atoms[j].y/muls->by,atoms[j].z/muls->c);
			*/
		} 

		fclose( writefp );
		return 1;

}


////////////////////////////////////////////////////
// This function siomply removes those atoms from the list
// whose x- y- or z-position is exactly 1.0:
int removeRedundantAtoms(atom *atoms,int natoms) {

	return 0;
}


#define MIN_EDGE_LENGTH 5.18 /* minimal allowed edge length in A
* going below this limit will crash 
* AtomEye.
*/

int writeCFG(atom *atoms,int natoms,char *fileName,MULS *muls) {
	FILE *cfgfp;
	int j;
	/*
	static char *elTable = {
	"H HeLiBeB C N O F NeNaMgAlSiP S Cl"
	"ArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBr"
	"KrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTe"
	"I XeCsBaLaCePrNdPmSmEuGdTbDyHoErTm"
	"YbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn"
	"FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLr"};
	*/
	char elem[16];
	double ax,by,cz;

	if (natoms < 1) {
		printf("Atom array empty - no file written\n");
		return 1;
	}

	cfgfp = fopen(fileName, "w");
	if (cfgfp == NULL) {
		printf("Cannot open file %s\n",fileName);
		return 0;
	}

	printf( "DEBUG: fileio_fftw2::writeCFG writing to file = %s \n", fileName );

	// ax = muls->ax > MIN_EDGE_LENGTH ? muls->ax : MIN_EDGE_LENGTH;
	// by = muls->by > MIN_EDGE_LENGTH ? muls->by : MIN_EDGE_LENGTH;
	// cz = muls->c  > MIN_EDGE_LENGTH ? muls->c  : MIN_EDGE_LENGTH;
	ax = muls->ax;
	by = muls->by;
	cz = muls->c;

	fprintf( cfgfp, "Number of particles = %d\n", natoms);
	fprintf( cfgfp, "A = 1.0 Angstrom (basic length-scale)\n");
	fprintf( cfgfp, "H0(1,1) = %g A\nH0(1,2) = 0 A\nH0(1,3) = 0 A\n", ax);
	fprintf( cfgfp, "H0(2,1) = 0 A\nH0(2,2) = %g A\nH0(2,3) = 0 A\n", by);
	fprintf( cfgfp, "H0(3,1) = 0 A\nH0(3,2) = 0 A\nH0(3,3) = %g A\n", cz);
	fprintf( cfgfp, ".NO_VELOCITY.\nentry_count = 6\n");
	printf("ax: %g, by: %g, cz: %g n: %d\n",muls->ax,muls->by,muls->c,natoms);


	elem[2] = '\0';
	elem[0] = elTable[2*atoms[0].Znum-2];
	elem[1] = elTable[2*atoms[0].Znum-1];
	// printf("ax: %g, by: %g, cz: %g n: %d\n",muls->ax,muls->by,muls->c,natoms);
	if (elem[1] == ' ') elem[1] = '\0';
	fprintf( cfgfp, "%g\n%s\n", 2.0*atoms[0].Znum, elem);
	fprintf( cfgfp, "%g %g %g %g %g %g\n", atoms[0].x / ax, atoms[0].y / by, atoms[0].z / cz,
		atoms[0].dw,atoms[0].occ,atoms[0].q);



	for (j=1;j<natoms;j++) {
		if (atoms[j].Znum != atoms[j-1].Znum) {
			elem[0] = elTable[2*atoms[j].Znum-2];
			elem[1] = elTable[2*atoms[j].Znum-1];
			if (elem[1] == ' ') elem[1] = '\0';
			fprintf(cfgfp, "%g\n%s\n", 2.0*atoms[j].Znum, elem);
			// printf("%d: %g\n%s\n",j,2.0*atoms[j].Znum,elem);
		}
		fprintf(cfgfp, "%g %g %g %g %g %g\n", atoms[j].x / ax, atoms[j].y / by, atoms[j].z / cz,
			atoms[j].dw,atoms[j].occ,atoms[j].q);
		// if (atoms[j].occ != 1) printf("Atom %d: occ = %g\n",j,atoms[j].occ);
	} 
	fclose(cfgfp);

	return 1;
}


// write CFG file using atomic positions stored in pos, Z's in Znum and DW-factors in dw
// the unit cell is assumed to be cubic
int writeCFGFractCubic(double *pos,int *Znum,double *dw,int natoms,char *fileName,
					   double a,double b,double c) {
						   FILE *fracfp;
						   int j;
						   /*
						   static char *elTable = {
						   "H HeLiBeB C N O F NeNaMgAlSiP S Cl"
						   "ArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBr"
						   "KrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTe"
						   "I XeCsBaLaCePrNdPmSmEuGdTbDyHoErTm"
						   "YbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn"
						   "FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLr"};
						   */
						   char elem[16];
						   double ax,by,cz;

						   if (natoms < 1) {
							   printf("Atom array empty - no file written\n");
							   return 1;
						   }

						   fracfp = fopen(fileName, "w");
						   if (fracfp == NULL) {
							   printf("Cannot open file %s\n",fileName);
							   return 0;
						   }

						   ax = a > MIN_EDGE_LENGTH ? a : MIN_EDGE_LENGTH;
						   by = b > MIN_EDGE_LENGTH ? b : MIN_EDGE_LENGTH;
						   cz = c > MIN_EDGE_LENGTH ? c : MIN_EDGE_LENGTH;

						   fprintf(fracfp, "Number of particles = %d\n", natoms);
						   fprintf(fracfp, "A = 1.0 Angstrom (basic length-scale)\n");
						   fprintf(fracfp, "H0(1,1) = %g A\nH0(1,2) = 0 A\nH0(1,3) = 0 A\n", ax);
						   fprintf(fracfp, "H0(2,1) = 0 A\nH0(2,2) = %g A\nH0(2,3) = 0 A\n", by);
						   fprintf(fracfp, "H0(3,1) = 0 A\nH0(3,2) = 0 A\nH0(3,3) = %g A\n", cz);
						   fprintf(fracfp, ".NO_VELOCITY.\nentry_count = 4\n");
						   printf("ax: %g, by: %g, cz: %g n: %d\n",ax,by,c,natoms);


						   elem[2] = '\0';
						   elem[0] = elTable[2*Znum[0]-2];
						   elem[1] = elTable[2*Znum[0]-1];
						   // printf("ax: %g, by: %g, cz: %g n: %d\n",muls->ax,muls->by,muls->c,natoms);
						   if (elem[1] == ' ') elem[1] = '\0';
						   fprintf(fracfp, "%g\n%s\n", 2.0*Znum[0], elem);
						   fprintf(fracfp, "%g %g %g %g\n", pos[0] * a / ax, pos[1] * b / by, pos[2] * c / cz, dw[0]);



						   for (j=1;j<natoms;j++) {
							   if (Znum[j] != Znum[j-1]) {
								   elem[0] = elTable[2*Znum[j]-2];
								   elem[1] = elTable[2*Znum[j]-1];
								   if (elem[1] == ' ') elem[1] = '\0';
								   fprintf(fracfp, "%g\n%s\n", 2.0*Znum[j], elem);
								   // printf("%d: %g\n%s\n",j,2.0*atoms[j].Znum,elem);
							   }
							   fprintf(fracfp, "%g %g %g %g\n", pos[3 * j + 0] * a / ax, pos[3 * j + 1] * b / by, pos[3 * j + 2] * c / cz, dw[j]);
						   } 
						   fclose(fracfp);

						   return 1;
}




/*******************************************************************************
* int phononDisplacement: 
* This function will calculate the phonon displacement for a given atom i of the
* unit cell, which has been replicated to the larger cell (icx,icy,icz)
* The phonon displacement is either defined by the phonon data file, or, 
* the Einstein model, if the appropriate flags in the muls struct are set
* The displacement will be given in fractional coordinates of a single unit cell.
*
* Input parameters:
* Einstein-mode:
* need only: dw, Znum, atomCount
* atomCount: give statistics report, if 0, important only for non-Einstein mode
* maxAtom: total number of atoms (will be called first, i.e. atomCount=maxAtoms-1:-1:0)
*
* Phonon-file mode:
* ...
*
********************************************************************************/ 
//  phononDisplacement(u,muls,jChoice,icx,icy,icz,j,atoms[jChoice].dw,*natom,jz);
//  j == atomCount
int phononDisplacement(double *u,MULS *muls,int id,int icx,int icy,
	int icz,int atomCount,double dw,int maxAtom,int ZnumIndex) {
	int ix,iy,idd; // iz;
	static FILE *fpPhonon = NULL;
	static int Nk, Ns;        // number of k-vectors and atoms per primitive unit cell
	static float *massPrim;   // masses for every atom in primitive basis
	static float **omega;     // array of eigenvalues for every k-vector 
	static fftwf_complex ***eigVecs;  // array of eigenvectors for every k-vector
	static float **kVecs;     // array for Nk 3-dim k-vectors
	static double **q1=NULL, **q2=NULL;
	int ik,lambda,icoord; // Ncells, nkomega;
	double kR,kRi,kRr,wobble;
	static double *u2=NULL,*u2T,ux=0,uy=0,uz=0; // u2Collect=0; // Ttotal=0;
	// static double uxCollect=0,uyCollect=0,uzCollect=0;
	static int *u2Count = NULL,*u2CountT,runCount = 1,u2Size = -1;
	static long iseed=0;
	static double **Mm=NULL,**MmInv=NULL;
	// static double **MmOrig=NULL,**MmOrigInv=NULL;
	static double *axCell,*byCell,*czCell,*uf,*b;
	static double wobScale = 0,sq3,scale=0;

	if (muls->tds == 0) return 0;

	/***************************************************************************
	* We will give a statistics report, everytime, when we find atomCount == 0
	***************************************************************************/

	if (atomCount == 0) {
							   for (ix=0;ix<muls->atomKinds;ix++) {
								   // u2Collect += u2[ix]/u2Count[ix];
								   // uxCollect += ux/maxAtom; uyCollect += uy/maxAtom; uzCollect += uz/maxAtom;
								   /*
								   printf("STATISTICS: sqrt(<u^2>): %g, CM: (%g %g %g) %d atoms, wob: %g\n"
								   "                         %g, CM: (%g %g %g) run %d\n",
								   sqrt(u2/u2Count),ux/u2Count,uy/u2Count,uz/u2Count,u2Count,scale*sqrt(dw*wobScale),
								   sqrt(u2Collect/runCount),uxCollect/runCount,uyCollect/runCount,uzCollect/runCount,
								   runCount);
								   */
								   // printf("Count: %d %g\n",u2Count[ix],u2[ix]);
								   u2[ix] /= u2Count[ix];
								   if (runCount > 0) 
									   muls->u2avg[ix] = sqrt(((runCount-1)*(muls->u2avg[ix]*muls->u2avg[ix])+u2[ix])/runCount);
								   else
									   muls->u2avg[ix] = sqrt(u2[ix]);

								   muls->u2[ix]    = sqrt(u2[ix]);

								   u2[ix]=0; u2Count[ix]=0;
							   }
							   runCount++;
							   ux=0; uy=0; uz=0; 
	}
	if (Mm == NULL) {
									   // MmOrig = double2D(3,3,"MmOrig");
									   // MmOrigInv = double2D(3,3,"MmOrigInv");
									   Mm = double2D(3,3,"Mm");
									   MmInv = double2D(3,3,"Mminv");
									   uf = double1D(3,"uf");
		b = double1D(3,"uf");

									   axCell=Mm[0]; byCell=Mm[1]; czCell=Mm[2];
		// memcpy(Mm[0],muls->Mm[0],3*3*sizeof(double));
		// We need to copy the transpose of muls->Mm to Mm.
		// we therefore cannot use the following command:
		// memcpy(Mm[0],muls->Mm[0],3*3*sizeof(double));
		// or Mm = muls->Mm;
		for (ix=0;ix<3;ix++) for (iy=0;iy<3;iy++) Mm[ix][iy]=muls->Mm[iy][ix];

		
									   // makeCellVectMuls(muls, axCell, byCell, czCell);
									   // memcpy(MmOrig[0],Mm[0],3*3*sizeof(double));
									   inverse_3x3(MmInv[0],Mm[0]);
	}
	if (ZnumIndex >= u2Size) {
							   
	   // printf("created phonon displacements %d!\n",ZnumIndex);								
	   u2 = (double *)realloc(u2,(ZnumIndex+1)*sizeof(double));
	   u2Count = (int *)realloc(u2Count,(ZnumIndex+1)*sizeof(int));
							   
		// printf("%d .... %d\n",ZnumIndex,u2Size);
	   if (u2Size < 1) {
		   for (ix=0;ix<=ZnumIndex;ix++) {
			   u2[ix] = 0;
			   u2Count[ix] = 0;
		   }
	   }
	   else {
		   for (ix=u2Size;ix<=ZnumIndex;ix++) {
			   u2[ix] = 0;
			   u2Count[ix] = 0;
		   }
	   }						  
		// printf("%d ..... %d\n",ZnumIndex,u2Size);

	   u2Size = ZnumIndex+1;
	}  


						   /***************************************************************************
						   * Thermal Diffuse Scattering according to accurate phonon-dispersion or 
						   * just Einstein model
						   *
						   * Information in the phonon file will be stored in binary form as follows:
						   * Nk (number of k-points: 32-bit integer)
						   * Ns (number of atomic species 32-bit integer)
						   * M_1 M_2 ... M_Ns  (32-bit floats)
						   * kx(1) ky(1) kz(1) (32-bit floats)
						   * w_1(1) q_11 q_21 ... q_(3*Ns)1    (32-bit floats)
						   * w_2(1) q_12 q_22 ... q_(3*Ns)2
						   * :
						   * w_(3*Ns)(1) q_1Ns q_2Ns ... q_(3*Ns)Ns
						   * kx(2) ky(2) kz(2) 
						   * :
						   * kx(Nk) ky(Nk) kz(Nk)
						   * :
						   * 
						   *
						   * Note: only k-vectors in half of the Brioullin zone must be given, since 
						   * w(k) = w(-k)  
						   * also: 2D arrays will be read slowly varying index = first index (i*m+j)
						   **************************************************************************/

						   if (wobScale == 0) {
							   wobScale = 1.0/(8*PID*PID);   
		sq3 = 1.0/sqrt(3.0);  /* sq3 is an additional needed factor which stems from
												   * int_-infty^infty exp(-x^2/2) x^2 dx = sqrt(pi)
												   * introduced in order to match the wobble factor with <u^2>
												   */
							   scale = (float) sqrt(muls->tds_temp/300.0) ;
							   iseed = -(long)(time(NULL));
						   }


						   if ((muls->Einstein == 0) && (fpPhonon == NULL)) {
							   if ((fpPhonon = fopen(muls->phononFile,"r")) == NULL) {
								   printf("Cannot find phonon mode file, will use random displacements!");
								   muls->Einstein = 1;
								   //muls->phononFile = NULL;
							   }
							   else {

								   if (2*sizeof(float) != sizeof(fftwf_complex)) {
									   printf("phononDisplacement: data type mismatch: fftw_complex != 2*float!\n");
									   exit(0);
								   }
								   fread(&Nk,sizeof(int),1,fpPhonon);
								   fread(&Ns,sizeof(int),1,fpPhonon);
								   massPrim =(float *)malloc(Ns*sizeof(float));  // masses for every atom in primitive basis
								   fread(massPrim,sizeof(float),Ns,fpPhonon);
								   kVecs = float32_2D(Nk,3,"kVecs");
								   omega = float32_2D(Nk,3*Ns,"omega");          /* array of eigenvalues for every k-vector 
																				 * omega is given in THz, but the 2pi-factor
																				 * is still there, i.e. f=omega/2pi
																				 */
								   eigVecs = complex3Df(Nk,3*Ns,3*Ns,"eigVecs"); // array of eigenvectors for every k-vector
								   for (ix=0;ix<Nk;ix++) {
									   fread(kVecs[ix],sizeof(float),3,fpPhonon);  // k-vector
									   for (iy=0;iy<3*Ns;iy++) {
										   fread(omega[ix]+iy,sizeof(float),1,fpPhonon);
										   fread(eigVecs[ix][iy],2*sizeof(float),3*Ns,fpPhonon);
									   }	
								   }
								   /*
								   printf("Masses: ");
								   for (ix=0;ix<Ns;ix++) printf(" %g",massPrim[ix]);
								   printf("\n");
								   for (ix=0;ix<3;ix++) {
								   printf("(%5f %5f %5f):  ",kVecs[ix][0],kVecs[ix][1],kVecs[ix][2]);
								   for (idd=0;idd<Ns;idd++) for (iy=0;iy<3;iy++)
								   printf("%6g ",omega[ix][iy+3*idd]);
								   printf("\n");
								   }
								   */      
								   /* convert omega into q scaling factors, since we need those, instead of true omega:    */
								   /* The 1/sqrt(2) term is from the dimensionality ((q1,q2) -> d=2)of the random numbers */
								   for (ix=0;ix<Nk;ix++) {
									   for (idd=0;idd<Ns;idd++) for (iy=0;iy<3;iy++) {
										   // quantize the energy distribution:
										   // tanh and exp give different results will therefore use exp
										   // nkomega = (int)(1.0/tanh(THZ_HBAR_2KB*omega[ix][iy+3*id]/muls->tds_temp));
										   // wobble  =      (1.0/tanh(THZ_HBAR_2KB*omega[ix][iy+3*id]/muls->tds_temp)-0.5);
										   // nkomega = (int)(1.0/(exp(THZ_HBAR_KB*omega[ix][iy+3*id]/muls->tds_temp)-1)+0.5);
										   if (omega[ix][iy+3*idd] > 1e-4) {
											   wobble = muls->tds_temp>0 ? (1.0/(exp(THZ_HBAR_KB*omega[ix][iy+3*idd]/muls->tds_temp)-1)):0;
											   // if (ix == 0) printf("%g: %d %g\n",omega[ix][iy+3*id],nkomega,wobble);
											   wobble = sqrt((wobble+0.5)/(2*PID*Nk*2*massPrim[idd]*omega[ix][iy+3*idd]* THZ_AMU_HBAR));  
										   }
										   else wobble = 0;
										   /* Ttotal += 0.25*massPrim[id]*((wobble*wobble)/(2*Ns))*
										   omega[ix][iy+3*id]*omega[ix][iy+3*id]*AMU_THZ2_A2_KB;
										   */
										   omega[ix][iy+3*idd] = wobble;
									   }  // idd
									   // if (ix == 0) printf("\n");
								   }
								   // printf("Temperature: %g K\n",Ttotal);
								   // printf("%d %d %d\n",(int)(0.4*(double)Nk/11.0),(int)(0.6*(double)Nk),Nk);
								   q1 = double2D(3*Ns,Nk,"q1");
								   q2 = double2D(3*Ns,Nk,"q2");

							   }
							   fclose(fpPhonon);    
						   }  // end of if phononfile

						   // 
						   // in the previous bracket: the phonon file is only read once.
						   /////////////////////////////////////////////////////////////////////////////////////
						   if ((muls->Einstein == 0) && (atomCount == maxAtom-1)) {
							   if (Nk > 800)
								   printf("Will create phonon displacements for %d k-vectors - please wait ...\n",Nk);
							   for (lambda=0;lambda<3*Ns;lambda++) for (ik=0;ik<Nk;ik++) {
								   q1[lambda][ik] = (omega[ik][lambda] * gasdev( &iseed ));
								   q2[lambda][ik] = (omega[ik][lambda] * gasdev( &iseed ));
							   }
							   // printf("Q: %g %g %g\n",q1[0][0],q1[5][8],q1[0][3]);
						   }
   /********************************************************************************
	* Do the Einstein model independent vibrations !!!
	*******************************************************************************/
	if (muls->Einstein) {	    
	   /* convert the Debye-Waller factor to sqrt(<u^2>) */
	   wobble = scale*sqrt(dw*wobScale);
	   u[0] = (wobble*sq3 * gasdev( &iseed ));
	   u[1] = (wobble*sq3 * gasdev( &iseed ));
	   u[2] = (wobble*sq3 * gasdev( &iseed ));
	   ///////////////////////////////////////////////////////////////////////
	   // Book keeping:
		u2[ZnumIndex] += u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
		ux += u[0]; uy += u[1]; uz += u[2];
		u2Count[ZnumIndex]++;


	  /* Finally we must convert the displacement for this atom back into its fractional
	   * coordinates so that we can add it to the current position in vector a
	   */
	   matrixProduct(&u,1,3,MmInv,3,3,&uf);
// test:
	   /*
	   matrixProduct(&uf,1,3,Mm,3,3,&b);
	   if (atomCount % 5 == 0) {
			printf("Z = %d, DW = %g, u=[%g %g %g]\n       u'=[%g %g %g]\n",muls->Znums[ZnumIndex],dw,u[0],u[1],u[2],b[0],b[1],b[2]);
			showMatrix(Mm,3,3,"Mm");
			showMatrix(MmInv,3,3,"MmInv");
	   }
	   */
// end test
	   memcpy(u,uf,3*sizeof(double));
	}
	else {
	   // id seems to be the index of the correct atom, i.e. ranges from 0 .. Natom
	   printf("created phonon displacements %d, %d, %d %d (eigVecs: %d %d %d)!\n",ZnumIndex,Ns,Nk,id,Nk,3*Ns,3*Ns);
	   /* loop over k and lambda:  */
	   memset(u,0,3*sizeof(double));
	   for (lambda=0;lambda<3*Ns;lambda++) for (ik=0;ik<Nk;ik++) {
		   // if (kVecs[ik][2] == 0){
		   kR = 2*PID*(icx*kVecs[ik][0]+icy*kVecs[ik][1]+icz*kVecs[ik][2]);
		   //  kR = 2*PID*(blat[0][0]*kVecs[ik][0]+blat[0][1]*kVecs[ik][1]+blat[0][2]*kVecs[ik][2]);
		   kRr = cos(kR); kRi = sin(kR);
		   for (icoord=0;icoord<3;icoord++) {
									   u[icoord] += q1[lambda][ik]*(eigVecs[ik][lambda][icoord+3*id][0]*kRr-
										   eigVecs[ik][lambda][icoord+3*id][1]*kRi)-
										   q2[lambda][ik]*(eigVecs[ik][lambda][icoord+3*id][0]*kRi+
										   eigVecs[ik][lambda][icoord+3*id][1]*kRr);
		   }
		}
		// printf("u: %g %g %g\n",u[0],u[1],u[2]);
	   /* Convert the cartesian displacements back to reduced coordinates
	   */ 
	   ///////////////////////////////////////////////////////////////////////
	   // Book keeping:
		u2[ZnumIndex] += u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
		ux += u[0]; uy += u[1]; uz += u[2];
		u[0] /= muls->ax;
		u[1] /= muls->by;
		u[2] /= muls->c;
		u2Count[ZnumIndex]++;

   } /* end of if Einstein */

	// printf("%d ... %d\n",ZnumIndex,u2Size);
   // printf("atomCount: %d (%d) %d %d\n",atomCount,muls->Einstein,ZnumIndex,u2Size);


	/*
	// Used for Debugging the memory leak on Aug. 20, 2010:
	if (_CrtCheckMemory() == 0) {
	   printf("Found bad memory check in phononDisplacement! %d %d\n",ZnumIndex,muls->atomKinds);
	}
	*/


	return 0;  
}


/***********************************************************************
* The following function returns the number of atoms in the specified
* CFG file and updates the cell parameters in the muls struct
***********************************************************************/
int readDATCellParams(MULS *muls, double **Mm, char *fileName) {
	int ncoord;
	char buf[256];
	double a,b,c,alpha=90.0,beta=90.0,gamma=90.0;

	// printf("Paramter File pointer (1): %d\n",(int)getFp());
	parFpPush(); /* push the current parameter file file pointer 
				 * on the stack to make this function totaly 
				 * transparent 
				 */
	if (!parOpen(fileName)) {
		printf("Could not open CFG input file %s\n",fileName);
		parFpPull();  /* restore old parameter file pointer */
		return 0;
	}
	// printf("Paramter File pointer (2): %d\n",(int)getFp());

	resetParamFile();  
	setComment('#');  
	if (readparam("Number of atoms =",buf,1)) sscanf(buf,"%d",&ncoord);

	if (readparam("a =",buf,1)) sscanf(buf,"%lf",&a);
	if (readparam("b =",buf,1)) sscanf(buf,"%lf",&b);
	if (readparam("c =",buf,1)) sscanf(buf,"%lf",&c);

	if (readparam("alpha =",buf,1)) sscanf(buf,"%lf",&alpha);
	if (readparam("beta =",buf,1)) sscanf(buf,"%lf",&beta);
	if (readparam("gamma =",buf,1)) sscanf(buf,"%lf",&gamma);

	setComment('%');
	parClose();   
	parFpPull();  /* restore old parameter file pointer */


	muls->ax = a;
	muls->by = b;
	muls->c  = c;
	muls->cGamma = alpha;
	muls->cBeta  = beta;
	muls->cAlpha = gamma;
	// construct the unit cell metric from the lattice parameters and angles:
	makeCellVectMuls(muls, Mm[0], Mm[1], Mm[2]);   
	if (ncoord < 1) {
		printf("Number of atoms in CFG file not specified!\n");
		ncoord = 0;
	}
	return ncoord;
}



/***********************************************************************
* The following function returns the number of atoms in the specified
* CFG file and updates the cell parameters in the muls struct
***********************************************************************/
int readCFGCellParams(MULS *muls, double **Mm, char *fileName) {
	int ncoord,i;
	char buf[256];
	double lengthScale;

	// printf("Paramter File pointer (1): %d\n",(int)getFp());
	parFpPush(); /* push the current parameter file file pointer 
				 * on the stack to make this function totaly 
				 * transparent 
				 */
	if (!parOpen(fileName)) {
		printf("Could not open CFG input file %s\n",fileName);
		parFpPull();  /* restore old parameter file pointer */
		return 0;
	}
	// printf("Paramter File pointer (2): %d\n",(int)getFp());

	resetParamFile();  
	setComment('#');  
	if (readparam("Number of particles =",buf,1)) sscanf(buf,"%d",&ncoord);
	if (readparam("A =",buf,1)) sscanf(buf,"%lf",&lengthScale);

	if (readparam("H0(1,1) =",buf,1)) sscanf(buf,"%lf",Mm[0]+0);
	if (readparam("H0(1,2) =",buf,1)) sscanf(buf,"%lf",Mm[0]+1);
	if (readparam("H0(1,3) =",buf,1)) sscanf(buf,"%lf",Mm[0]+2);

	if (readparam("H0(2,1) =",buf,1)) sscanf(buf,"%lf",Mm[0]+3);
	if (readparam("H0(2,2) =",buf,1)) sscanf(buf,"%lf",Mm[0]+4);
	if (readparam("H0(2,3) =",buf,1)) sscanf(buf,"%lf",Mm[0]+5);

	if (readparam("H0(3,1) =",buf,1)) sscanf(buf,"%lf",Mm[0]+6);
	if (readparam("H0(3,2) =",buf,1)) sscanf(buf,"%lf",Mm[0]+7);
	if (readparam("H0(3,3) =",buf,1)) sscanf(buf,"%lf",Mm[0]+8);
	/*
	if (readparam(".NO_VELOCITY.",buf,1)) noVelocityFlag = 1; 
	if (readparam("entry_count =",buf,1)) sscanf(buf,"%lf",&entryCount);
	if (!noVelocityFlag) entryCount+=3;
	*/
	setComment('%');
	parClose();   
	parFpPull();  /* restore old parameter file pointer */

	for (i=0;i<9;i++) Mm[0][i] *= lengthScale;

	muls->ax = sqrt(Mm[0][0]*Mm[0][0]+Mm[0][1]*Mm[0][1]+Mm[0][2]*Mm[0][2]);
	muls->by = sqrt(Mm[1][0]*Mm[1][0]+Mm[1][1]*Mm[1][1]+Mm[1][2]*Mm[1][2]);
	muls->c  = sqrt(Mm[2][0]*Mm[2][0]+Mm[2][1]*Mm[2][1]+Mm[2][2]*Mm[2][2]);
	muls->cGamma = atan2(Mm[1][1],Mm[1][0]);
	muls->cBeta = acos(Mm[2][0]/muls->c);
	muls->cAlpha = acos(Mm[2][1]*sin(muls->cGamma)/muls->c+cos(muls->cBeta)*cos(muls->cGamma));
	muls->cGamma /= (float)PI180;
	muls->cBeta  /= (float)PI180;
	muls->cAlpha /= (float)PI180;
	if (ncoord < 1) {
		printf("Number of atoms in CFG file not specified!\n");
		ncoord = 0;
	}
	return ncoord;
}

/***********************************************************************
* The following function returns the number of atoms in the specified
* CFG file and updates the cell parameters in the muls struct
***********************************************************************/
int readCSSRCellParams(MULS *muls, double **Mm, char *fileName) {
	FILE *cellfp;
	char s1[32],s2[32],s3[32],buf[NCMAX];
	int spaceGrp = 1,ncoord;

	cellfp = fopen(fileName, "r");
	if (cellfp == NULL) {
		printf("Cannot open file %s\n",fileName);
		return 0;
	}
	/* read the cell parameters from first and secons line */
	ReadLine(cellfp, buf, NCMAX, "in ReadXYZcoord");
	sscanf( buf, " %s %s %s",s1,s2,s3);
	muls->ax = atof(s1); muls->by = atof(s2); muls->c = atof(s3);

	ReadLine(cellfp, buf, NCMAX, "in ReadXYZcoord");
	sscanf( buf, " %s %s %s",s1,s2,s3);
	muls->cAlpha = atof(s1); muls->cBeta  = atof(s2); muls->cGamma = atof(s3);
	makeCellVectMuls(muls, Mm[0], Mm[1], Mm[2]);   

	/* check the space group: */
	spaceGrp = atoi(strstr(buf,"SPGR =")+strlen("SPGR ="));
	if (spaceGrp != 1) {
		printf("cannot interpret space group %d\n",spaceGrp);
		exit(0);
	}
	ReadLine(cellfp, buf, NCMAX, "in ReadXYZcoord");
	ncoord = atoi(buf);
	fclose(cellfp);
	return ncoord;
}

/*******************************************************************************
* This function reads the atomic position and element data for a single atom
* from a .dat file.  The atomic positions are given in reduced coordinates.
* Calling this function several times will read one atom after another from the 
* file.
* This function will return -1, if the end of file is reached prematurely, 0 otherwise.
*******************************************************************************/

int readNextDATAtom(atom *newAtom, int flag, char *fileName) {
	int printFlag = 0;
	static FILE *atomfp = NULL;
	static int noVelocityFlag = 1,entryCount = 3,element = 1;
	static char buf[NCMAX];
	static double *atomData = NULL;
	char *str,elementStr[3];
	int j;

	if (flag < 0) {
		if (atomfp != NULL) {
			parClose();   
			parFpPull();  /* restore old parameter file pointer */      
			atomfp = NULL;
			setComment('%');
		}
		return -1;
	}

	if (atomfp == NULL) {
		parFpPush();  /* save old parameter file pointer */      
		if (!parOpen(fileName)) {
			printf("Could not open DAT input file %s\n",fileName);
			parFpPull();  /* restore old parameter file pointer */
			return -1;
		}
		resetParamFile();  
		// advance file pointer to last (known) header line
		readparam("gamma =",buf,1);
		atomfp = getFp();  /* get the file pointer from the parameter file routines */
		atomData = (double *)malloc(entryCount*sizeof(double));
	}

	element = 0;
	do {
		if (fgets(buf, NCMAX, atomfp) == NULL) {
			if (printFlag) printf("Found end of file!\n");
			return -1;
		}
		/* check, if this is a new element name */	
		// str = strnext(buf," \t");
		if (buf != NULL) {
			memcpy(elementStr,buf,2);
			element = getZNumber(elementStr);
		}
	} while (element == 0);
	if (printFlag) printf("Found Z=%d\n",element);
	if (element > 0) {
		str = buf+2;
		// skip leading spaces:
		while (strchr(" \t",*str) != NULL) str++; 
		for (j=0;j<entryCount;j++) {
			if (str==NULL) {
				printf("readNextCFGatom: Error: incomplete data line: >%s<\n",buf);
				return -1;
			}
			atomData[j] = atof(str); str=strnext(str," \t");
		}		
	}

	newAtom->Znum = element;
	newAtom->x    = atomData[0];
	newAtom->y    = atomData[1];
	newAtom->z    = atomData[2];
	newAtom->dw   = 0.45*28.0/(double)(2.0*element);	
	newAtom->occ  = 1.0;
	newAtom->q    = 0.0;
	printf("Atom: %d (%g %g %g), occ=%g, q=%g\n",newAtom->Znum,newAtom->x,newAtom->y,newAtom->z,newAtom->occ,newAtom->q);	  

	return 0;
}




/*******************************************************************************
* This function reads the atomic position and element data for a single atom
* from a .cfg file.  The atomic positions are given in reduced coordinates.
* Calling this function several times will read one atom after another from the 
* file.
* This function will return -1, if the end of file is reached prematurely, 0 otherwise.
*******************************************************************************/

int readNextCFGAtom(atom *newAtom, int flag, char *fileName) {
	static FILE *fpNextCFG=NULL;
	static int noVelocityFlag = 1,entryCount = 3,element = 1;
	static char buf[NCMAX];
	static double *atomData = NULL;
	static double mass = 28;
	char *str = NULL;
	int j;


	if (flag < 0) {
		if (fpNextCFG != NULL) {
			parClose();   
			parFpPull();  /* restore old parameter file pointer */
			fpNextCFG = NULL;
			setComment('%');
		}
		return -1;
	}

	if (fpNextCFG == NULL) {
		parFpPush();  /* save old parameter file pointer */      
		if (!parOpen(fileName)) {
			printf("Could not open CFG input file %s\n",fileName);
			parFpPull();  /* restore old parameter file pointer */
			return -1;
		}
		resetParamFile();  
		if (readparam(".NO_VELOCITY.",buf,1)) noVelocityFlag = 1; 
		else noVelocityFlag = 0;
		if (readparam("entry_count =",buf,1)) sscanf(buf,"%d",&entryCount);
		if (!noVelocityFlag) entryCount+=3;
		fpNextCFG = getFp();  /* get the file pointer from the parameter file routines */
		atomData = (double *)malloc((entryCount+1)*sizeof(double));
	}

	if (fgets(buf, NCMAX, fpNextCFG) == NULL) return -1;
	/* check, if this is a new mass number */
	str = strnext(buf," \t");
	if ((atof(buf) >= 1.0) && ((str==NULL) || (*str == '#'))) {
		mass = atof(buf);
		// printf("nV: %d, eC: %d (%g)\n",noVelocityFlag, entryCount,atof(buf));
		if (fgets(buf, NCMAX, fpNextCFG) == NULL) return -1;
		element = getZNumber(buf); 
		// printf("*** found element %d (%s %d) ***\n",element,buf,strlen(buf));
		if (fgets(buf, NCMAX, fpNextCFG) == NULL) return -1;
	}
	str = buf;
	// skip leading spaces:
	while (strchr(" \t",*str) != NULL) str++; 
	for (j=0;j<entryCount;j++) {
		if (str==NULL) {
			printf("readNextCFGatom: Error: incomplete data line: >%s<\n",buf);
			return -1;
		}
		atomData[j] = atof(str); str=strnext(str," \t");
	}


	newAtom->Znum = element;
	newAtom->x    = atomData[0];
	newAtom->y    = atomData[1];
	newAtom->z    = atomData[2];
	// newAtom->dw   = 0.45*28.0/((double)(2*element));	
	// printf("Element: %d, mass=%g\n",element,mass);
	newAtom->dw   = 0.45*28.0/mass;	
	newAtom->occ  = 1.0;
	newAtom->q    = 0.0;
	// read the DW-factor
	if (entryCount > 3+3*(1-noVelocityFlag)) 
		newAtom->dw = atomData[3+3*(1-noVelocityFlag)];
	// read the atom's occupancy:
	if (entryCount > 4+3*(1-noVelocityFlag)) 
		newAtom->occ = atomData[4+3*(1-noVelocityFlag)];
	// read the atom's charge:
	if (entryCount > 5+3*(1-noVelocityFlag)) 
		newAtom->q = atomData[5+3*(1-noVelocityFlag)];
	// printf("Atom: %d (%g %g %g), occ=%g, q=%g\n",newAtom->Znum,newAtom->x,newAtom->y,newAtom->z,newAtom->occ,newAtom->q);	


	return 0;
}

/*******************************************************************************
* This function reads the atomic position and element data for a single atom
* from a .cssr file.  The atomic positions are given in reduced coordinates.
* Calling this function several times will read one atom after another from the 
* file.
* A ReadLine error will occur, if the end of file is reached prematurely.
*******************************************************************************/
int readNextCSSRAtom(atom *newAtom,int flag, char *fileName) {
	static FILE *fpNextCSSR = NULL;
	static char buf[NCMAX];
	int count;
	static char element[32],s1[32],s2[32],s3[32];
	double dw;

	if (flag < 0) {
		if (fpNextCSSR != NULL) fclose(fpNextCSSR);
		fpNextCSSR = NULL;
		return 0;
	}

	if ( fpNextCSSR == NULL) 
	{
		fpNextCSSR = fopen(fileName, "r");

		if (fpNextCSSR == NULL) 
		{
			printf("Cannot open file %s\n",fileName);
			exit( 0 );
		}
		// skip the header:
		ReadLine( fpNextCSSR, buf, NCMAX, "in ReadXYZcoord");
		ReadLine( fpNextCSSR, buf, NCMAX, "in ReadXYZcoord");
		ReadLine( fpNextCSSR, buf, NCMAX, "in ReadXYZcoord");
		ReadLine( fpNextCSSR, buf, NCMAX, "in ReadXYZcoord");
	}

	ReadLine(fpNextCSSR, buf, NCMAX, "in ReadXYZcoord()");
	/* for Si */
	/*    dw = 0.444; */
	sscanf(buf,"%d %s %s %s %s %d %d %d %d %d %d %d %d %lf",
		&count,element,s1,s2,s3,&count,&count,
		&count,&count,&count,&count,&count,&count,&dw);

	newAtom->x = atof(s1);
	newAtom->y = atof(s2);
	newAtom->z = atof(s3);
	newAtom->occ = 1.0;
	newAtom->Znum = getZNumber(element); 
	newAtom->dw = dw;

	/*
	switch (newAtom->Znum) {
	case 14:  newAtom->dw = 0.45f;
	break;
	case 29:  newAtom->dw = 0.21f;
	break;
	default: newAtom->dw = 0.21f;
	}
	*/
	return 0;
}

// #define NCINMAX 500
// #define NPARAM	64    /* number of parameters */

//MCS - why do we have this function and initMuls in stem3.cpp?
MULS initMu() {
	MULS muls;
	int sCount,i,slices = 2;
	char waveFile[32];
	char *waveFileBase = "w";

	muls.slices = slices;

	/* general setup: */
	muls.lpartl = 0;

	// muls.wave = NULL;
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
	muls.cz = NULL;  // (float_t *)malloc(muls.slices*sizeof(float_t));

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
	return muls;
}
// #undef NCINMAX 500
// #undef NPARAM	64    /* number of parameters */

int readCubicCFG(double **pos,double **dw, int **Znums, double *ax,double *by,double *cz, 
				 double ctiltx, double ctilty) {
					 atom *atoms;
					 int Natom;
					 int j;
                  
					 MULS mu = initMu();

					 mu.atomKinds = 0;
					 mu.Znums = NULL;
					 mu.tds = 0;
					 mu.u2 = NULL;
					 mu.nCellX = 1;
					 mu.nCellY = 1;
					 mu.nCellZ = 1;
					 mu.ctiltx = ctiltx; mu.ctilty = ctilty; mu.ctiltz = 0;
					 mu.cubex = 0; mu.cubey = 0; mu.cubez = 0;

					 atoms = readUnitCell(&Natom,"demo.cfg",&mu,1);
					 /*  
					 *Znums = (int *)fftwf_malloc(Natom*sizeof(int));
					 *pos   = (double *)fftwf_malloc(Natom*3*sizeof(double));
					 *dw    = (double *)fftwf_malloc(Natom*sizeof(double));
					 */
					 *Znums = (int *)malloc(Natom*sizeof(int));
					 *pos   = (double *)malloc(Natom*3*sizeof(double));
					 *dw    = (double *)malloc(Natom*sizeof(double));
					 *ax = mu.ax;
					 *by = mu.by;
					 *cz = mu.c;
					 for (j=0;j<Natom;j++) {
						 (*Znums)[j] = atoms[j].Znum;
						 (*dw)[j] = 0; //  atoms[j].dw;
						 (*pos)[3*j+0] = atoms[j].y/mu.by;
						 (*pos)[3*j+1] = atoms[j].x/mu.ax;
						 (*pos)[3*j+2] = 1.0-atoms[j].z/mu.c;
						 if (j<10) {
							 printf("%2d: (%g %g %g) %d [%g]\n",j,(*pos)[3*j],(*pos)[3*j+1],(*pos)[3*j+2],(*Znums)[j],(*dw)[j]);
						 }
					 }
					 free(atoms);
					 return Natom;
}

/*******************************************************
* This function reads in a .cssr or .pdb file and fills a 
* atom struct array with the files data
* .cssr: Cerius data format, also supported by Cerius2
* .cfg:  MD Simulations Extended Configuration format, 
* supported by AtomEye
* atomic positions have to be in FRACTIONAL coordinates!!!
******************************************************/
#define FORMAT_UNKNOWN 0
#define FORMAT_CSSR 1
#define FORMAT_CFG 2
#define FORMAT_PDB 3
#define FORMAT_XYZ 4
#define FORMAT_DAT 5

////////////////////////////////////////////////////////////////////////
// replicateUnitCell
// 
// Replicates the unit cell NcellX x NCellY x NCellZ times
// applies phonon displacement and removes vacancies and atoms appearing
// on same position:
// ncoord is the number of atom positions that has already been read.
// memory for the whole atom-array of size natom has already been allocated
// but the sites beyond natom are still empty.
void replicateUnitCell(int ncoord,int *natom,MULS *muls,atom* atoms,int handleVacancies) {
	int i,j,i2,jChoice,ncx,ncy,ncz,icx,icy,icz,jz,jCell,jequal,jVac;
	int 	atomKinds = 0;
	double totOcc;
	double choice,lastOcc;
	double *u;
	// seed for random number generation
	static long idum = -1;

	ncx = muls->nCellX;
	ncy = muls->nCellY;
	ncz = muls->nCellZ;
	u = (double *)malloc(3*sizeof(double));

	atomKinds = muls->atomKinds;
	//////////////////////////////////////////////////////////////////////////////
	// Look for atoms which share the same position:
	jVac = 0;  // no atoms have been removed yet
	for (i=ncoord-1;i>=0;) {

		////////////////
		if ((handleVacancies) && (atoms[i].Znum > 0)) {
			totOcc = atoms[i].occ;
			for (jequal=i-1;jequal>=0;jequal--) {
				// if there is anothe ratom that comes close to within 0.1*sqrt(3) A we will increase 
				// the total occupany and the counter jequal.
				if ((fabs(atoms[i].x-atoms[jequal].x) < 1e-6) && (fabs(atoms[i].y-atoms[jequal].y) < 1e-6) && (fabs(atoms[i].z-atoms[jequal].z) < 1e-6)) {
					totOcc += atoms[jequal].occ;
				}
				else break;
			} // jequal-loop
		}
		else {
			jequal = i-1;
			totOcc = 1;
			// Keep a record of the kinds of atoms we are reading
		}
		if (jequal == i-1) {
			for (jz=0;jz<atomKinds;jz++)	if (muls->Znums[jz] == atoms[i].Znum) break;
		}

		////////////////
		memset(u,0,3*sizeof(double));
		/* replicate unit cell ncx,y,z times: */
		/* We have to start with the last atoms first, because once we added the displacements 
		* to the original unit cell (icx=icy=icz=0), we cannot use those positions			
		* as unit cell coordinates for the other atoms anymore
		*/
		// printf("Will start phonon displacement (%f)\n",muls->tds,muls->temperature);
		// for (jz=0;jz<muls->atomKinds;jz++)	if (atoms[i].Znum == muls->Znums[jz]) break;

		for (icx=ncx-1;icx>=0;icx--) {
			for (icy=ncy-1;icy>=0;icy--) {
				for (icz=ncz-1;icz>=0;icz--) {
					jCell = (icz+icy*ncz+icx*ncy*ncz)*ncoord;
					j = jCell+i;
					/* We will also add the phonon displacement to the atomic positions now: */
					atoms[j].dw = atoms[i].dw;
					atoms[j].occ = atoms[i].occ;
					atoms[j].q = atoms[i].q;
					atoms[j].Znum = atoms[i].Znum; 

					// Now is the time to remove atoms that are on the same position or could be vacancies:
					// if we encountered atoms in the same position, or the occupancy of the current atom is not 1, then
					// do something about it:
					jChoice = i;
					if ((totOcc < 1) || (jequal < i-1)) { // found atoms at equal positions or an occupancy less than 1!
						// ran1 returns a uniform random deviate between 0.0 and 1.0 exclusive of the endpoint values. 
						// 
						// if the total occupancy is less than 1 -> make sure we keep this
						// if the total occupancy is greater than 1 (unphysical) -> rescale all partial occupancies!
						if (totOcc < 1.0) choice = ran1(&idum);   
						else choice = totOcc*ran1(&idum);
						// printf("Choice: %g %g %d, %d %d\n",totOcc,choice,j,i,jequal);
						lastOcc = 0;
						for (i2=i;i2>jequal;i2--) {
							atoms[jCell+i2].dw = atoms[i2].dw;
							atoms[jCell+i2].occ = atoms[i2].occ;
							atoms[jCell+i2].q = atoms[i2].q;
							atoms[jCell+i2].Znum = atoms[i2].Znum; 

							// if choice does not match the current atom:
							// choice will never be 0 or 1(*totOcc) 
							if ((choice <lastOcc) || (choice >=lastOcc+atoms[i2].occ)) {
								// printf("Removing atom %d, Z=%d\n",jCell+i2,atoms[jCell+i2].Znum);
								atoms[jCell+i2].Znum =  0;  // vacancy
								jVac++;
							}
							else {
								jChoice = i2;
							}
							lastOcc += atoms[i2].occ;
						}

						// Keep a record of the kinds of atoms we are reading
						for (jz=0;jz<atomKinds;jz++) {
							if (muls->Znums[jz] == atoms[jChoice].Znum) break;
						}
					}
					// printf("i2=%d, %d (%d) [%g %g %g]\n",i2,jequal,jz,atoms[jequal].x,atoms[jequal].y,atoms[jequal].z);

					// this function does nothing, if muls->tds == 0
					// if (j % 5 == 0) printf("atomKinds: %d (jz = %d, %d)\n",atomKinds,jz,atoms[jChoice].Znum);
					phononDisplacement(u,muls,jChoice,icx,icy,icz,j,atoms[jChoice].dw,*natom,jz);
					// printf("atomKinds: %d (jz = %d, %d)\n",atomKinds,jz,atoms[jChoice].Znum);

					for (i2=i;i2>jequal;i2--) {
						atoms[jCell+i2].x = atoms[i2].x+icx+u[0];
						atoms[jCell+i2].y = atoms[i2].y+icy+u[1];
						atoms[jCell+i2].z = atoms[i2].z+icz+u[2];
					}
				}  // for (icz=ncz-1;icz>=0;icz--)
			} // for (icy=ncy-1;icy>=0;icy--) 
		} // for (icx=ncx-1;icx>=0;icx--)
		i=jequal;
	} // for (i=ncoord-1;i>=0;)
	if ((jVac > 0 ) &&(muls->printLevel)) printf("Removed %d atoms because of occupancies < 1 or multiple atoms in the same place\n",jVac);

}


// #define printf mexPrintf
//
// This function reads the atomic positions from fileName and also adds 
// Thermal displacements to their positions, if muls.tds is turned on.
atom *readUnitCell(int *natom,char *fileName,MULS *muls, int handleVacancies) {
	int printFlag = 1;
	// char buf[NCMAX], *str,element[16];
	// FILE *fp;
	// float_t alpha,beta,gamma;
	int ncoord=0,ncx,ncy,ncz,icx,icy,icz,jz;
	// float_t dw,occ,dx,dy,dz,r;
	int i,i2,j,format=FORMAT_UNKNOWN,ix,iy,iz,atomKinds=0;
	// char s1[16],s2[16],s3[16];
	double boxXmin=0,boxXmax=0,boxYmin=0,boxYmax=0,boxZmin=0,boxZmax=0;
	double boxCenterX,boxCenterY,boxCenterZ,boxCenterXrot,boxCenterYrot,boxCenterZrot,bcX,bcY,bcZ;
	double x,y,z,totOcc;
	double choice,lastOcc;
	double *u = NULL;
	double **Mm = NULL;
	static atom *atoms = NULL;
	static int ncoord_old = 0;

	printFlag = muls->printLevel;

	if (Mm == NULL) {
		Mm = double2D(3,3,"Mm");
		memset(Mm[0],0,9*sizeof(double));
		muls->Mm = Mm;
		u = (double *)malloc(3*sizeof(double));
	}

	/* figure out, whether we have  cssr, pdb, or cfg */
	if (strstr(fileName,".cssr") == fileName+strlen(fileName)-5) {
		format = FORMAT_CSSR;
		ncoord = readCSSRCellParams(muls,Mm,fileName);
	}
	if (strstr(fileName,".pdb") == fileName+strlen(fileName)-4) {
		format = FORMAT_PDB;
		printf("Cannot read Protein Data Bank format yet - sorry!\n");
		return NULL;
	}
	if (strstr(fileName,".xyz") == fileName+strlen(fileName)-4) {
		format = FORMAT_XYZ;
		printf("Cannot read XYZ format yet - sorry!\n");
		return NULL;
	}
	if (strstr(fileName,".cfg") == fileName+strlen(fileName)-4) {
		format = FORMAT_CFG;
		ncoord = readCFGCellParams(muls,Mm,fileName);
		// return readCFGUnitCell(natom,fileName,muls);
	}
	if (strstr(fileName,".dat") == fileName+strlen(fileName)-4) {
		format = FORMAT_DAT;
		ncoord = readDATCellParams(muls,Mm,fileName);
		// return readCFGUnitCell(natom,fileName,muls);
	}
	if (format == FORMAT_UNKNOWN) {
		printf("Cannot read anything else than .cssr, .cfg, or .pdb files (%s)!\n",fileName);
		return NULL;
	}

	if (ncoord == 0) {
		printf("Error reading configuration file %s - ncoord =0\n",fileName);
		return NULL;
	}


	ncx = muls->nCellX;
	ncy = muls->nCellY;
	ncz = muls->nCellZ;

	if (printFlag) {
		printf("Lattice parameters: ax=%g by=%g cz=%g (%d atoms)\n",
			muls->ax,muls->by,muls->c,ncoord);

		if ((muls->cubex == 0) || (muls->cubey == 0) || (muls->cubez == 0))	 
			printf("Size of Super-lattice: ax=%g by=%g cz=%g (%d x %d x %d)\n",
			muls->ax*ncx,muls->by*ncy,muls->c*ncz,ncx,ncy,ncz);
		else
			printf("Size of Cube: ax=%g by=%g cz=%g\n",
			muls->cubex,muls->cubey,muls->cubez);
	}
	/************************************************************
	* now that we know how many coordinates there are
	* allocate the arrays 
	************************************************************/
	/*
	printf("%d atoms (%d %d %d) (%d %d %d)\n",ncoord,ncx,ncy,ncz,
	muls->nCellX,muls->nCellY,muls->nCellZ);
	*/
	if (ncoord_old != ncoord) {
		if (atoms != NULL) free(atoms);
		atoms = (atom *)malloc(ncoord*sizeof(atom)*ncx*ncy*ncz);
		ncoord_old = ncoord;
	}
	if (atoms == NULL) {
		printf("Could not allocate memory for atoms!\n");
		return NULL;
	}
	*natom = ncoord*ncx*ncy*ncz;

	/*
	if (muls->atomKinds < 1) {
	muls->atomKinds = 1;
	muls->Znums = (int *)malloc(muls->atomKinds*sizeof(int));
	}
	*/
	atomKinds = 0;

	/***********************************************************
	* Read actual Data
	***********************************************************/
	for(jz=0,i=ncoord-1; i>=0; i--) {



		switch (format) {
		case FORMAT_CFG: 
		if (readNextCFGAtom(atoms+i,0,fileName) < 0) {
			printf("number of atoms does not agree with atoms in file!\n");
			return NULL;
		}
		break;
		case FORMAT_DAT: 

		if (readNextDATAtom(atoms+i,0,fileName) < 0) {
			printf("number of atoms does not agree with atoms in file!\n");
			return NULL;
		}
		break;

		case FORMAT_CSSR:
		readNextCSSRAtom(atoms+i,0,fileName);
		break;
		default: return NULL;
		}

		if((atoms[i].Znum < 1 ) || (atoms[i].Znum > NZMAX)) {
			/* for (j=ncoord-1;j>=i;j--)
			printf("%2d: %d (%g,%g,%g)\n",j,atoms[j].Znum,atoms[j].x,
			atoms[j].y,atoms[j].z);
			*/
			printf("Error: bad atomic number %d in file %s (atom %d [%d: %g %g %g])\n",
				atoms[i].Znum,fileName,i,atoms[i].Znum,atoms[i].x,atoms[i].y,atoms[i].z);
			return NULL;
		}

		// Keep a record of the kinds of atoms we are reading
		for (jz=0;jz<atomKinds;jz++)	if (muls->Znums[jz] == atoms[i].Znum) break;
		// allocate more memory, if there is a new element
		if (jz == atomKinds) {
			atomKinds++;
			if (atomKinds > muls->atomKinds) {
				muls->Znums = (int *)realloc(muls->Znums,atomKinds*sizeof(int));
				muls->atomKinds = atomKinds;
				// printf("%d kinds (%d)\n",atomKinds,atoms[i].Znum);
			}  
			muls->Znums[jz] = atoms[i].Znum;
		}


	} // for 1=ncoord-1:-1:0  - we've just read all the atoms.
	if (muls->tds) {
		if (muls->u2 == NULL) {
			// printf("AtomKinds: %d\n",muls->atomKinds);
			muls->u2 = (double *)malloc(atomKinds*sizeof(double));
			memset(muls->u2,0,atomKinds*sizeof(double));
		}
		if (muls->u2avg == NULL) {
			muls->u2avg = (double *)malloc(muls->atomKinds*sizeof(double));
			memset(muls->u2avg,0,muls->atomKinds*sizeof(double));
		}
	}

		////////////////////////////////////////////////////////////////
	// Close the file for further reading, and restore file pointer 
	switch (format) {
		case FORMAT_CFG: 
			readNextCFGAtom(NULL,-1,NULL);
			break;
		case FORMAT_DAT: 
			readNextDATAtom(NULL,-1,NULL);
			break;
		case FORMAT_CSSR:
			readNextCSSRAtom(NULL,-1,NULL);
			break;
	}

	// First, we will sort the atoms by position:
	if (handleVacancies) {
		qsort((void *)atoms,ncoord,sizeof(atom),atomCompareZYX);
	}


	/////////////////////////////////////////////////////////////////
	// Compute the phonon displacement and remove atoms which appear 
	// twice or have some probability to be vacancies:
	if ((muls->cubex > 0) && (muls->cubey > 0) && (muls->cubez > 0)) {
		/* at this point the atoms should have fractional coordinates */
		// printf("Entering tiltBoxed\n");
		atoms = tiltBoxed(ncoord,natom,muls,atoms,handleVacancies);
		// printf("ncoord: %d, natom: %d\n",ncoord,*natom);
	}
	else {  // work in NCell mode
		// atoms are in fractional coordinates so far, we need to convert them to 
		// add the phonon displacement in this condition, because there we can 
		// actually do the correct Eigenmode treatment.
		// but we will probably just do Einstein vibrations anyway:
		replicateUnitCell(ncoord,natom,muls,atoms,handleVacancies);
		/**************************************************************
		* now, after we read all of the important coefficients, we
		* need to decide if this is workable
		**************************************************************/
		*natom = ncoord*ncx*ncy*ncz;
		if (1) { // ((Mm[0][0]*Mm[1][1]*Mm[2][2] == 0) || (Mm[0][1]!=0)|| (Mm[0][2]!=0)|| (Mm[1][0]!=0)|| (Mm[1][2]!=0)|| (Mm[2][0]!=0)|| (Mm[2][1]!=0)) {
				// printf("Lattice is not orthogonal, or rotated\n");
				for(i=0;i<*natom;i++) {
					/*
					x = Mm[0][0]*atoms[i].x+Mm[0][1]*atoms[i].y+Mm[0][2]*atoms[i].z;
					y = Mm[1][0]*atoms[i].x+Mm[1][1]*atoms[i].y+Mm[1][2]*atoms[i].z;
					z = Mm[2][0]*atoms[i].x+Mm[2][1]*atoms[i].y+Mm[2][2]*atoms[i].z;
					*/

					// This converts also to cartesian coordinates
					x = Mm[0][0]*atoms[i].x+Mm[1][0]*atoms[i].y+Mm[2][0]*atoms[i].z;
					y = Mm[0][1]*atoms[i].x+Mm[1][1]*atoms[i].y+Mm[2][1]*atoms[i].z;
					z = Mm[0][2]*atoms[i].x+Mm[1][2]*atoms[i].y+Mm[2][2]*atoms[i].z;

					atoms[i].x = x;
					atoms[i].y = y;
					atoms[i].z = z;

				}      
		}
		/**************************************************************
		* Converting to cartesian coordinates
		*************************************************************/
		else { 
			for(i=0;i<*natom;i++) {
				atoms[i].x *= muls->ax; 
				atoms[i].y *= muls->by; 
				atoms[i].z *= muls->c;
			}		 
		}
		// Now we have all the cartesian coordinates of all the atoms!
		muls->ax *= ncx;
		muls->by *= ncy;
		muls->c  *= ncz;

		/***************************************************************
		* Now let us tilt around the center of the full crystal
		*/   
			
		bcX = ncx/2.0;
		bcY = ncy/2.0;
		bcZ = ncz/2.0;
		u[0] = Mm[0][0]*bcX+Mm[1][0]*bcY+Mm[2][0]*bcZ;
		u[1] = Mm[0][1]*bcX+Mm[1][1]*bcY+Mm[2][1]*bcZ;
		u[2] = Mm[0][2]*bcX+Mm[1][2]*bcY+Mm[2][2]*bcZ;
		boxCenterX = u[0];
		boxCenterY = u[1];
		boxCenterZ = u[2];
		
		// rotateVect(u,u,muls->ctiltx,muls->ctilty,muls->ctiltz);  // simply applies rotation matrix
		// boxCenterXrot = u[0]; boxCenterYrot = u[1];	boxCenterZrot = u[2];
		
		// Determine the size of the (rotated) super cell
		for (icx=0;icx<=ncx;icx+=ncx) for (icy=0;icy<=ncy;icy+=ncy) for (icz=0;icz<=ncz;icz+=ncz) {
			u[0] = Mm[0][0]*(icx-bcX)+Mm[1][0]*(icy-bcY)+Mm[2][0]*(icz-bcZ);
			u[1] = Mm[0][1]*(icx-bcX)+Mm[1][1]*(icy-bcY)+Mm[2][1]*(icz-bcZ);
			u[2] = Mm[0][2]*(icx-bcX)+Mm[1][2]*(icy-bcY)+Mm[2][2]*(icz-bcZ);
			rotateVect(u,u,muls->ctiltx,muls->ctilty,muls->ctiltz);  // simply applies rotation matrix
			// x = u[0]+boxCenterXrot; y = u[1]+boxCenterYrot; z = u[2]+boxCenterZrot;
			x = u[0]+boxCenterX; y = u[1]+boxCenterY; z = u[2]+boxCenterZ;
			if ((icx == 0) && (icy == 0) && (icz == 0)) {
				boxXmin = boxXmax = x;
				boxYmin = boxYmax = y;
				boxZmin = boxZmax = z;
			}
			else {
				boxXmin = boxXmin>x ? x : boxXmin; boxXmax = boxXmax<x ? x : boxXmax; 
				boxYmin = boxYmin>y ? y : boxYmin; boxYmax = boxYmax<y ? y : boxYmax; 
				boxZmin = boxZmin>z ? z : boxZmin; boxZmax = boxZmax<z ? z : boxZmax; 
			}
		}

		// printf("(%f, %f, %f): %f .. %f, %f .. %f, %f .. %f\n",muls->ax,muls->by,muls->c,boxXmin,boxXmax,boxYmin,boxYmax,boxZmin,boxZmax);


		if ((muls->ctiltx != 0) || (muls->ctilty != 0) || (muls->ctiltz != 0)) {			
			for(i=0;i<(*natom);i++) {

				u[0] = atoms[i].x-boxCenterX; 
				u[1] = atoms[i].y-boxCenterY; 
				u[2] = atoms[i].z-boxCenterZ; 
				rotateVect(u,u,muls->ctiltx,muls->ctilty,muls->ctiltz);  // simply applies rotation matrix
				u[0] += boxCenterX;
				u[1] += boxCenterY; 
				u[2] += boxCenterZ; 
				atoms[i].x = u[0];
				atoms[i].y = u[1]; 
				atoms[i].z = u[2]; 
				// boxXmin = boxXmin>u[0] ? u[0] : boxXmin; boxXmax = boxXmax<u[0] ? u[0] : boxXmax; 
				// boxYmin = boxYmin>u[1] ? u[1] : boxYmin; boxYmax = boxYmax<u[1] ? u[1] : boxYmax; 
				// boxZmin = boxZmin>u[2] ? u[2] : boxZmin; boxZmax = boxZmax<u[2] ? u[2] : boxZmax; 
			}
		} /* if tilts != 0 ... */

		for(i=0;i<(*natom);i++) {
			atoms[i].x-=boxXmin; 
			atoms[i].y-=boxYmin; 
			atoms[i].z-=boxZmin; 
		}
		muls->ax = boxXmax-boxXmin;
		muls->by = boxYmax-boxYmin;
		muls->c  = boxZmax-boxZmin;

		// printf("(%f, %f, %f): %f .. %f, %f .. %f, %f .. %f\n",muls->ax,muls->by,muls->c,boxXmin,boxXmax,boxYmin,boxYmax,boxZmin,boxZmax);
		/*******************************************************************
		* If one of the tilts was more than 30 degrees, we will re-assign 
		* the lattice constants ax, by, and c by boxing the sample with a box 
		******************************************************************/

	// Offset the atoms in x- and y-directions:
	// Do this after the rotation!
	if ((muls->xOffset != 0) || (muls->yOffset != 0)) {
		for(i=0;i<*natom;i++) {
			atoms[i].x += muls->xOffset; 
			atoms[i].y += muls->yOffset; 
		}		 
	}
	} // end of Ncell mode conversion to cartesian coords and tilting.
	// printf("Offset: (%f, %f)\n",muls->xOffset,muls->yOffset);

	/*
	printf("Atom box: (%g %g %g) .. (%g %g %g)\n",boxXmin,boxYmin,boxZmin,boxXmax,boxYmax,boxZmax);
	printf("ax: %g, bx: %g, c: %g\n",muls->ax,muls->by,muls->c);

	printf("%d atoms read from file <%s>, %d atoms in model (Tilt: x: %g mrad, y: %g mrad).\n",
	ncoord,fileName,*natom,muls->ctiltx,muls->ctilty);
	*/

	// initialize vibration amplitude counters:
	//////////////////////////////////////////////////////////////////////////////////////////



	return atoms;
} // end of readUnitCell



atom *tiltBoxed(int ncoord,int *natom, MULS *muls,atom *atoms,int handleVacancies) {
	int atomKinds = 0;
	int iatom,jVac,jequal,jChoice,i2,ix,iy,iz,atomCount = 0,atomSize;
	static double *axCell,*byCell,*czCell=NULL;
	static double **Mm = NULL, **Mminv = NULL, **MpRed = NULL, **MpRedInv = NULL;
	static double **MbPrim = NULL, **MbPrimInv = NULL, **MmOrig = NULL,**MmOrigInv=NULL;
	static double **a = NULL,**aOrig = NULL,**b= NULL,**bfloor=NULL,**blat=NULL;
	static double *uf;
	static int oldAtomSize = 0;
	double x,y,z,dx,dy,dz; 
	double totOcc,lastOcc,choice;
	atom *unitAtoms,newAtom;
	int nxmin,nxmax,nymin,nymax,nzmin,nzmax,jz;
	// static FILE *fpPhonon = NULL;
	// FILE *fpu2;
	// static int Nk, Ns;     // number of k-vectors and atoms per primitive unit cell
	// static float **bPrim;   // basis of primitive unit cell
	// static float **posPrim; // atomic positions per primitive basis
	// static float *massPrim;  // masses for every atom in primitive basis
	// static float **omega;  // array of eigenvalues for every k-vector 
	// static fftw_complex ***eigVecs;  // array of eigenvectors for every k-vector
	// static float **kVecs;    // array for Nk 3-dim k-vectors
	// static double **q1=NULL, **q2=NULL;
	int Ncells;
	// double kR,kRi,kRr;
	// FILE *fp;
	//static double u2=0;
	//static int u2Count = 0;
	// static long iseed=0;
	static double *u;
	static long idum = -1;


	// if (iseed == 0) iseed = -(long) time( NULL );
	Ncells = muls->nCellX * muls->nCellY * muls->nCellZ;

	/* calculate maximum length in supercell box, which is naturally 
	* the room diagonal:

	maxLength = sqrt(muls->cubex*muls->cubex+
	muls->cubey*muls->cubey+
	muls->cubez*muls->cubez);
	*/

	if (Mm == NULL) {
		MmOrig		= double2D(3,3,"MmOrig");
		MmOrigInv	= double2D(3,3,"MmOrigInv");
		MbPrim		= double2D(3,3,"MbPrim");	// double version of primitive lattice basis 
		MbPrimInv	= double2D(3,3,"MbPrim"); // double version of inverse primitive lattice basis 
		MpRed		= double2D(3,3,"MpRed");    /* conversion lattice to obtain red. prim. coords 
												* from reduced cubic rect.
												*/
		MpRedInv	= double2D(3,3,"MpRedInv");    /* conversion lattice to obtain red. cub. coords 
												   * from reduced primitive lattice coords
												   */
		Mm			= double2D(3,3,"Mm");
		Mminv		= double2D(3,3,"Mminv");
		axCell = Mm[0]; byCell = Mm[1]; czCell = Mm[2];
		a			= double2D(1,3,"a");
		aOrig		= double2D(1,3,"aOrig");
		b			= double2D(1,3,"b");
		bfloor		= double2D(1,3,"bfloor");
		blat		= double2D(1,3,"blat");
		uf			= (double *)malloc(3*sizeof(double));
		u			= (double *)malloc(3*sizeof(double));
	}


	dx = 0; dy = 0; dz = 0;
	dx = muls->xOffset;
	dy = muls->yOffset;
	/* find the rotated unit cell vectors .. 
	* muls does still hold the single unit cell vectors in ax,by, and c
	*/
	// makeCellVectMuls(muls, axCell, byCell, czCell);
	// We need to copy the transpose of muls->Mm to Mm.
	// we therefore cannot use the following command:
	// memcpy(Mm[0],muls->Mm[0],3*3*sizeof(double));
	for (ix=0;ix<3;ix++) for (iy=0;iy<3;iy++) Mm[ix][iy]=muls->Mm[iy][ix];

	memcpy(MmOrig[0],Mm[0],3*3*sizeof(double));
	inverse_3x3(MmOrigInv[0],MmOrig[0]);
	/* remember that the angles are in rad: */
	rotateMatrix(Mm[0],Mm[0],muls->ctiltx,muls->ctilty,muls->ctiltz);
	/*
	// This is wrong, because it implements Mrot*(Mm'):
	rotateVect(axCell,axCell,muls->ctiltx,muls->ctilty,muls->ctiltz);
	rotateVect(byCell,byCell,muls->ctiltx,muls->ctilty,muls->ctiltz);
	rotateVect(czCell,czCell,muls->ctiltx,muls->ctilty,muls->ctiltz);
	*/
	inverse_3x3(Mminv[0],Mm[0]);  // computes Mminv from Mm!
	/* find out how far we will have to go in unit of unit cell vectors.
	* when creating the supercell by checking the number of unit cell vectors 
	* necessary to reach every corner of the supercell box.
	*/
	// showMatrix(MmOrig,3,3,"Morig");
	// printf("%d %d\n",(int)Mm, (int)MmOrig);
	memset(a[0],0,3*sizeof(double));
	// matrixProduct(a,1,3,Mminv,3,3,b);
	matrixProduct(Mminv,3,3,a,3,1,b);
	// showMatrix(Mm,3,3,"M");
	// showMatrix(Mminv,3,3,"M");
	nxmin = nxmax = (int)floor(b[0][0]-dx); 
	nymin = nymax = (int)floor(b[0][1]-dy); 
	nzmin = nzmax = (int)floor(b[0][2]-dz);
	for (ix=0;ix<=1;ix++) for (iy=0;iy<=1;iy++)	for (iz=0;iz<=1;iz++) {
		a[0][0]=ix*muls->cubex-dx; a[0][1]=iy*muls->cubey-dy; a[0][2]=iz*muls->cubez-dz;

		// matrixProduct(a,1,3,Mminv,3,3,b);
		matrixProduct(Mminv,3,3,a,3,1,b);

		// showMatrix(b,1,3,"b");
		if (nxmin > (int)floor(b[0][0])) nxmin=(int)floor(b[0][0]);
		if (nxmax < (int)ceil( b[0][0])) nxmax=(int)ceil( b[0][0]);
		if (nymin > (int)floor(b[0][1])) nymin=(int)floor(b[0][1]);
		if (nymax < (int)ceil( b[0][1])) nymax=(int)ceil( b[0][1]);
		if (nzmin > (int)floor(b[0][2])) nzmin=(int)floor(b[0][2]);
		if (nzmax < (int)ceil( b[0][2])) nzmax=(int)ceil( b[0][2]);	  
	}

	// nxmin--;nxmax++;nymin--;nymax++;nzmin--;nzmax++;
	unitAtoms = (atom *)malloc(ncoord*sizeof(atom));
	memcpy(unitAtoms,atoms,ncoord*sizeof(atom));
	atomSize = (1+(nxmax-nxmin)*(nymax-nymin)*(nzmax-nzmin)*ncoord);
	if (atomSize != oldAtomSize) {
		atoms = (atom *)realloc(atoms,atomSize*sizeof(atom));
		oldAtomSize = atomSize;
	}
	// showMatrix(Mm,3,3,"Mm");
	// showMatrix(Mminv,3,3,"Mminv");
	// printf("Range: (%d..%d, %d..%d, %d..%d)\n",
	// nxmin,nxmax,nymin,nymax,nzmin,nzmax);

	atomCount = 0;  
	jVac = 0;
	memset(u,0,3*sizeof(double));
	for (iatom=0;iatom<ncoord;) {
		// printf("%d: (%g %g %g) %d\n",iatom,unitAtoms[iatom].x,unitAtoms[iatom].y,
		//   unitAtoms[iatom].z,unitAtoms[iatom].Znum);
		memcpy(&newAtom,unitAtoms+iatom,sizeof(atom));
		for (jz=0;jz<muls->atomKinds;jz++)	if (muls->Znums[jz] == newAtom.Znum) break;
		// allocate more memory, if there is a new element
		/*
		if (jz == atomKinds) {
			atomKinds++;
			if (atomKinds > muls->atomKinds) {
				muls->Znums = (int *)realloc(muls->Znums,atomKinds*sizeof(int));
				muls->atomKinds = atomKinds;
				// printf("%d kinds (%d)\n",atomKinds,atoms[i].Znum);
			}  
			muls->Znums[jz] = newAtom.Znum;
		}
		*/
		/////////////////////////////////////////////////////
		// look for atoms at equal position
		if ((handleVacancies) && (newAtom.Znum > 0)) {
			totOcc = newAtom.occ;
			for (jequal=iatom+1;jequal<ncoord;jequal++) {
				// if there is anothe ratom that comes close to within 0.1*sqrt(3) A we will increase 
				// the total occupany and the counter jequal.
				if ((fabs(newAtom.x-unitAtoms[jequal].x) < 1e-6) && (fabs(newAtom.y-unitAtoms[jequal].y) < 1e-6) && (fabs(newAtom.z-unitAtoms[jequal].z) < 1e-6)) {
					totOcc += unitAtoms[jequal].occ;
				}
				else break;
			} // jequal-loop
		}
		else {
			jequal = iatom+1;
			totOcc = 1;
		}



		// printf("%d: %d\n",atomCount,jz);
		for (ix=nxmin;ix<=nxmax;ix++) {
			for (iy=nymin;iy<=nymax;iy++) {
				for (iz=nzmin;iz<=nzmax;iz++) {
					// atom position in cubic reduced coordinates: 
					aOrig[0][0] = ix+newAtom.x; aOrig[0][1] = iy+newAtom.y; aOrig[0][2] = iz+newAtom.z;

					// Now is the time to remove atoms that are on the same position or could be vacancies:
					// if we encountered atoms in the same position, or the occupancy of the current atom is not 1, then
					// do something about it:
					// All we need to decide is whether to include the atom at all (if totOcc < 1
					// of which of the atoms at equal positions to include
					jChoice = iatom;  // This will be the atom we wil use.
					if ((totOcc < 1) || (jequal > iatom+1)) { // found atoms at equal positions or an occupancy less than 1!
						// ran1 returns a uniform random deviate between 0.0 and 1.0 exclusive of the endpoint values. 
						// 
						// if the total occupancy is less than 1 -> make sure we keep this
						// if the total occupancy is greater than 1 (unphysical) -> rescale all partial occupancies!
						if (totOcc < 1.0) choice = ran1(&idum);   
						else choice = totOcc*ran1(&idum);
						// printf("Choice: %g %g %d, %d %d\n",totOcc,choice,j,i,jequal);
						lastOcc = 0;
						for (i2=iatom;i2<jequal;i2++) {
							// atoms[atomCount].Znum = unitAtoms[i2].Znum; 
							// if choice does not match the current atom:
							// choice will never be 0 or 1(*totOcc) 
							if ((choice <lastOcc) || (choice >=lastOcc+unitAtoms[i2].occ)) {
								// printf("Removing atom %d, Z=%d\n",jCell+i2,atoms[jCell+i2].Znum);
								// atoms[atomCount].Znum =  0;  // vacancy
								jVac++;
							}
							else {
								jChoice = i2;
							}
							lastOcc += unitAtoms[i2].occ;
						}
						// printf("Keeping atom %d (%d), Z=%d\n",jChoice,iatom,unitAtoms[jChoice].Znum);
					}
					// if (jChoice != iatom) memcpy(&newAtom,unitAtoms+jChoice,sizeof(atom));
					if (jChoice != iatom) {
						for (jz=0;jz<muls->atomKinds;jz++)	if (muls->Znums[jz] == unitAtoms[jChoice].Znum) break;
					}




					// here we need to call phononDisplacement:
					// phononDisplacement(u,muls,iatom,ix,iy,iz,atomCount,atoms[i].dw,*natom,atoms[i].Znum);
					if (muls->Einstein == 1) {
						// phononDisplacement(u,muls,iatom,ix,iy,iz,1,newAtom.dw,10,newAtom.Znum);
						if (muls->tds) {
							phononDisplacement(u,muls,jChoice,ix,iy,iz,1,unitAtoms[jChoice].dw,atomSize,jz);
							a[0][0] = aOrig[0][0]+u[0]; a[0][1] = aOrig[0][1]+u[1]; a[0][2] = aOrig[0][2]+u[2];
						}
						else {
							a[0][0] = aOrig[0][0]; a[0][1] = aOrig[0][1]; a[0][2] = aOrig[0][2];
						}
					}
					else {
						printf("Cannot handle phonon-distribution mode for boxed sample yet - sorry!!\n");
						exit(0);
					}
					// matrixProduct(aOrig,1,3,Mm,3,3,b);
					matrixProduct(Mm,3,3,aOrig,3,1,b);

					// if (atomCount < 2) {showMatrix(a,1,3,"a");showMatrix(b,1,3,"b");}
					// b now contains atom positions in cartesian coordinates */
					x  = b[0][0]+dx; 
					y  = b[0][1]+dy; 
					z  = b[0][2]+dz; 
					if ((x >= 0) && (x <= muls->cubex) &&
						(y >= 0) && (y <= muls->cubey) &&
						(z >= 0) && (z <= muls->cubez)) {
							// matrixProduct(a,1,3,Mm,3,3,b);
							matrixProduct(Mm,3,3,a,3,1,b);
							atoms[atomCount].x		= b[0][0]+dx; 
							atoms[atomCount].y		= b[0][1]+dy; 
							atoms[atomCount].z		= b[0][2]+dz; 
							atoms[atomCount].dw		= unitAtoms[jChoice].dw;
							atoms[atomCount].occ	= unitAtoms[jChoice].occ;
							atoms[atomCount].q		= unitAtoms[jChoice].q;
							atoms[atomCount].Znum	= unitAtoms[jChoice].Znum;
								
							atomCount++;	
							/*
							if (unitAtoms[jChoice].Znum > 22)
								printf("Atomcount: %d, Z = %d\n",atomCount,unitAtoms[jChoice].Znum);
							*/
					}

				} /* iz ... */
			} /* iy ... */
		} /* ix ... */
		iatom = jequal;
	} /* iatom ... */
	if (muls->printLevel > 2) printf("Removed %d atoms because of multiple occupancy or occupancy < 1\n",jVac);
	muls->ax = muls->cubex;
	muls->by = muls->cubey;
	muls->c  = muls->cubez;
	*natom = atomCount;
	// call phononDisplacement again to update displacement data:
	phononDisplacement(u,muls,iatom,ix,iy,iz,0,newAtom.dw,*natom,jz);


	free(unitAtoms);
	return atoms;
}  // end of 'tiltBoxed(...)'


















/*******************************************************
* This function reads in a .cssr or .pdb file and fills a 
* atom struct array with the files data
* .cssr: Cerius data format, also supported by Cerius2
* .pdb:  Protein Data Bank format, supported by AtomEye
* atomic positions will be in FRACTIONAL coordinates!!!
******************************************************/

atom *readUnitCell_old(int *natom,char *fileName,MULS *muls) {
	char buf[NCMAX], *str,element[16];
	FILE *fpUnit;
	real alpha,beta,gamma;
	int ncoord,ncx,ncy,ncz,icx,icy,icz;
	real dw,occ,dx,dy,dz,r;
	int spaceGrp,i,count,j,type;
	static atom *atoms = NULL;
	char s1[32],s2[32],s3[32],dummy[32];
	static int ncoord_old = 0;
	int format = FORMAT_UNKNOWN;

	/* figure out, whether we have xyz, or cssr */
	if (strstr(fileName,".cssr") == fileName+strlen(fileName)-5)
		format = FORMAT_CSSR;
	if (strstr(fileName,".pdb") == fileName+strlen(fileName)-4)
		format = FORMAT_PDB;
	if (format == FORMAT_UNKNOWN) {
		printf("Cannot read anything else than .cssr or .pdb files (%s)!\n",
			fileName);
		exit(0);
	}

	printf("Format of input file: %d (%s)\n",format,format==FORMAT_CSSR ? "cssr" : "pdb");

	ncx = muls->nCellX;
	ncy = muls->nCellY;
	ncz = muls->nCellZ;
	/******************************************************************/
	fpUnit = fopen( fileName, "r" );
	if ( fpUnit == NULL ) {
		printf("Cannot open file %s\n",fileName);
		exit( 0 );
	}
	if (format == FORMAT_CSSR) {

		/* skip first line in first pass */
		ReadLine( fpUnit, buf, NCMAX, "in ReadXYZcoord" );
		ReadLine( fpUnit, buf, NCMAX, "in ReadXYZcoord" );

		/************************************************************
		* check angles, of cssr file (must be 90 degrees)
		*************************************************************/
		str = strchr(buf,'9');
		alpha = (real)atof(str++);
		beta = (real)atof(strchr(str++,'9'));
		gamma = (real)atof(strchr(str,'9'));
		muls->cAlpha = alpha;
		muls->cBeta  = beta;
		muls->cGamma = gamma;
		/*
		if ((alpha != 90.0) || (beta != 90.0) || (gamma!= 90.0)) {
		printf("angles are not 90 degrees! (%g, %g, %g) \n",
		alpha,beta,gamma);
		exit(0);
		}
		*/
		/* check the space group: */
		spaceGrp = atoi(strstr(buf,"SPGR =")+strlen("SPGR ="));
		if (spaceGrp != 1) {
			printf("cannot interpret space group %d\n",spaceGrp);
			exit(0);
		}
		ReadLine( fpUnit, buf, NCMAX, "in ReadXYZcoord" );
		ncoord = atoi(buf);
	}
	else if (format == FORMAT_PDB) {
		/*****************************************************************
		* File Format of pdb files:
		* HEADER    libAtoms:Config_save_as_pdb; 388 atoms
		* CRYST1   15.638   15.638   42.547  90.00  90.00  90.00 P 1           1
		* ATOM      1  C           1      11.729   7.819   0.000
		* ATOM      2  C           1      11.483   9.183  42.182
		*****************************************************************
		* We will look for the first string that is a non-zero number in 
		* the first line and make it the number of atoms.
		* The second line will determine a,b,c,alpha,beta,gamma
		****************************************************************/

		ReadLine( fpUnit, buf, NCMAX, "in ReadXYZcoord" );
		str = buf;
		ncoord = 0;
		do {
			ncoord = atoi(str);
			if (ncoord > 0) break;
		}
		while(((str=strtok(str," \t")) != NULL) && (ncoord == 0));
		if (ncoord == 0) {
			printf("Could not determine number of atmoms - quit\n");
			exit(0);
		}
		ReadLine( fpUnit, buf, NCMAX, "in ReadXYZcoord" );

		/************************************************************
		* check angles, of pdb file (must be 90 degrees)
		* This cannot work yet - please check !!!
		*************************************************************/
		str = buf;
		muls->ax = (real)atof(str);
		muls->by = (real)atof(str = strtok(buf," \t"));
		muls->c  = (real)atof(str = strtok(buf," \t"));
		alpha = (real)atof(str = strtok(buf," \t"));
		beta = (real)atof(str = strtok(buf," \t"));
		gamma = (real)atof(str = strtok(buf," \t"));
		muls->cAlpha = alpha;
		muls->cBeta  = beta;
		muls->cGamma = gamma;
		/*
		if ((alpha != 90.0) || (beta != 90.0) || (gamma!= 90.0)) {
		printf("angles are not 90 degrees! (%g, %g, %g) \n",
		alpha,beta,gamma);
		exit(0);
		}
		*/
	}
	else {
		printf("Wrong file format!\n");
		exit(0);
	}

	/************************************************************
	* now that we know how many coordinates there are
	* allocate the arrays 
	************************************************************/
	if (muls->nCellX == 0) {
		printf("Number of unit cell repetitions in X is 0, adjusted to 1/n");
		muls->nCellX = 1;
	}
	if (muls->nCellY == 0) {
		printf("Number of unit cell repetitions in Y is 0, adjusted to 1/n");
		muls->nCellY = 1;
	}
	if (muls->nCellZ == 0) {
		printf("Number of unit cell repetitions in Z is 0, adjusted to 1/n");
		muls->nCellZ = 1;
	}

	if (ncoord_old != ncoord) {
		if (atoms != NULL)
			free(atoms);
		atoms = (atom *)malloc(ncoord*sizeof(atom)*ncx*ncy*ncz);
		ncoord_old = ncoord;
	}
	if (atoms == NULL) {
		printf("Could not allocate memory for atoms!\n");
		exit(0);
	}
	if (format == FORMAT_CSSR) {
		fclose( fpUnit );
		fpUnit = fopen( fileName, "r" );
		if ( fpUnit == NULL ) {
			printf("Cannot open file %s\n",fileName);
			exit( 0 );
		}

		/***********************************************************
		* Read Header for .cssr and .pdb
		***********************************************************/
		ReadLine( fpUnit, buf, NCMAX, "in ReadXYZcoord" );
		/* printf(buf); */
		sscanf( buf, " %s %s %s",s1,s2,s3);
		/* printf("\n%s\n",s2); */
		muls->ax = atof(s1);
		muls->by = atof(s2);
		muls->c = atof(s3);
		ReadLine( fpUnit, buf, NCMAX, "in ReadXYZcoord" );
		ReadLine( fpUnit, buf, NCMAX, "in ReadXYZcoord" );
		ReadLine( fpUnit, buf, NCMAX, "in ReadXYZcoord" );
	} /* end of if FORMAT_CSSR */
	printf("Lattice parameters: ax=%g(%dx%g) by=%g(%dx%g) cz=%g(%dx%g)\n",
		muls->ax*ncx,ncx,muls->ax,muls->by*ncy,ncy,muls->by,
		muls->c*ncz,ncz,muls->c);

	/***********************************************************
	* Read actual Data
	***********************************************************/
	for( i=0; i<ncoord; i++) {
		ReadLine( fpUnit, buf, NCMAX, "in ReadXYZcoord()" );
		/* for Si */
		/*    dw = 0.444; */
		occ = 1.0F;
		if (format == FORMAT_CSSR) {
			sscanf(buf,"%d %s %s %s %s %d %d %d %d %d %d %d %d %g",
				&count,element,s1,s2,s3,&count,&count,
				&count,&count,&count,&count,&count,&count,&dw);
		}
		else if (format==FORMAT_PDB) {
			/* ATOM      1  C 1      11.729   7.819   0.000 */
			sscanf(buf,"%s %d %s %d %s %s %s",dummy,&count,element,&type,s1,s2,s3);
		}
		atoms[i].x = atof(s1);
		atoms[i].y = atof(s2);
		atoms[i].z = atof(s3);
		atoms[i].occ = occ;
		atoms[i].q = 0;
		atoms[i].Znum = getZNumber(element); 
		switch (atoms[i].Znum) {
	case 14:  atoms[i].dw = 0.45f;
		break;
	case 29:  atoms[i].dw = 0.21f;
		break;
	default: atoms[i].dw = 0.21f;
		}
		/* get wobble factor from debye waller factor 
		(*wobble)[i] = sqrt(*(*wobble+i)/(8*PI*PI));
		*/    
		if((atoms[i].Znum < 1 ) || (atoms[i].Znum > NZMAX)) {
			printf("Error: bad atomic number %d in file %s\n",
				atoms[i].Znum,fileName);
			exit(0);
		}
		if (format == FORMAT_CSSR) {
			/* replicate unit cell ncx,y,z times: */
			for (icx=0;icx<ncx;icx++) {
				for (icy=0;icy<ncy;icy++) {
					for (icz=0;icz<ncz;icz++) {
						j = (icz+icy*ncz+icx*ncy*ncz)*ncoord+i;
						atoms[j].x   = atoms[i].x+icx;
						atoms[j].y   = atoms[i].y+icy;
						atoms[j].z   = atoms[i].z+icz;
						atoms[j].dw  = atoms[i].dw;
						atoms[j].occ = atoms[i].occ;
						atoms[j].q   = atoms[i].q;
						atoms[j].Znum = atoms[i].Znum; 	  
					}
				}
			}
		}
		else {
			/* replicate unit cell ncx,y,z times: */
			for (icx=0;icx<ncx;icx++) {
				for (icy=0;icy<ncy;icy++) {
					for (icz=0;icz<ncz;icz++) {
						j = (icz+icy*ncz+icx*ncy*ncz)*ncoord+i;
						atoms[j].x = atoms[i].x+icx*muls->ax;
						atoms[j].y = atoms[i].y+icy*muls->by;
						atoms[j].z = atoms[i].z+icz*muls->c;
						atoms[j].dw = atoms[i].dw;
						atoms[j].occ = atoms[i].occ;
						atoms[j].q   = atoms[i].q;
						atoms[j].Znum = atoms[i].Znum; 	  
					}
				}
			}
		}
	}
	fclose( fpUnit );

	*natom = ncoord*ncx*ncy*ncz;

	/**************************************************************
	* Converting to cartesian coordinates
	if (format == FORMAT_CSSR) {
	for(i=0;i<*natom;i++) {
	atoms[i].x *= muls->ax; 
	atoms[i].y *= muls->by; 
	atoms[i].z *= muls->c;
	}	
	}	 
	*************************************************************/
	muls->ax *= ncx;
	muls->by *= ncy;
	muls->c  *= ncz;

	/***************************************************************
	* Now let us tilt around the center of the full crystal
	*/

	if ((muls->ctiltx != 0) || (muls->ctilty)!=0) {
		for(i=0;i<*natom;i++) {
			/* apply tilt around x-axis: */
			if (muls->ctiltx != 0) {
				dy =(atoms[i].y-muls->by/2.0);
				dz = (atoms[i].z-muls->c/2.0);
				r = sqrt(dz*dz+dy*dy);
				alpha = atan2(dz,dy)+muls->ctiltx;
				atoms[i].y = muls->by/2.0+r*cos(alpha); 
				atoms[i].z = muls->c/2.0+r*sin(alpha); 
			}
			if (muls->ctilty != 0) {
				dx =(atoms[i].x-muls->ax/2.0);
				dz = (atoms[i].z-muls->c/2.0);
				r = sqrt(dz*dz+dx*dx);
				alpha = atan2(dz,dx)+muls->ctilty;
				atoms[i].x = muls->ax/2.0+r*cos(alpha); 
				atoms[i].z = muls->c/2.0+r*sin(alpha); 
			}
		}
	}



	/*  natom = ReadXYZcoord(filein,1,1,1,
	&ax, &by, &cz, &Znum, &x, &y, &z, &occ, &wobble,
	description, NCMAX );
	*/



	printf("%d atoms read from file <%s>, %d atoms in model (Tilt: x: %g mrad, y: %g mrad).\n",
		ncoord,fileName,*natom,muls->ctiltx,muls->ctilty);
	return atoms;
}

/*--------------------- ReadLine() -----------------------*/
/*
read a full line from a file and 
return length of line read

to bad this looks like Pascal but its the easiest
way to read just whole line because fscanf() ignores
end of line characters

fpread = pointer to file
cMax = length of data buffer cRead
cRead = char[] buffer to read into
mesg = error message to print if not successful
*/
int ReadLine( FILE* fpRead, char* cRead, int cMax, const char *mesg )
{
	if( fgets( cRead, cMax, fpRead) == NULL ) {
		return 0;
		/*   printf("error reading input file: %s\n", mesg);
		exit( 0 );
		*/
	}
	return( strlen( cRead ) );

}  /* end ReadLine */



int getZNumber(char *element) {
	char *elem;
#ifndef WIN32
	char widows = '\n';
	widows = (char)((int)widows + 3);

	// fprintf(stderr,"getZnumber : element = *%s* \n",element); 
	element[2] = '\0'; 
	if ((atoi(element+1) != 0) || (element[1] == '\n')|| (element[1] == '\0') || element[1] == widows){
		element[1] = ' ';
		// fprintf(stderr,"getZnumber : element = *%s*  conversion \n",element); 
	} 
#else
	element[2] = '\0';
	if ((atoi(element+1) != 0) || (element[1] == '\n')|| (element[1] == '\0'))
		element[1] = ' ';
#endif
	if ((elem = strstr(elTable,element)) == NULL)
		return 0;	
	else
		return (int)(elem-elTable)/2+1;

}

#define CHARGE 0.0

void writeFrameWork(FILE *fpFrame,superCellBox superCell) {
	int i,id = 0,newId = 1;
	double charge=0.0;

	if ( fpFrame != NULL ) fprintf( fpFrame, "frame1 1 framework\n" );
	for (i=0;i<superCell.natoms;i++) {
		if (i==0) {
			if (idArray == NULL) {
				idArraySize = 2;
				idArray=(int *)malloc(idArraySize*sizeof(int));
				id = 0;
				idArrayPtr = 0;
			}
			else {
				/* check, if this Znum is already somewhere in the list, and if not,
				* add it to the list, extending its size, if necessary
				*/
				for (id=0;id<=idArrayPtr;id++) if (superCell.atoms[i].Znum == idArray[id]) break;
				if (id>idArrayPtr) {
					if (idArrayPtr >= idArraySize-1) {
						idArraySize *=2;
						idArray = (int *)realloc(idArray,idArraySize*sizeof(int));
					}
					idArray[++idArrayPtr] = superCell.atoms[i].Znum;
					newId = 1;
				}           
			}
			idArray[id] = superCell.atoms[i].Znum;
		}
		else {
			newId = 0;
			if (superCell.atoms[i].Znum != superCell.atoms[i-1].Znum) {
				/* check, if this Znum is already somewhere in the list, and if not,
				* add it to the list, extending its size, if necessary
				*/
				for (id=0;id<=idArrayPtr;id++) if (superCell.atoms[i].Znum == idArray[id]) break;
				if (id>idArrayPtr) {
					if (idArrayPtr >= idArraySize-1) {
						idArraySize *=2;
						idArray = (int *)realloc(idArray,idArraySize*sizeof(int));
					}
					idArray[++idArrayPtr] = superCell.atoms[i].Znum;
					newId = 1;
				}
			}
		}
		if (newId) {
			switch (superCell.atoms[i].Znum) {
	  case 14: charge = 4.0; break;
	  case 7:  charge = -3.0; break;
	  case 8:  charge = -2.0; break;
	  case 39: charge = 3.0; break;
	  default: charge = 0.0; 
			}

			if ( fpFrame != NULL ) fprintf( fpFrame, "%d %7.3f %7.3f %7.3f %6.3f %6.3f %c%c\n", id + 1, superCell.atoms[i].x,
				superCell.atoms[i].y,superCell.atoms[i].z,
				massArray[superCell.atoms[i].Znum-1],charge,
				elTable[2*superCell.atoms[i].Znum-2],elTable[2*superCell.atoms[i].Znum-1]);
		}
		else
		if ( fpFrame != NULL ) fprintf( fpFrame, "%d %7.3f %7.3f %7.3f\n", id + 1, superCell.atoms[i].x,
			superCell.atoms[i].y,superCell.atoms[i].z);
	}
}


#define A 0.0     // 975.857
#define B 43684.0   // 425379 
#define C 3.4483  // 5.04748

/* This function will write the amorphous data to the MD input file starting at 
* atoms nstart and ending just before atom nstop.   
*/
void writeAmorphous(FILE *fpAmorph,superCellBox superCell,int nstart,int nstop) {
	int i,j,id;
	int *idCountArray = NULL;
	double charge,b,x,y,z;
	// char elem[8];

	memset(chargeTable,0,MAX_MASS_INDEX*sizeof(double));

	printf("amorph: %d ..%d-1\n",nstart,nstop);

	if (idArray == NULL) {
		idArraySize = 2;
		idArray=(int *)malloc(idArraySize*sizeof(int));
		id = 0;
		idArrayPtr = 0;
		idArray[0] = superCell.atoms[nstart].Znum;
	}
	else {
		/* check, if this Znum is already somewhere in the list, and if not,
		* add it to the list, extending its size, if necessary
		*/
		for (id=0;id<=idArrayPtr;id++) if (superCell.atoms[nstart].Znum == idArray[id]) break;
		if (id>idArrayPtr) {
			if (idArrayPtr >= idArraySize-1) {
				idArraySize *=2;
				idArray = (int *)realloc(idArray,idArraySize*sizeof(int));
			}
			idArray[++idArrayPtr] = superCell.atoms[nstart].Znum;
		}           
	}
	idCountArray = (int *)malloc(idArraySize*sizeof(int));
	memset(idCountArray,0,idArraySize*sizeof(int));
	idCountArray[id] = 1;

	for (i=nstart+1;i<nstop;i++) {
		if (superCell.atoms[i].Znum != superCell.atoms[i-1].Znum) {
			/* check, if this Znum is already somewhere in the list, and if not,
			* add it to the list, extending its size, if necessary
			*/

			for (id=0;id<=idArrayPtr;id++) {
				if (superCell.atoms[i].Znum == idArray[id]) break;
			}
			if (id>idArrayPtr) {
				if (idArrayPtr >= idArraySize-1) {
					idArraySize *=2;
					idArray = (int *)realloc(idArray,idArraySize*sizeof(int));
					idCountArray = (int *)realloc(idCountArray,idArraySize*sizeof(int));
				}
				idArray[++idArrayPtr] = superCell.atoms[i].Znum;
			}
		}
		idCountArray[id] += 1;
	}

	/**************************************************************************************
	* Now that we found all the species present in the amorphous phase we can go ahead and 
	* make a list of them, together with the number of atoms present of each kind
	*/
	for (id=0;id<=idArrayPtr;id++) {
		if (idCountArray[id] > 0) { 
			if (idArray[id] >= MAX_MASS_INDEX) {
				printf("mass exceeds array!!! - extend array!\n");
				exit(0);
			}
			switch (idArray[id]) {
	  case 14: charge = 4.0; break;
	  case 7:  charge = -3.0; break;
	  case 8:  charge = -2.0; break;
	  case 39: charge = 3.0; break;
	  default: charge = 0.0; 
			}
			if ( fpAmorph != NULL ) {
				fprintf( fpAmorph, "%c%c %d\n", elTable[2 * idArray[id] - 2], elTable[2 * idArray[id] - 1], idCountArray[id] );
				fprintf( fpAmorph, "%d %7.3f %7.3f %7.3f %6.3f %6.3f %c%c\n", id + 1, 0.0, 0.0, 0.0,
					massArray[idArray[id]-1],charge,
					elTable[2*idArray[id]-2],elTable[2*idArray[id]-1]);    
			}
		}
	}
	if ( fpAmorph != NULL ) fprintf( fpAmorph, "end\n" );
	/***********************************************************************
	* The list of the potential parameters is next
	*/
	if ( fpAmorph != NULL ) {
		fprintf( fpAmorph, "buckingham\n" );
	for (id=0;id<=idArrayPtr;id++)
		for (j=id;j<=idArrayPtr;j++) {
			// The following parameters are from S. Garofalini, J. Am. Cer. Soc. 67, 133 (1987)
			switch (10*(id+1) + (j+1)) {
	  case 11: b = 0.0; break;  // N -N
	  case 12: b = 4.5; break;  // N -Si
	  case 13: b = 4.5; break;  // N -Y
	  case 14: b = 0.0; break;  // N -O
	  case 22: b = 1.8770; break;  // Si-Si
	  case 23: b = 4.4200; break;  // Si-Y
	  case 24: b = 2.9620; break;  // Si-O
	  case 33: b = 9.9706; break;  // Y -Y
	  case 34: b = 6.0802; break;  // Y -O
	  case 44: b = 0.7254; break;  // O -O		
	  default : b=0.0;
			}
			b = b*6.022e4;  // *1e-16*6.022e23*1e-3
			fprintf( fpAmorph, "%d %d %g %g %g\n", id + 1, j + 1, A, b, C );

		}
		fprintf( fpAmorph, "end\n" );

		/*************************************************************************
		* We can now write down the size of the MD cell:
		*/
		fprintf( fpAmorph, "%g %g %g 90 90 90 1 1 1\n", superCell.ax, superCell.by, superCell.cz );

		/****************************************************************************
		* now list the position of each one of the atoms in the amorphous phase
		* in FRACTIONAL COORDINATES !!!
		*/
		for(i=nstart;i<nstop;i++) {   
			x = superCell.atoms[i].x/superCell.ax;
			y = superCell.atoms[i].y/superCell.by;
			z = superCell.atoms[i].z/superCell.cz;
			if (fabs(x-1.0) < 1e-5) x = 0.0;
			if (fabs(y-1.0) < 1e-5) y = 0.0;
			if (fabs(z-1.0) < 1e-5) z = 0.0;
			fprintf( fpAmorph, "%c%c %7.5f %7.5f %7.5f\n", elTable[2 * superCell.atoms[i].Znum - 2],
				elTable[2*superCell.atoms[i].Znum-1],x,y,z);
		}
		if (nstart > 0)
			fprintf( fpAmorph, "frame1 %7.5f %7.5f %7.5f 1 0 0 0\n", superCell.cmx, superCell.cmy, superCell.cmz );
		fprintf( fpAmorph, "end\n" );
	}  // end of if fp != NULL
}


#define IA 16807 
#define IM 2147483647 
#define AM (1.0/IM) 
#define IQ 127773 
#define IR 2836 
#define MASK 123459876

/********************************************************************
* Minimal  random number generator of Park and Miller. 
* Returns a uniform random deviate between 0.0 and 1.0.
* Set idum to any integer, except MASK, to initialize the sequence.
*******************************************************************/
float ran(long *idum) { 
	static long k; 
	static float ans; 

	*idum ^= MASK; // XORing with MASK allows use of zero and other . 
	k=(*idum)/IQ;  // simple bit patterns for idum.
	*idum=IA*(*idum-k*IQ)-IR*k; // Compute idum=(IA*idum) % IM without over-  ows by Schrage s method. 
	if (*idum < 0) *idum += IM; 
	ans=AM*(*idum); // Convert idum to a  floating result. 
	*idum ^= MASK; // Unmask before return. 
	return ans; 
}


#define IA 16807 
#define IM 2147483647 
#define AM (1.0/IM) 
#define IQ 127773 
#define IR 2836 
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB) 
#define EPS 1.2e-7 
#define RNMX (1.0-EPS)
/**********************************************************************
*  Minimal  random number generator of Park and Miller with Bays-Durham shuffle 
* and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 
* (exclusive of the endpoint values). Call with idum a negative integer to initialize; 
* thereafter, do not alter idum between successive deviates in a sequence. 
* RNMX should approximate the largest  floating value that is less than 1.
*/
double ran1(long *idum) { 
	int j; 
	long k; 
	static long iy=0; 
	static long iv[NTAB]; 
	double temp; 
	if (*idum <= 0 || !iy) { // Initialize. 
		if (-(*idum) < 1) *idum=1; // Be sure to prevent  idum = 0. 
		else *idum = -(*idum); 
		for (j=NTAB+7;j>=0;j--) { // Load the shu e table (after 8 warm-ups). 
			k=(*idum)/IQ; 
			*idum=IA*(*idum-k*IQ)-IR*k; 
			if (*idum < 0) *idum += IM; 
			if (j < NTAB) iv[j] = *idum; 
		} 
		iy=iv[0]; 
	} 
	k=(*idum)/IQ; // Start here when not initializing. 
	*idum=IA*(*idum-k*IQ)-IR*k; // Compute idum=(IA*idum) % IM without overflows by Schrage s method. 
	if (*idum < 0) *idum += IM; 
	j=iy/NDIV; // Will be in the range 0..NTAB-1. 
	iy=iv[j];  // Output previously stored value and re ll the shu e table. 
	iv[j] = *idum; 
	if ((temp=AM*iy) > RNMX) return RNMX; // Because users don t expect endpoint values. 
	else return temp; 
}


/*****************************************************************
* Gaussian distribution with unit variance
* idum must be initailized to a negative integer 
****************************************************************/
double gasdev(long *idum) 
/* Returns a normally distributed deviate with zero mean and unit variance, 
* using ran1(idum) as the source of uniform deviates. */
{ 
	// float ran1(long *idum); 
	static int iset=0; 
	static float gset; 
	double fac,rsq,v1,v2; 
	if (*idum < 0) {
		iset=0; // Reinitialize. 
		//    printf("reinit gasdev\n");
	}
	if (iset == 0) { 
		/* We don t have an extra deviate handy, so 
		* pick two uniform numbers in the square extending from -1 to +1 in each direction, */
		do { 
			v1=2.0*ran1(idum)-1.0;  
			v2=2.0*ran1(idum)-1.0; 
			rsq=v1*v1+v2*v2;  // see if they are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again. 
		fac=sqrt(-2.0*log(rsq)/rsq); 
		/* Now make the Box-Muller transformation to get two normal deviates. 
		*  Return one and save the other for next time. 
		*/
		gset=v1*fac; 
		iset=1; // Set flag. 
		return v2*fac; 
	} 
	else { //We have an extra deviate handy,so unset the flag, and return it. 
		iset=0;  
		return gset; 
	} 
}

void writeSTEMinput(char* stemFile,char *cfgFile,MULS *muls) {
	FILE *fpSTEM;
	char folder[64];


	strcpy(folder,cfgFile);
	folder[strlen(folder)-4] = 0;  // cut off .cfg ending

	if ( (fpSTEM = fopen( stemFile, "w" )) == NULL ) {
		printf("writeSTEMinput: Error opening file %s\n",stemFile);
		exit(0);
	}

	fprintf( fpSTEM, "%% STEM input file for tomography series, cerated automatically by STEM\n\n" );

	fprintf( fpSTEM, "mode: CBED\n" );
	fprintf( fpSTEM, "print level: 1\nsave level: %d\n", muls->saveLevel );
	fprintf( fpSTEM, "filename: %s\n", cfgFile );
	fprintf( fpSTEM, "NCELLX: %d\nNCELLY: %d\nNCELLZ: %d/%d\n",
		muls->nCellX,muls->nCellY,muls->nCellZ,muls->cellDiv);
	fprintf( fpSTEM, "v0: %g\n", muls->v0 );
	fprintf( fpSTEM, "tds: %s\n", muls->tds ? "yes" : "no" );
	fprintf( fpSTEM, "temperature: %g\n", muls->tds_temp );
	fprintf( fpSTEM, "slice-thickness: %g\n", muls->sliceThickness );
	fprintf( fpSTEM, "periodicXY: no\n" );
	fprintf( fpSTEM, "periodicZ: no\n" );


	fprintf( fpSTEM, "nx: %d\n", muls->nx );
	fprintf( fpSTEM, "ny: %d\n", muls->ny );
	fprintf( fpSTEM, "Cs: %g\n", muls->Cs*1e-7 );
	fprintf( fpSTEM, "C5: %g\n", muls->C5 );
	fprintf( fpSTEM, "Cc: %g\n", muls->Cc*1e-7 );
	fprintf( fpSTEM, "dV/V: %g\n", muls->dE_E );
	fprintf( fpSTEM, "alpha: %g\n", muls->alpha );
	fprintf( fpSTEM, "AIS aperture: %g %% A\n", muls->aAIS );
	fprintf( fpSTEM, "smooth: %s\n", muls->ismoth ? "yes" : "no" );
	fprintf( fpSTEM, "defocus: %g \n", muls->df0 );
	fprintf( fpSTEM, "Source Size (diameter): %g \n", 2 * muls->sourceRadius );
	fprintf( fpSTEM, "gaussian: %s\n", muls->gaussFlag ? "yes" : "no" );
	fprintf( fpSTEM, "potential3D: %s \n", muls->potential3D ? "yes" : "no" );
	fprintf( fpSTEM, "atom radius: %g \n", muls->atomRadius );
	fprintf( fpSTEM, "plot V(r)*r: yes \n" );
	fprintf( fpSTEM, "bandlimit f_trans: yes\n" );
	fprintf( fpSTEM, "save potential: no	\n" );
	fprintf( fpSTEM, "one time integration: yes \n" );

	fprintf( fpSTEM, "tomo tilt: %lf\n", muls->tomoTilt );


	fprintf( fpSTEM, "Display Gamma: 0 \n" );
	fprintf( fpSTEM, "Folder: %s\n", folder );
	fprintf( fpSTEM, "Runs for averaging: %d\n", muls->avgRuns );
	fprintf( fpSTEM, "Structure Factors: DT  \n" );
	fprintf( fpSTEM, "show Probe: %s \n", muls->showProbe ? "yes" : "no" );
	fprintf( fpSTEM, "propagation progress interval: 10 \n" );
	fprintf( fpSTEM, "potential progress interval: 1000 \n" );
	fprintf( fpSTEM, "beams: n \n" );
	fprintf( fpSTEM, "sequence: 1 1\n" );

	fclose( fpSTEM );
}
