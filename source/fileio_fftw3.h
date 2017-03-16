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

#ifndef FILEIO_H
#define FILEIO_H

#include "data_containers.h"
#include "stemtypes_fftw3.h"

atom *readUnitCell(int *natom,char *fileName,MULS *muls,int handleVacancies);
void replicateUnitCell(int ncoord,int *natom,MULS *muls,atom* atoms,int handleVacancies);
atom *tiltBoxed(int ncoord,int *natom, MULS *muls,atom *atoms,int handleVacancies);
int writePDB(atom *atoms,int natoms,char *fileName,MULS *muls);
int writeCFG(atom *atoms,int natoms,char *fileName,MULS *muls);
// write CFG file using atomic positions stored in pos, Z's in Znum and DW-factors in dw
// the unit cell is assumed to be cubic
int writeCFGFractCubic(double *pos,int *Znum,double *dw,int natoms,char *fileName,
		       double a,double b,double c);
int readCubicCFG(double **pos,double **dw, int **Znum, double *ax,double *by,double *cz,
		 double ctiltx, double ctilty);

void writeSTEMinput(char* stemFile,char *cfgFile,MULS *muls);

/* Helper functions for above functions: */
int ReadLine( FILE* fpRead, char* cRead, int cMax, const char *mesg );
int getZNumber(char *element);
int readCFGCellParams(MULS *muls, double **Mm, char *fileName);
int readCSSRCellParams(MULS *muls, double **Mm, char *fileName);

void writeFrameWork(FILE *fp,superCellBox superCell);
void writeAmorphous(FILE *fp,superCellBox superCell,int nstart,int nstop);

double gasdev(long *idum); 
double ran1(long *idum);
float ran(long *idum);
int atomCompareZYX(const void *atPtr1,const void *atPtr2);
int atomCompareZnum(const void *atPtr1,const void *atPtr2);
#endif /* FILEIO_H */
