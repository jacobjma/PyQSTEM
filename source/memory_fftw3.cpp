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

#include <stdlib.h>
#include <stdio.h>
#include "fftw3.h"
#include "memory_fftw3.h"
 
#ifndef WIN32
#include <stdint.h>
#endif
/*
#define PRINT_MESSAGE
*/




/*---------------------------- float1D() -------------------------------*/
/*
	1D array allocator for type float
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
float_tt *float1D( int n, const char *message )
{
	float_tt *m;
	
	m = (float_tt*) fftw_malloc( n * sizeof( float_tt) );
	if( m == NULL ) {
		printf("float1D() cannot allocate memory size=%d: %s\n",
		       n, message);
		exit( 0 );
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif

	return( m );

}  /* end float1D() */

/*---------------------------- double1D() -------------------------------*/
/*
	1D array allocator for type double
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
double* double1D( int n, const char *message )
{
	double *m;
	
	m = (double*) fftw_malloc( n * sizeof( double ) );
	if( m == NULL ) {
		printf("double1D() cannot allocate memory size=%d: %s\n",
		       n, message);
		exit( 0 );
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif

	return( m );

} /* end double1D() */


/*---------------------------- short2D() -------------------------------*/
/*
	2D array allocator for type short
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
short **short2D( int nx, int ny, const char *message )
{	short **m;
	int i;

	m = (short**) fftw_malloc( nx * sizeof( short* ) ); 
	if( m == NULL ) {
		printf("short2D cannot allocate pointers, size=%d : %s\n",
		       nx, message );
		exit(0);
	}
	m[0] = (short *) fftw_malloc( ny *nx* sizeof(short) );
	if( m[0] == NULL ){
	  printf("long2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = &(m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif

	return m;
}  /* end short2d() */

/*---------------------------- int2D() -------------------------------*/
/*
	2D array allocator for type int
	make space for m[0...(nx-1)][0..(ny-1)]

*/
int **int2D( int nx, int ny, const char *message )
{	int **m;
	int i;

	m = (int**) fftw_malloc( nx * sizeof(int* ) ); 
	if( m == NULL ) {
		printf("int2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (int *) fftw_malloc( ny *nx* sizeof(int) );
	if( m[0] == NULL ){
	  printf("int2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = (int *)(&m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s (int) = %d\n",message,(int)m);
#endif

	return m;

}  /* end int2D() */

/*---------------------------- long2D() -------------------------------*/
/*
	2D array allocator for type long
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
long **long2D( int nx, int ny, const char *message )
{	long **m;
	int i;

	m = (long**) fftw_malloc( nx * sizeof( long* ) ); 
	if( m == NULL ) {
		printf("long2D cannot allocate pointers, size=%d : %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (long *) fftw_malloc( ny *nx* sizeof(long) );
	if( m[0] == NULL ){
	  printf("long2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = &(m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif

	return m;

}  /* end long2d() */

/*---------------------------- float32_2D() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/

float **float32_2D( int nx, int ny, const char *message )
{	
	float **m;
	int i;

	m = (float**) fftw_malloc( nx * sizeof( float* ) ); 
	if( m == NULL ) {
		printf("float2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (float *) fftw_malloc( ny *nx* sizeof( float ) );
	if( m[0] == NULL ){
	  printf("float2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = (float *)(&m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s (float_tt) = %d\n",message,(int)m);
#endif

	return m;

}  /* end float2D() */


/*---------------------------- float2D() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/
float_tt **float2D( int nx, int ny, const char *message )
{	
	float_tt **m;
	int i;

	m = (float_tt**) fftw_malloc( nx * sizeof( float_tt* ) ); 
	if( m == NULL ) {
		printf("float2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (float_tt *) fftw_malloc( ny *nx* sizeof( float_tt ) );
	if( m[0] == NULL ){
	  printf("float2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = (float_tt *)(&m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s (float_tt) = %d\n",message,(int)m);
#endif

	// initialize array to 0
	for (int ix=0;ix<nx;ix++) for (int iy=0;iy < ny; iy++)
	{
		m[ix][iy] = 0.0f;
	}
	//memset(m, 0, sizeof(float_tt)*ny*nx);

	return m;

}  /* end float2D() */
/*---------------------------- float3D() -------------------------------*/
/*
	3D array allocator for type float_tt
	make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

*/
float_tt ***float3D( int nx, int ny,int nz, const char *message )
{	
  float_tt ***m;
  int i,j;
  
  m = (float_tt***) fftw_malloc( nx * sizeof(float_tt**) ); 
  if( m == NULL ) {
    printf("float3D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (float_tt**)fftw_malloc(ny*sizeof(float_tt*));
    if (m[i] == NULL) {
      printf("float3D cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = (float_tt*) fftw_malloc(nz*ny*nx* sizeof(float_tt) );
  if( m[0] == NULL ){
    printf("float2D cannot allocate arrays, size=%d: %s\n",
	   ny*nx, message );
    exit(0);
  }
  for (i=0; i<nx; i++) for (j=0;j<ny;j++) {
    m[i][j] = (float_tt*)(&(m[0][0][nz*(i*ny+j)]));
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (float_tt) = %d\n",message,(int)m);
#endif
  
  return m;
    
}  /* end float3D() */


/*---------------------------- float32_3D() -------------------------------*/
/*
	3D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

*/
float ***float32_3D( int nx, int ny,int nz, const char *message )
{	
  float ***m;
  int i,j;
  
  m = (float***) fftw_malloc( nx * sizeof(float**) ); 
  if( m == NULL ) {
    printf("float3D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (float **)fftw_malloc(ny*sizeof(float*));
    if (m[i] == NULL) {
      printf("float3D cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = (float*) fftw_malloc(nz*ny*nx* sizeof(float) );
  if( m[0] == NULL ){
    printf("float32_3D cannot allocate arrays, size=%d: %s\n",
	   ny*nx, message );
    exit(0);
  }
  for (i=0; i<nx; i++) for (j=0;j<ny;j++) {
    m[i][j] = (float*)(&(m[0][0][nz*(i*ny+j)]));
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (float32) = %d\n",message,(int)m);
#endif
  
  return m;
    
}  /* end float32_3D() */




/*---------------------------- double2D() -------------------------------*/
/*
	2D array allocator for type doubel
	make space for m[0...(nx-1)][0..(ny-1)]

*/
double **double2D( int nx, int ny, const char *message )
{	double **m;
	int i;

	m = (double**) fftw_malloc( nx * sizeof(double* ) ); 
	if( m == NULL ) {
		printf("double2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (double *) fftw_malloc( ny *nx* sizeof(double) );
	if( m[0] == NULL ){
	  printf("double2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = &(m[0][i*ny]);
	}

	// initialize array to 0
	for (int ix=0;ix<nx;ix++) for (int iy=0;iy < ny; iy++)
	{
		m[ix][iy] = 0.0;
	}
	//memset(m, 0, sizeof(float_tt)*ny*nx);

#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif

	return m;

}  /* end double2D() */

/*---------------------------- complex2D() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/
fftw_complex **complex2D( int nx, int ny, const char *message)
{
  fftw_complex **m;
  int i;
  
  m = (fftw_complex**) fftw_malloc( nx * sizeof(fftw_complex*) ); 
  if( m == NULL ) {
    printf("float2D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  
  m[0] = (fftw_complex*) fftw_malloc( ny *nx* sizeof(fftw_complex) );
  if( m[0] == NULL ){
    printf("float2D cannot allocate arrays, size=%d: %s\n",
	   ny*nx, message );
    exit(0);
  }
  for (i=1; i<nx; i++){
    m[i] = (fftw_complex*)(&m[0][i*ny]);
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif
  
  return m;
  
}  /* end complex2D() */


/*---------------------------- complex2Df() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/
fftwf_complex **complex2Df( int nx, int ny, const char *message)
{	fftwf_complex **m;
	int i;

	m = (fftwf_complex**) fftw_malloc( nx * sizeof(fftwf_complex*) ); 
	if( m == NULL ) {
		printf("float2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (fftwf_complex*) fftw_malloc( ny *nx* sizeof(fftwf_complex) );
	if( m[0] == NULL ){
	  printf("float2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = (fftwf_complex*)(&m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif

	return m;

}  /* end complex2Df() */

/*---------------------------- complex3D() -------------------------------*/
/*
	3D array allocator for type fftw_complex
	make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

*/
fftw_complex ***complex3D( int nx, int ny,int nz, const char *message)
{	
  fftw_complex ***m;
  int i,j;
  
  m = (fftw_complex***)fftw_malloc( nx * sizeof(fftw_complex**) ); 
  if( m == NULL ) {
    printf("complex3D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (fftw_complex**)fftw_malloc(ny*sizeof(fftw_complex*));
    if (m[i] == NULL) {
      printf("complex3D cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = (fftw_complex*) fftw_malloc(nz*ny*nx* sizeof(fftw_complex) );
  if( m[0] == NULL ){
    printf("float2D cannot allocate arrays, size=%d: %s\n",
	   ny*nx, message );
    exit(0);
  }
  for (i=0; i<nx; i++) for (j=0;j<ny;j++) {
    m[i][j] = (fftw_complex*)(&(m[0][0][nz*(i*ny+j)]));
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif
  
  return m;
    
}  /* end complex3D() */

/*---------------------------- complex3Df() -------------------------------*/
/*
	3D array allocator for type fftwf_complex
	make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

*/
fftwf_complex ***complex3Df( int nx, int ny,int nz, const char *message)
{	
  fftwf_complex ***m;
  int i,j;
  
  m = (fftwf_complex***)fftwf_malloc( nx * sizeof(fftwf_complex**) ); 
  if( m == NULL ) {
    printf("complex3Df cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (fftwf_complex**)fftwf_malloc(ny*sizeof(fftwf_complex*));
    if (m[i] == NULL) {
      printf("complex3Df cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = (fftwf_complex*) fftwf_malloc(nz*ny*nx* sizeof(fftwf_complex) );
  if( m[0][0] == NULL ){
    printf("complex3Df cannot allocate consecutive memory %d MB (for array %s)\n",
	    nz*ny*nx* sizeof(fftwf_complex)/(1024*1024),message);
    exit(0);
  }
  for (i=0; i<nx; i++) for (j=0;j<ny;j++) {
    m[i][j] = (fftwf_complex*)(&(m[0][0][nz*(i*ny+j)]));
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif
  
  return m;
    
}  /* end complex3D() */



/*---------------------------- any2D() -------------------------------*/
/*
	2D array allocator for any type of size 'size'
	make space for m[0...(nx-1)][0..(ny-1)]

*/
void **any2D( int nx, int ny,int size, const char *message )
{	void **m;
	int i;

	m = (void **)fftw_malloc( nx * sizeof(void *)); 
	if( m == NULL ) {
		printf("any2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] =  fftw_malloc( ny *nx* size );
	if( m[0] == NULL ){
	  printf("any2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = (void *)((intptr_t)(m[0])+(i*ny*size));
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif

	return m;

}  /* end any2D() */


/*---------------------------- any3D() -------------------------------*/
/*
	3D array allocator for any type 
	make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

*/
void ***any3D( int nx, int ny,int nz,int size, const char *message )
{	
  void ***m;
  int i,j;
  
  m = (void***) fftw_malloc( nx * sizeof(void **)); 
  if( m == NULL ) {
    printf("any3D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (void **)fftw_malloc(ny*sizeof(void *));
    if (m[i] == NULL) {
      printf("any3D cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = fftw_malloc(nz*ny*nx* size );
  if( m[0][0] == NULL ){
	printf("any3Df cannot allocate consecutive memory %d MB (for array %s)\n",
	    nz*ny*nx*size/(1024*1024),message);
    exit(0);
  }
  for (i=0; i<nx; i++) for (j=0;j<ny;j++) {
    m[i][j] =(void *)((intptr_t)(m[0][0])+size*nz*(i*ny+j));
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif
  
  return m;
    
}  /* end any3D() */


