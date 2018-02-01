#define NUMMAPCODES 5

#define KAPCODE 0
#define KSZCODE 1
#define TAUCODE 2
#define DTBCODE 3
#define CIBCODE 4

#ifndef ALLVARS_H
#define ALLVARS_H

#include <mpi.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#ifdef NOTYPEPREFIX_FFTW
#include        <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	/* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif

#include "tables.h"
#include "globaltablevars.h"
#include "chealpix.h"

// IO
extern FILE *logfile;

// MPI Variables
extern int myid, nproc;

// FFTW Variables
extern rfftwnd_mpi_plan plan, iplan; 
extern int local_nx, local_x_start, 
           local_ny_after_transpose, local_y_start_after_transpose, 
           total_local_size;

// Slab sizes
extern long Nlocal;             // Local slab dimension 
extern long int size;           // Local slab size 
extern long int size_fftw;
extern long mapsize, tmapsize;

// Arrays
extern fftw_real *delta, *delta1, *delta2, *sx1, *sy1, *sz1, *sx2, *sy2, *sz2;
extern float *zmask;
extern float *halomask;
extern float *taumap,  *kszmap,  *dtbmap,  *kapmap,  *cibmap;
extern float *taumapl, *kszmapl, *dtbmapl, *kapmapl, *cibmapl;

// Halos
extern struct Halos{
  int   N; 
  float RTHmax, redshift;
  float *x,*y,*z,*vx,*vy,*vz,*xL,*yL,*zL,*RTH;
} halos;

// Parameters
extern struct Parameter{
  float Omegam, Omegab, Omegal, h, ns, Sigma8, w;
  float BoxSize, zInit, fov, theta, phi;
  float nu1,nu2;
  int Nnu;
  float InitialRedshift, FinalRedshift;  
  int N, NSide, NPixels, NRedshifts;
  char DeltaFile[256];
  int DoMap[NUMMAPCODES], CurrentCode;
} Parameters;

extern struct clParameter{
  int verbose, N, zmask, halomask, deltain, mapcode, lptcode, evolve, Nchunk;
  float BoxCenter[3], BoxSize, Periodicity, zInit, zKappa;
  char ParamFile[256], DeltaFile[256], RedshiftFile[256], 
    HaloFile[256], FluxFile[256], BaseOut[256];  
} clParameters;

// Redshift table
extern float dztable;

// Derived parameters
extern float DInit, ovrt, oram;

#endif

