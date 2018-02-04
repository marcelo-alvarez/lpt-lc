#include "allvars.h"
#include "proto.h"
#include <math.h>

/////////////////////////////////////////////////////////////////////////
// THIS PROGRAM MAKES MAPS FROM A LINEAR DENSITY FIELD
//
//                                                AUTHOR: MARCELO ALVAREZ  
//                                             LAST EDIT:        09.08.17
/////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // Parse command line  
  CommandLine(argc, &argv[0]);

  // Read parameter file
  ReadParameterFile();
  int N=clParameters.N;
  int Nchunk = clParameters.Nchunk;
  Parameters.N = N;

  // Growth factor at initial time
  DInit=growth(clParameters.zInit,Parameters.Omegam,Parameters.Omegal,
	       Parameters.w);

  // Make FFTW plans
  plan  = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,N, N, N,
				  FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  iplan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,N, N, N,
				  FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
  
  rfftwnd_mpi_local_sizes(plan, &local_nx, &local_x_start,
			  &local_ny_after_transpose,
			  &local_y_start_after_transpose,
			  &total_local_size);

  ovrt=0;
  oram=0;
  GetOverhead(total_local_size,&ovrt,&oram);
  ReportMemory("after plan creation",total_local_size,ovrt,oram);

  // Allocate arrays
  AllocateArrays();
  ReportMemory("after array allocation",size_fftw,ovrt,oram);

  // Read input files
  if(clParameters.zmask==1)
    ReadGridFromFile( zmask, clParameters.RedshiftFile);
  if(clParameters.halomask==1)
    ReadHaloFile(clParameters.HaloFile);
  if(clParameters.deltain==1)
    ReadGridFromFile( delta1, clParameters.DeltaFile);
  MPI_Barrier(MPI_COMM_WORLD);

  // LPT Displacements
  if(clParameters.lptcode>1){
    Displace_2LPT(delta1, delta2, sx1, sy1, sz1, sx2, sy2, sz2);
    Displace_1LPT(delta2, sx2, sy2, sz2);
    if(clParameters.deltain==1)
      ReadGridFromFile( delta1, clParameters.DeltaFile);
  }
  if(clParameters.lptcode>0) Displace_1LPT(delta1, sx1, sy1, sz1);  
  MPI_Barrier(MPI_COMM_WORLD);

  double chunksize = floor(Parameters.Nnu/Nchunk);
  int Nnu_tot = Parameters.Nnu;
  float nu1 = Parameters.nu1;
  float nu2 = Parameters.nu2;
  float dnu = (nu2-nu1)/Nnu_tot;

  for (int ci=0; ci<Nchunk; ci++){
    // Make maps
       Parameters.nu1 = nu1 + ci*(chunksize)*dnu;
       Parameters.nu2 = Parameters.nu1 + (chunksize-1)*dnu;
       if (Parameters.nu2 >= nu2){
           Parameters.nu2 = nu2;
           Parameters.Nnu = (int)((nu2 - Parameters.nu1)/dnu);
       }
       else{
           Parameters.Nnu = chunksize;
       }
      if (myid==0){
       printf("\n Chunk %d \n",ci);
       printf("nu1=%f, nu2=%f, dnu=%f, nu0 = %f, numax=%f , Nnu=%d, chunksize=%f \n",Parameters.nu1,Parameters.nu2,dnu, nu1, nu2,Parameters.Nnu, chunksize);
      } // myid==0
    MakeMaps();
    MPI_Barrier(MPI_COMM_WORLD);

    // Write maps
    WriteMaps();
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // Finalize and return
  MPI_Finalize();  if(myid==0) printf("\n\n"); 
  
  fclose(stdout);
  exit(0);

}

