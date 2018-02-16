#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "proto.h"
#include "allvars.h"
#include "chealpix.h"

#define FLOAT 1
#define DOUBLE 2
#define STRING 3
#define INT 4
#define MAXTAGS 300

void WriteSingleMap(float *, char *);
void WriteSingletMap(float *, char *);

void ReadParameterFile()
{

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  char fname[256];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];

  sprintf(fname,clParameters.ParamFile);

  nt=0;
  
  strcpy(tag[nt], "Omegam");
  addr[nt] = &Parameters.Omegam;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Omegab");
  addr[nt] = &Parameters.Omegab;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Omegal");
  addr[nt] = &Parameters.Omegal;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "h");
  addr[nt] = &Parameters.h;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "ns");
  addr[nt] = &Parameters.ns;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Sigma8");
  addr[nt] = &Parameters.Sigma8;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "w");
  addr[nt] = &Parameters.w;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "BoxSize");
  addr[nt] = &Parameters.BoxSize;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "zInit");
  addr[nt] = &Parameters.zInit;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "DeltaFile");
  addr[nt] = &Parameters.DeltaFile;
  id[nt++] = STRING;  

  strcpy(tag[nt], "fov");
  addr[nt] = &Parameters.fov;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "theta");
  addr[nt] = &Parameters.theta;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "phi");
  addr[nt] = &Parameters.phi;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "nu1");
  addr[nt] = &Parameters.nu1;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "nu2");
  addr[nt] = &Parameters.nu2;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Nnu");
  addr[nt] = &Parameters.Nnu;
  id[nt++] = INT;

  strcpy(tag[nt], "NPixels");
  addr[nt] = &Parameters.NPixels;
  id[nt++] = INT;

  strcpy(tag[nt], "NSide");
  addr[nt] = &Parameters.NSide;
  id[nt++] = INT;  

  strcpy(tag[nt], "NRedshifts");
  addr[nt] = &Parameters.NRedshifts;
  id[nt++] = INT;

  strcpy(tag[nt], "InitialRedshift");
  addr[nt] = &Parameters.InitialRedshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "FinalRedshift");
  addr[nt] = &Parameters.FinalRedshift;
  id[nt++] = FLOAT;

  if((fd = fopen(fname, "r")))
    {
      sprintf(buf, "%s%s", fname, "-usedvalues");
      if(!(fdout = fopen(buf, "w")))
	{
	  printf("error opening file '%s' \n", buf);
	}
      else
	{
	  while(!feof(fd))
	    {
	      *buf = 0;
	      fgets(buf, 200, fd);
	      if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		continue;
	      
	      if(buf1[0] == '%' || buf1[0] == '#')
		continue;
	      
	      for(i = 0, j = -1; i < nt; i++)
		if(strcmp(buf1, tag[i]) == 0)
		  {
		    j = i;
		    tag[i][0] = 0;
		    break;
		  }
	      
	      if(j >= 0)
		{
		  switch (id[j])
		    {
		    case FLOAT:
		      *((float *) addr[j]) = atof(buf2);
		      fprintf(fdout, "%-35s%f\n", buf1, *((float *) addr[j]));
		      break;
		    case DOUBLE:
		      *((float *) addr[j]) = atof(buf2);
		      fprintf(fdout, "%-35s%f\n", buf1, *((float *) addr[j]));
		      break;
		    case STRING:
		      strcpy((char *)addr[j], buf2);
		      fprintf(fdout, "%-35s%s\n", buf1, buf2);
		      break;
		    case INT:
		      *((int *) addr[j]) = atoi(buf2);
		      fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
		      break;
		    }
		}
	      else
		{
		  fprintf(stdout, 
			  "Error in %s: Tag '%s' not allowed or multiple defined.\n",
			  fname, buf1);
		}
	    }
	  fclose(fd);
	  fclose(fdout);
	  
	}
    }
  else
    {
      printf("\nParameter file %s not found.\n\n", fname);
    }
  
   if(myid==0) printf("\n Parameter file read...");

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
  
}

void ReadHaloFile(char *fname){

  FILE *fd = fopen(fname,"rb");
  if(fd==NULL){
    printf("could not open halo file %s\n",fname);
    exit(0);
    MPI_Finalize();
  }
  
  fread(&halos.N, 4,1,fd);
  fread(&halos.RTHmax, 4,1,fd);
  fread(&halos.redshift, 4,1,fd);

  int Nhalos = halos.N;

  if(myid==0){
    printf("\n\n");
    printf(" Nhalos = %d\n",Nhalos);
    printf(" RTHmax = %f\n",halos.RTHmax);
  }

  halos.x   = new float[Nhalos];
  halos.y   = new float[Nhalos];
  halos.z   = new float[Nhalos];
  halos.vx  = new float[Nhalos];
  halos.vy  = new float[Nhalos];
  halos.vz  = new float[Nhalos];
  halos.RTH = new float[Nhalos];
  halos.xL  = new float[Nhalos];
  halos.yL  = new float[Nhalos];
  halos.zL  = new float[Nhalos];

  for(long i=0;i<Nhalos;i++){
    fread(&halos.x[i],  4,1,fd);
    fread(&halos.y[i],  4,1,fd);
    fread(&halos.z[i],  4,1,fd);
    fread(&halos.vx[i], 4,1,fd);
    fread(&halos.vy[i], 4,1,fd);
    fread(&halos.vz[i], 4,1,fd);
    fread(&halos.RTH[i],4,1,fd);
    fread(&halos.zL[i], 4,1,fd);
    fread(&halos.yL[i], 4,1,fd);
    fread(&halos.xL[i], 4,1,fd);
  }

  if(myid==0) printf("\n halo data read...");

  /*
  halos.N = 1;
  halos.xL[0] = 0;
  halos.yL[0] = 0;
  halos.zL[0] = 4000;
  halos.RTH[0] = 500;
  */
  fclose(fd);

}

void ReadGridFromFile(float *grid, char *fname)
{
    
  int N=clParameters.N;

  float *slab = new float[N*N];

  long offset = (long)myid*(long)N*(long)N*(long)Nlocal*(long)sizeof(float);
  long bsize = N*N*sizeof(float);

  FILE *fd = fopen(fname,"rb");
  if(fd==NULL){
    printf("could not open grid file %s\n",fname);
    exit(0);
    MPI_Finalize();
  }

  for(long i=0;i<Nlocal;i++){
    // Read this slab from the input file

    long int offset_local = i*N*N*sizeof(float) + offset;
    parallel_read(fname, bsize, offset_local, slab);

    for(long j=0;j<N;j++){
      for(long k=0;k<N;k++){
	long index = i*N*(N+2)+j*(N+2)+k;
	grid[index] = slab[j*N+k]; 
      }
    }

  }

  delete[] slab;

  if(myid==0) printf("\n grid data read...");

}

void WriteLPT()
{
  
  if(myid==0) printf("\n 0 myid, Nlocal = %d %d\n",myid,Nlocal);

  int N=clParameters.N;
  float BoxSize=Parameters.BoxSize;

  if(myid==0) printf("\n BoxSize, Nlocal = %f %d\n",BoxSize,Nlocal);

  float *slab = new float[N*N];

  long blocksize = (long)N * (long)N * (long)N * (long)(sizeof(float));
  long offset_inblock = (long)myid*(long)N*(long)N*(long)Nlocal*(long)(sizeof(float));
  long bsize = N*N*sizeof(float);
  long offset;

  char fnamelpt[256];

  // Growth factor at z=0 relative to zInit
  float D0  = growth((float)0.0,Parameters.Omegam,
		    Parameters.Omegal, Parameters.w);
  float D01 = D0 / growth(clParameters.zInit,Parameters.Omegam,
	   Parameters.Omegal, Parameters.w);

  printf("D0, D01, clParameters.zInit = %f %f %f\n",D0,D01,clParameters.zInit);
  float D02 = 3./7. * D01 * D01;

  sprintf(fnamelpt,"%s_lpt.bin",clParameters.BaseOut);
  printf("\n myid, Nlocal = %d %d\n",myid,Nlocal);

  offset = 0*blocksize + offset_inblock;
  for(long i=0;i<Nlocal;i++){for(long j=0;j<N;j++){for(long k=0;k<N;k++){
	long index = i*N*(N+2)+j*(N+2)+k; slab[j*N+k] = sx1[index]*D01;}}
    parallel_write(fnamelpt,bsize, offset, slab);offset+=bsize;}
  printf("myid, offset = %d %ld\n",myid,offset);

  offset = 1*blocksize + offset_inblock;
  for(long i=0;i<Nlocal;i++){for(long j=0;j<N;j++){for(long k=0;k<N;k++){
  long index = i*N*(N+2)+j*(N+2)+k; slab[j*N+k] = sy1[index]*D01;}}
    parallel_write(fnamelpt,bsize, offset, slab);offset+=bsize;}
  printf("myid, offset = %d %ld\n",myid,offset);

  offset = 2*blocksize + offset_inblock;
  for(long i=0;i<Nlocal;i++){for(long j=0;j<N;j++){for(long k=0;k<N;k++){
  long index = i*N*(N+2)+j*(N+2)+k; slab[j*N+k] = sz1[index]*D01;}}
    parallel_write(fnamelpt,bsize, offset, slab);offset+=bsize;}
  printf("myid, offset = %d %ld\n",myid,offset);

  offset = 3*blocksize + offset_inblock;
  for(long i=0;i<Nlocal;i++){for(long j=0;j<N;j++){for(long k=0;k<N;k++){
  long index = i*N*(N+2)+j*(N+2)+k; slab[j*N+k] = sx2[index]*D02;}}
    parallel_write(fnamelpt,bsize, offset, slab);offset+=bsize;}
  printf("myid, offset = %d %ld\n",myid,offset);

  offset = 4*blocksize + offset_inblock;
  for(long i=0;i<Nlocal;i++){for(long j=0;j<N;j++){for(long k=0;k<N;k++){
  long index = i*N*(N+2)+j*(N+2)+k; slab[j*N+k] = sy2[index]*D02;}}
    parallel_write(fnamelpt,bsize, offset, slab);offset+=bsize;}
  printf("myid, offset = %d %ld\n",myid,offset);

  offset = 5*blocksize + offset_inblock;
  for(long i=0;i<Nlocal;i++){for(long j=0;j<N;j++){for(long k=0;k<N;k++){
  long index = i*N*(N+2)+j*(N+2)+k; slab[j*N+k] = sz2[index]*D02;}}
    parallel_write(fnamelpt,bsize, offset, slab);offset+=bsize;}
  printf("myid, offset = %d %ld\n",myid,offset);

  delete[] slab;

  if(myid==0) printf("\n 2LPT data written...");
}

void WriteMaps()
{

  if(myid==0) {
    printf("\n Writing maps...\n");
  } else {
    MPI_Barrier(MPI_COMM_WORLD);
    return;
  }

  if(myid==0){
    
    if(Parameters.DoMap[KAPCODE]) WriteSingleMap( kapmap,"kap");
    if(Parameters.DoMap[KSZCODE]) WriteSingleMap (kszmap,"ksz");
    if(Parameters.DoMap[TAUCODE]) WriteSingleMap (taumap,"tau");
    if(Parameters.DoMap[CIBCODE]) WriteSingleMap (cibmap,"cib");
    if(Parameters.DoMap[DTBCODE]) WriteSingletMap(dtbmap,"dtb");

    if(myid==0) printf("\n Maps written...");

  }

  MPI_Barrier(MPI_COMM_WORLD);
  return;

}

void WriteSingleMap(float *map, char *base){

  char fname[256], coord[1];
  FILE *fout;
  int  mapsize = Parameters.NSide*Parameters.NSide*12;

  // binary format
  sprintf(fname,"%s_%s.bin",clParameters.BaseOut,base);
  fout = fopen(fname,"wb");   
  fwrite(map,4,mapsize,fout);
  fclose(fout);

  // fits format
  sprintf(fname,"!%s_%s.fits",clParameters.BaseOut,base);          
  sprintf(coord,"C");
  if(clParameters.binary_only==0) write_healpix_map(map, Parameters.NSide, fname, 1, coord);

}

void WriteSingletMap(float *map, char *base){

  char fname[256], coord[1];
  FILE *fout;
  int  mapsize = Parameters.NSide*Parameters.NSide*12;

  float nu, dnu, nu1, nu2;
  int   Nnu;

  nu1 = Parameters.nu1; nu2 = Parameters.nu2; Nnu = Parameters.Nnu; 
  dnu = (nu2 - nu1) / Nnu ;
  for(int inu=1;inu<Nnu-1;inu++){
    // Exclude the first and last frequencies. Ensure chunks overlap.
    nu = nu1 + (inu+0.5) * dnu ;
    
    // binary format
    sprintf(fname,"%s_%s_%6.2f.bin",clParameters.BaseOut,base,nu);
    fout = fopen(fname,"wb");   
    fwrite(&map[inu*mapsize],4,mapsize,fout);
    fclose(fout);
    
    // fits format
    sprintf(fname,"!%s_%s_%6.2f.fits",clParameters.BaseOut,base,nu);
    sprintf(coord,"C");
    if(clParameters.binary_only==0) write_healpix_map(&map[inu*mapsize], Parameters.NSide, fname, 1, coord);

  }
}
