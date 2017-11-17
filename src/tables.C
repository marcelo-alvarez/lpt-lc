#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tables.h"
#include <mpi.h>

int   RedshiftLnM_nz, RedshiftLnM_nm;
float RedshiftLnM_zinitial, RedshiftLnM_zfinal; 
float RedshiftLnM_minitial, RedshiftLnM_mfinal; 
float OmegaHI(float, float, float);

//////////////// [Redshift, ln M] to float //////////////////

float RedshiftLnM2Float(float redshift, float lnM, double *table)
{

  float dztable = ((float)RedshiftLnM_zfinal - RedshiftLnM_zinitial) / RedshiftLnM_nz;
  float dmtable = ((float)RedshiftLnM_zfinal - RedshiftLnM_minitial) / RedshiftLnM_nm;

  int zbin = (redshift - RedshiftLnM_zinitial) / dztable ;
  int mbin = (lnM      - RedshiftLnM_minitial) / dmtable ;

  float fz1 = redshift - (RedshiftLnM_zinitial + zbin*dztable) ;
  float fm1 = lnM      - (RedshiftLnM_minitial + mbin*dmtable) ;

  if(zbin<0 || zbin>RedshiftLnM_nz-1 || mbin<0 || mbin>RedshiftLnM_nm-1 ){
    printf("Error: z, m, zmin, zmaz, mmin, mmax = %f %f %f %f %f %f\n",
	   redshift,lnM,RedshiftLnM_zinitial,RedshiftLnM_zfinal,
	                RedshiftLnM_minitial,RedshiftLnM_mfinal);
    MPI_Finalize();
    exit(0);
  }

  int bin00 = (zbin+0) * RedshiftLnM_nm + (mbin+0);
  int bin01 = (zbin+0) * RedshiftLnM_nm + (mbin+1);
  int bin10 = (zbin+1) * RedshiftLnM_nm + (mbin+0);
  int bin11 = (zbin+1) * RedshiftLnM_nm + (mbin+1);

  float fz0 = 1-fz1;
  float fm0 = 1-fm1;

  float value = fz0 * fm0 * table[bin00] + 
                fz0 * fm1 * table[bin01] + 
                fz1 * fm0 * table[bin10] + 
                fz1 * fm1 * table[bin11] ;

  return value;

}

//////////////// Redshift to float //////////////////

float Float2Float(float var, int nvart, float varmint, float varmaxt, double **table,
		  float defval)
{

  double *tablel = *table;
  float dvart = (varmaxt - varmint) / (nvart-1);


  int bin = (var - varmint) / dvart ; 
  float f = (var - (varmint + bin*dvart) ) / dvart ; 

  if(var<varmint || var>varmaxt) return defval;

  return (1-f)*tablel[bin] + f*tablel[bin+1];

}

//////////////// Redshift to float //////////////////

float Redshift2Float(float redshift, double *table)
{

  float dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;


  int bin = (redshift - ZTABLE_INITIAL) / dztable ;
  float f = redshift - (ZTABLE_INITIAL + bin*dztable) ;

  if(bin<0 || bin>NZTABLE-1){
    printf("Error: z, zmin, zmaz = %f %f %f\n",redshift,ZTABLE_INITIAL,ZTABLE_FINAL);
    MPI_Finalize();
    exit(0);
  }
  float value = (1-f)*table[bin] + f*table[bin+1];

  return value;

}

//////////////// Radius to float //////////////////

float Radius2Float(float radius, double *table)
{

  float drtable = ((float)RTABLE_FINAL - RTABLE_INITIAL) / NRTABLE;

  int bin = (radius - RTABLE_INITIAL) / drtable ;
  float f = radius - (RTABLE_INITIAL + bin*drtable) ;
  float value = (1-f)*table[bin] + f*table[bin+1];

  return value;

}

//////////////// Wavenumber to float //////////////////

float Wavenumber2Float(int nk, float k, double *table, double *wntable){

  for(int i=0;i<nk-1;i++){
    if(k>=wntable[i] && k<=wntable[i+1]){
      float f = (log(k) - log(wntable[i]))/(log(wntable[i+1])-log(wntable[i]));
      float p = (1-f)*table[i] + f*table[i+1];
      return p;
    }
  }

  printf("\n Error: entry not found in wavenumber table for k = %e\n",k);
  for(int i=0;i<nk-1;i++) printf("%d %e %e\n",i,wntable[i],table[i]);

  MPI_Finalize();
  exit(0);

}

//////////////// Wavenumber to float (uniform log spacing in k) //////////////////

float Wavenumber2FloatLogSpace(int nk, float kmin, float kmax, float k, double *table){


  k    = log(k);
  kmin = log(kmin);
  kmax = log(kmax);

  float dk = (kmax-kmin)/nk;


  int bin = (k-kmin) / dk ;
  float f = k - (kmin + bin*dk) ;
  float value = (1-f)*table[bin] + f*table[bin+1];

  return value;

}

//////////////// Radius table //////////////////

void SetRedshift2RadiusTable(float h, float Omegam, float Omegal, double *table){

  double H0 = h*100.; // in km/sec/Mpc
  double c  = 2.99792458e5;      // in km/sec
  float r0  = c/H0;              // in Mpc 

  float dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;
  table[0]=0;

  for(int bin=1;bin<NZTABLE;bin++){
    double z  = ZTABLE_INITIAL + (bin-0.5)*dztable;
    double dr = r0 * dztable / sqrt(Omegam*pow((1+z),3)+Omegal);
    table[bin]=table[bin-1]+dr;
  }

}

//////////////// Redshift table //////////////////

void SetRadius2RedshiftTable(float h, float Omegam, float Omegal, double *table){

  double H0 = h*100.; // in km/sec/Mpc
  double c  = 2.99792458e5;      // in km/sec
  float r0  = c/H0;              // in Mpc 

  float drtable = ((float)RTABLE_FINAL - RTABLE_INITIAL) / NRTABLE;
  table[0]=0;

  for(int bin=1;bin<NRTABLE;bin++){
    double r  = RTABLE_INITIAL + (bin-0.5)*drtable;
    double z  = table[bin-1];
    double dz = drtable * sqrt(Omegam*pow((1+z),3)+Omegal) / r0;
    table[bin]=table[bin-1]+dz;
  }

}

//////////////// History table //////////////////

void SetRedshift2HistoryTable(int N, int Nlocal, float *zreion, double *table){

  float dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;

  double *history  = new double[NZTABLE]();

  for(int i=0;i<NZTABLE;i++) history[i]=0;

  for(int i=0;i<Nlocal;i++){
    for(int j=0;j<N;j++){
      for(int k=0;k<N;k++){
	
	int index=i*(N+2)*N+j*(N+2)+k;
	float zcur=zreion[index];
	int bin=(int)((zcur-ZTABLE_INITIAL)/dztable);
	if(bin>=0 && bin<NZTABLE) history[bin]++;

      }
    }
  }

  MPI_Allreduce(history, table, NZTABLE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  

  for(int i=0;i<NZTABLE;i++) table[i]/=pow(N,3);
  for(int i=NZTABLE-2;i>=0;i--) table[i]+=table[i+1];

}

//////////////// Dg table //////////////////

void SetRedshift2DgTable(double *history, double *table){

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  float dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;

  for(int i=1;i<NZTABLE-1;i++){
    float x1    = history[i-1];
    float x2    = history[i+1];
    float x     = history[i];
    double dxdz = (x2-x1)/2/dztable;
    float z     = ZTABLE_INITIAL + i*dztable; 
    table[i] = 1.5*pow(1+z,0.5)*x + pow(1+z,1.5)*dxdz;
  }

  table[0]=table[1];
  table[NZTABLE-1]=table[NZTABLE-2];

  return;

}

//////////////// Tau table //////////////////

void SetRedshift2TauTable(float h, float Omegab, float Omegam, double *table, double *Redshift2HistoryTable){

  float dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;

  float thompson = 6.65e-25;
  float YHe = 0.25;
  float hubble0 = 100./3.086e19*h;
  float ne0=Omegab*h*h*1.88e-29/1.67e-24;
  float c=3e10;
  
  float dtau0=thompson*ne0*c/hubble0;
  
  float NHe;
  for(int i=1;i<NZTABLE-1;i++){
    
    float z = ZTABLE_INITIAL + (i+0.5)*dztable;
    float x = Redshift2HistoryTable[i];

    NHe=1;
    if(z<3) NHe=2;

    float dtau = x*dtau0*pow((1+z),2)/sqrt(Omegam*pow((1+z),3)+1-Omegam)*
           (1-YHe+NHe*YHe/4)*dztable;

    table[i]+=table[i-1]+dtau;

  }

}

//////////////// Lensing kappa table //////////////////

void SetRedshift2WKappaTable(float h, float Omegam, float Omegal, float SourceRedshift, double *table, double *Redshift2RadiusTable){

  // Get distance to LSS at z=1100
  double H0 = h*100.; // in km/sec/Mpc
  double  c = 2.99792458e5;      // in km/sec
  float  r0 = c/H0;              // in Mpc 
  int    nz = 100000; 
 
  float  dztable = SourceRedshift / nz;
  float  chi0    = 0.;

  for(int bin=1;bin<nz;bin++){
    double z  = (bin-0.5)*dztable;
    double dr = r0 * dztable / sqrt(Omegam*pow((1+z),3)+Omegal);
    chi0 += dr;
  }

  // Make the table

  dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;

  for(int i=0;i<NZTABLE-1;i++){

    float z   = ZTABLE_INITIAL + (i+0.5)*dztable;
    float chi = Redshift2RadiusTable[i];
    // wkappa = dkappa/dchi = dkappa/dz * H(z) / c
    float wkappa = 3./2.*Omegam*(1+z)*chi*(1-chi/chi0)/r0/r0; 
    table[i] = 0;
    if(chi<chi0) table[i] = wkappa;

  }

} 


//////////////// 21cm kernel table //////////////////

void SetRedshift2WdtbTable(float h, float Omegab, float Omegam, float Omegal, double *table, int evol_flag){
  
  double H0 = h*100.;       // in km/sec/Mpc
  double  c = 2.99792458e5; // in km/sec
  float  r0 = c/H0;         // in Mpc 
  float nu0 = 1.4e3;        // in MHz
  float  T0 = 1.8e2 * h;    // in mK
  float YHe = 0.25;

  // Make the table

  float dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;

  float Tbmean = 1;
  for(int i=0;i<NZTABLE-1;i++){
    

    float z   = ZTABLE_INITIAL + (i+0.5)*dztable;
    if(evol_flag > 0) Tbmean = T0 * OmegaHI(z, Omegab, YHe);
    float wdtb = Tbmean * nu0 / r0 * sqrt(Omegam*pow((1+z),3)+Omegal) ;
    table[i] = wdtb;
  }

}

//////////////// P3D table //////////////////

void SetWavenumber2P3DTable(int nk, double *table, double *wntable, double *ps1, double *ps2, double *ps3, double *ook2){

  // Make the table

  for(int i=0;i<nk;i++){
    
    table[i]=ps1[i]+ps2[i]+ps3[i];
    wntable[i]=sqrt(1/ook2[i]);

  }

}

//////////////// P1D table //////////////////

void SetWavenumber2P1DTable(int nk, double *table, double *wntable, double *ps, double *ook2){

  // Make the table

  for(int i=0;i<nk;i++){
    
    table[i]=ps[i];
    wntable[i]=sqrt(1/ook2[i]);

  }

}

//////////////// Read RedshiftLnM2Float table //////////////////
 
void ReadRedshiftLnM2FloatTable(char *filename, double *table){

  FILE *f = fopen(filename,"rb");

  fread(&RedshiftLnM_nm,      sizeof(int),  1,f);
  fread(&RedshiftLnM_nz,      sizeof(int),  1,f);
  fread(&RedshiftLnM_zinitial,sizeof(float),1,f);
  fread(&RedshiftLnM_zfinal,  sizeof(float),1,f);
  fread(&RedshiftLnM_minitial,sizeof(float),1,f);
  fread(&RedshiftLnM_mfinal,  sizeof(float),1,f);
  
  int ntable = RedshiftLnM_nm * RedshiftLnM_nz;
  
  table = new double[ntable];

  fread(table, sizeof(double), ntable, f);

  fclose(f);

}

//////////////// Read RedshiftLnM2Float table //////////////////
 
void ReadFloat2FloatTable(char *filename, int *nvart, float *varmint, float *varmaxt, double **table){

  FILE *f = fopen(filename,"rb");

  // header
  fread(nvart,   sizeof(int),  1,f);
  fread(varmint, sizeof(float),1,f);
  fread(varmaxt, sizeof(float),1,f);

  // allocate table
  *table = (double *)malloc(*nvart*sizeof(double));

  // read table
  fread(*table, sizeof(double), *nvart, f);

  fclose(f);

}

//////////////// OmegaHI //////////////////
float OmegaHI(float z, float Omegab, float YHe){
  if(z<6){
    return 4e-4 * pow((1+z),0.6);
  } else{
    return (1-YHe)*Omegab;
  }
}
