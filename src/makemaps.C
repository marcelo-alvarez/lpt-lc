#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto.h"
#include "allvars.h"

void  ReportMapStatistics(float *, int, char *);
int   SlabShellOverlap(float **, float, float);
float Nu2Redshift(float);
float Redshift2Nu(float);

/* ************* ************* ************* ************* ************* ************* */
/* ************* ************* ************* ************* ************* ************* */

void MakeMaps()
{
  
  int i,j,k,index;
  int i1,j1,k1,i2,j2,k2;
  
  float xm,ym,zm;    // position in map coordinates
  float x,y,z;       // position in cube coordinates
  
  double theta,phi;   // angle in cube coordinates
  float deg2rad=2.*M_PI/360.;

  float Omegam = Parameters.Omegam;
  float Omegab = Parameters.Omegab;
  float Omegal = 1-Omegam;
  float      h = Parameters.h;
  float  zInit = clParameters.zInit;
  float zKappa = clParameters.zKappa;

  float Tcmb = 2.726e6; // Tcmb in muK
  float hoc = h * 1e2 / 3e5; // H0/c in units of 1/Mpc

  float thompson = 6.65e-25;
  float YHe = 0.25;
  float hubble0 = 100./3.086e19*h;
  float ne0=Omegab*h*h*1.88e-29/1.67e-24;
  float c=3e10;

  float dtau0=thompson*ne0*c/hubble0;

  float NHe;

  // Set Radius to Redshift Table
  
  Radius2RedshiftTable = new double[NRTABLE];
  SetRadius2RedshiftTable(h, Omegam, Omegal, Radius2RedshiftTable);

  // Set Redshift to Radius Table
  
  Redshift2RadiusTable = new double[NZTABLE];
  SetRedshift2RadiusTable(h, Omegam, Omegal, Redshift2RadiusTable);

  // Set Redshift to WKappa Table
  
  Redshift2WKappaTable = new double[NZTABLE];
  SetRedshift2WKappaTable(h, Omegam, Omegal, zKappa, Redshift2WKappaTable, 
			  Redshift2RadiusTable);

  // Set Redshift to Wdtb Table
  
  Redshift2WdtbTable = new double[NZTABLE];
  SetRedshift2WdtbTable(h, Omegab, Omegam, Omegal, Redshift2WdtbTable, clParameters.evolve); 

  // Read Redshift to Flux Table
  double *Redshift2FluxPerChiTable;
  int    FluxTable_nz;
  float  FluxTable_zmin, FluxTable_zmax;
  if(Parameters.DoMap[CIBCODE]==1){
  ReadFloat2FloatTable(clParameters.FluxFile,
		       &FluxTable_nz, &FluxTable_zmin, &FluxTable_zmax,
		       &Redshift2FluxPerChiTable);
  }

  // Local copies of parameters

  float BoxSize     = clParameters.BoxSize;
  float Periodicity = clParameters.Periodicity;

  float x0,y0,z0; // position of corner of box so it spans [x0,y0,z0] to [x0+boxsize,y0+boxsize,z0+boxsize], in Mpc
  float xc0,yc0,zc0; // position of first box center
  float xc1,yc1,zc1,xc2,yc2,zc2;

  // set shell and periodicity
  float nu1  = Parameters.nu1;
  float nu2  = Parameters.nu2;
  int   Nnu  = Parameters.Nnu;
  float dnu  = (nu2-nu1)/Nnu;

  float zmin = Parameters.InitialRedshift;
  float zmax = Parameters.FinalRedshift;
  if (clParameters.mapcode == 8){
    zmax = Nu2Redshift(nu1);
    zmin = Nu2Redshift(nu2);
    if(myid==0) printf("zmin=%f zmax=%f nu1=%f nu2=%f\n",zmin,zmax,nu1,nu2);
  }
  float rmin = Redshift2Float(zmin,Redshift2RadiusTable);
  float rmax = Redshift2Float(zmax,Redshift2RadiusTable);
  int nperiodic = (int)(2*rmax / BoxSize + 2);
  int tiledbox = nperiodic * BoxSize ; 
  if(myid==0) printf("rmin = %f rmax = %f nperiodic = %d BoxSize = %f tmapsize=%d \n",rmin,rmax,nperiodic,BoxSize,tmapsize);

  x0  = - tiledbox / 2     ; y0  = - tiledbox / 2     ; z0  = - tiledbox / 2     ;
  xc0 = x0 + BoxSize / 2   ; yc0 = y0 + BoxSize / 2   ; zc0 = z0 + BoxSize / 2 ; 
  
  /*
  xc0 = clParameters.BoxCenter[0] ;
  yc0 = clParameters.BoxCenter[1] ;
  zc0 = clParameters.BoxCenter[2] ;
 
  x0 = clParameters.BoxCenter[0] - BoxSize / 2;
  y0 = clParameters.BoxCenter[1] - BoxSize / 2;
  z0 = clParameters.BoxCenter[2] - BoxSize / 2;
  */
  
  if(myid==0) printf("x0,y0,z0 = %f,%f,%f\n",x0,y0,z0);

  int N=clParameters.N;

  // if (nperiodic <= 2) rmax = fminf(rmax,BoxSize);
  
  int Nxmin = Nlocal * myid;
  float slabsize = BoxSize / nproc;
  float slab_xoffset = slabsize * myid ;

  // local copies of maps

  MPI_Barrier(MPI_COMM_WORLD);

  float *kapmapl, *kszmapl, *taumapl, *dtbmapl, *cibmapl;
  //Local copies of maps
  if(Parameters.DoMap[KAPCODE]==1) kapmapl = (float *)malloc( mapsize*sizeof(float));
  if(Parameters.DoMap[KSZCODE]==1) kszmapl = (float *)malloc( mapsize*sizeof(float));
  if(Parameters.DoMap[TAUCODE]==1) taumapl = (float *)malloc( mapsize*sizeof(float));
  if(Parameters.DoMap[CIBCODE]==1) cibmapl = (float *)malloc( mapsize*sizeof(float));
  if(Parameters.DoMap[DTBCODE]==1) dtbmapl = (float *)malloc(tmapsize*sizeof(float));

  ReportMemory("before map projection",total_local_size,ovrt,oram);

  if(myid==0){
    float dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;    
    FILE* wfile = fopen("weights.dat","w");
    for(int i=0; i<NZTABLE; i++){
      float zcur  = ZTABLE_INITIAL + i*dztable;
      float wkapt = Redshift2WKappaTable[i];
      float wkap = Redshift2Float(zcur,Redshift2WKappaTable);
      float dcur  = growth(zcur,Parameters.Omegam,Parameters.Omegal, Parameters.w)/DInit;
      float wcib = 0;
      if(Parameters.DoMap[CIBCODE]==1)
	wcib = Float2Float(zcur,FluxTable_nz,FluxTable_zmin,FluxTable_zmax,
			   &Redshift2FluxPerChiTable,0.0);
      fprintf(wfile,"%e %e %e %e %e\n",zcur,wcib,wkap,wkapt,dcur);
    }
    fclose(wfile);
  }

  float CellSize = BoxSize / N ;
  float CellVolume = CellSize*CellSize*CellSize;

  // before looping over periodic slabs, find out maximum number of 
  // images in each dimension use 15 Gpc as largest possible radius

  double *vec = new double[3];
  double corner[3];
  float **bb = new float*[2];
  bb[0] = new float[3];
  bb[1] = new float[3];

  // slab corners

  float xs1,ys1,zs1,xs2,ys2,zs2;

  MPI_Barrier(MPI_COMM_WORLD);

  double t1 = MPI_Wtime();

  int iLPT = clParameters.lptcode;

  int Nhalos = halos.N;

  // Each process loops over all periodic slab images overlapping with 
  // integration region. At the end, the maps from all the processes 
  // are added to obtain the final map

  for(int ip=0;ip<nperiodic;ip++){
  for(int jp=0;jp<nperiodic;jp++){
  for(int kp=0;kp<nperiodic;kp++){

    // set corners of octant
    xc1 = xc0 + Periodicity * ( ip - 0.5 );
    yc1 = yc0 + Periodicity * ( jp - 0.5 );
    zc1 = zc0 + Periodicity * ( kp - 0.5 );

    xc2 = xc0 + Periodicity * ( ip + 0.5 );
    yc2 = yc0 + Periodicity * ( jp + 0.5 );
    zc2 = zc0 + Periodicity * ( kp + 0.5 );

    // set corners of slab

    xs1 = x0 + ip*Periodicity + slab_xoffset;
    ys1 = y0 + jp*Periodicity;
    zs1 = z0 + kp*Periodicity;

    xs2 = xs1 + slabsize;
    ys2 = ys1 + BoxSize;
    zs2 = zs1 + BoxSize;

    bb[0][0]=xs1; bb[0][1]=ys1; bb[0][2]=zs1;
    bb[1][0]=xs2; bb[1][1]=ys2; bb[1][2]=zs2;
    
    // if slab does not intersect shell, go to next slab
    if (SlabShellOverlap(bb, rmin, rmax)==0){
//        if (myid==0) printf("Skipping %d %d %d\n",ip,jp,kp);
	continue;
    }
    // Set mask to false
    for(int ic=0;ic<Nlocal;ic++){
      for(int jc=0;     jc<N;jc++){
	for(int kc=0;     kc<N;kc++){
	  int index_dv = ic*N*(N+2) + jc*(N+2) + kc;
	  halomask[index_dv]=0;
	}
      }
    }
    if(myid==0) printf("\n Image %d %d %d\n",ip,jp,kp);
    //MPI_Barrier(MPI_COMM_WORLD);
    if(myid==0) printf("  Masking %d halos\n",Nhalos);
    // Loop over halos and mask
    float fcoll=0;
    float HaloVolumeTot=0;
    for(int ih=0;ih<Nhalos;ih++){

      float xhL = halos.xL[ih]; float yhL = halos.yL[ih]; float zhL = halos.zL[ih]; 
      float xhE = halos.x[ih];  float yhE = halos.y[ih];  float zhE = halos.z[ih]; 
      float RTH = halos.RTH[ih];
      float R200 = RTH / pow(200.,1./3.);

      float x1  = xhL - RTH ;
      float y1  = yhL - RTH ;
      float z1  = zhL - RTH ;

      float x2  = xhL + RTH ;
      float y2  = yhL + RTH ;
      float z2  = zhL + RTH ;

      float HaloVolume = 4./3.*3.14159*RTH*RTH*RTH;

      if( ((x1>xs1 && x1<xs2) || (x2>xs1 && x2<xs2)) &&
	  ((y1>ys1 && y1<ys2) || (y2>ys1 && y2<ys2)) &&
	  ((z1>zs1 && z1<zs2) || (z2>zs1 && z2<zs2))   ){

	HaloVolumeTot += HaloVolume;

	int nshift = 500;
	float big = CellSize*nshift;

	x1+=big;y1+=big;z1+=big;
	x2+=big;y2+=big;z2+=big;

	i1 = (int)((x1-xs1)/CellSize) - nshift;
	j1 = (int)((y1-ys1)/CellSize) - nshift;
	k1 = (int)((z1-zs1)/CellSize) - nshift;

	i2 = (int)((x2-xs1)/CellSize) - nshift;
	j2 = (int)((y2-ys1)/CellSize) - nshift;
	k2 = (int)((z2-zs1)/CellSize) - nshift;

	x1-=big;y1-=big;z1-=big;
	x2-=big;y2-=big;z2-=big;

	// Count overlapping cells
	int OverlapNumber=0;
	for (int i=i1; i<i2; i++){
	  float x = xs1 + (i+0.5)*CellSize;
	  for (int j=j1; j<j2; j++){
	    float y = ys1 + (j+0.5)*CellSize;
	    for (int k=k1; k<k2; k++){
	      float z = zs1 + (k+0.5)*CellSize;
	      
	      float R = pow((x-xhL)*(x-xhL)+(y-yhL)*(y-yhL)+(z-zhL)*(z-zhL),0.5);
	      if(R<RTH+CellSize/1.5) OverlapNumber++;
	    }
	  }
	}

	if(OverlapNumber == 0) printf("ERROR: Halo must overlap at least one cell\n");

	if(i1 <       0) i1 = 0;        if(j1 <  0) j1 = 0  ; if(k1 <  0) k1 = 0  ; 
	if(i2 >= Nlocal) i2 = Nlocal-1; if(j2 >= N) j2 = N-1; if(k2 >= N) k2 = N-1; 

	float CellFill = HaloVolume / CellVolume / OverlapNumber;

	// Loop over grid cells that might overlap 
	for (int i=i1; i<i2; i++){
	  float x = xs1 + (i+0.5)*CellSize;
	  for (int j=j1; j<j2; j++){
	    float y = ys1 + (j+0.5)*CellSize;
	    for (int k=k1; k<k2; k++){
	      float z = zs1 + (k+0.5)*CellSize;
	      
	      float R = pow((x-xhL)*(x-xhL)+(y-yhL)*(y-yhL)+(z-zhL)*(z-zhL),0.5);

	      if(R<RTH+CellSize/1.5){
		float frad = pow(R/RTH,3);
		int index_mask = i*N*(N+2) + j*(N+2) + k;
		sx1[index_mask] = (xhE-x) * frad;
		sy1[index_mask] = (yhE-y) * frad;
		sz1[index_mask] = (zhE-z) * frad;
		halomask[index_mask] += CellFill;
		fcoll += CellFill;
	      }
	    }
	  }
	}

      }      
    }

    if(myid==0) printf("  Projecting LPT field\n");
    // Note that x and z are switched
    for(int ic=0;ic<Nlocal;ic++){
      float xL = xs1 + (ic+0.5)*CellSize;
      if(xL>xc1 && xL<xc2){
    for(int jc=0;jc<N;jc++){
      float yL = ys1 + (jc+0.5)*CellSize;
      if(yL>yc1 && yL<yc2){
    for(int kc=0;kc<N;kc++){
      float zL = zs1 + (kc+0.5)*CellSize;
      if(zL>zc1 && zL<zc2){
      
      float r = sqrt(xL*xL + yL*yL + zL*zL);

      if(r    < rmin || r > rmax) continue;

      float zcur   = Radius2Float(r,Radius2RedshiftTable);

      int index_dv = ic*N*(N+2) + jc*(N+2) + kc;

      // Use zreion field to determine if point in Lagrangian space is ionized
      float xHI = 1;
      if(clParameters.zmask == 1 && zmask[index_dv]>zcur) xHI=0;
      
      // CMB Lensing redshift factor
      float kapfac;
      float Wkap;
      if(Parameters.DoMap[KAPCODE]==1){
	Wkap = Redshift2Float(zcur,Redshift2WKappaTable);
	kapfac = Wkap * 
	  pow(CellSize,3) / pow(r,2) * mapsize / 4. / 3.14159;
      }

      // CIB redshift factor
      float cibfac;
      float Wcib;
      if(Parameters.DoMap[CIBCODE]==1){
	Wcib = Float2Float(zcur,FluxTable_nz,FluxTable_zmin,FluxTable_zmax,
			   &Redshift2FluxPerChiTable,0.0);
	cibfac = Wcib * 
	  pow(CellSize,3) / pow(r,2) * mapsize / 4. / 3.14159;
      }

      // 21cm redshift factors
      float dtbfac;
      float   Wdtb;
      float     nu;
      int      inu;
      if(Parameters.DoMap[DTBCODE]==1){
	Wdtb   = Redshift2Float(zcur,Redshift2WdtbTable);
	nu     = Redshift2Nu(zcur);
	dtbfac = Wdtb * xHI * 
	  pow(CellSize,3) / pow(r,2) * mapsize / 4. / 3.14159 / dnu;
	inu    = (int)((nu-nu1)/dnu);
	if(inu<0 || inu>= Nnu){ inu=0; dtbfac = 0; }
      }
      
      float D      = growth(zcur,Parameters.Omegam,Parameters.Omegal, Parameters.w)/DInit;
      if(clParameters.evolve == 0) D = 1;
      float D2     = 3. / 7. * D * D;

      if(iLPT < 1) D  = 0;
      if(iLPT < 2) D2 = 0;

      long pixel, tpixel;
      
      float xE,yE,zE;
      // Displacements

      float sx,sy,sz;

      sx = D * sx1[index_dv] + D2 * sx2[index_dv] ;
      sy = D * sy1[index_dv] + D2 * sy2[index_dv] ;
      sz = D * sz1[index_dv] + D2 * sz2[index_dv] ;

      if(halomask[index_dv]>0){
	sx = sx1[index_dv] ;
	sy = sy1[index_dv] ;
	sz = sz1[index_dv] ;
      }

      xE = xL + sx ;
      yE = yL + sy ; 
      zE = zL + sz ;

      // Subtract contribution from Lagrangian point for lensing convergence
      if(Parameters.DoMap[KAPCODE]==1){
	vec[2] = xL; vec[1] = yL; vec[0] = zL;
	vec2pix_nest(Parameters.NSide, vec, &pixel);
	kapmapl[pixel] -= kapfac ;
      }
      
      // Add contribution to Eulerian point      
      vec[2] = xE; vec[1] = yE; vec[0] = zE; // x-z flipped
      vec2pix_nest(Parameters.NSide, vec, &pixel); tpixel = pixel + inu * mapsize;      
      if(Parameters.DoMap[KAPCODE]==1) kapmapl[ pixel] += kapfac ; // * (1-halomask[index_dv]);
      if(Parameters.DoMap[CIBCODE]==1) cibmapl[ pixel] += cibfac ; // * (1-halomask[index_dv]);
      if(Parameters.DoMap[DTBCODE]==1) dtbmapl[tpixel] += dtbfac ; // * (1-halomask[index_dv]);
      
    }
    }
    }    

    }
    }
    }    

  }
  }
  }    
  MPI_Barrier(MPI_COMM_WORLD);
  double dt = MPI_Wtime() - t1;
  if(myid==0) printf("\n Projection took %le seconds\n",dt);

  // sum process contributions

  if(Parameters.DoMap[KAPCODE]==1)
    MPI_Reduce(kapmapl, kapmap, mapsize, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(Parameters.DoMap[KSZCODE]==1)
    MPI_Reduce(kszmapl, kszmap, mapsize, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(Parameters.DoMap[TAUCODE]==1)
    MPI_Reduce(taumapl, taumap, mapsize, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(Parameters.DoMap[CIBCODE]==1)
    MPI_Reduce(cibmapl, cibmap, mapsize, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(Parameters.DoMap[DTBCODE]==1)
    MPI_Reduce(dtbmapl, dtbmap,tmapsize, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(myid==0) printf("\n Sum process contributions complete\n");

  // report statistics

  if(Parameters.DoMap[KAPCODE]==1 && myid==0) ReportMapStatistics(kapmap, mapsize," kappa    ");
  if(Parameters.DoMap[KSZCODE]==1 && myid==0) ReportMapStatistics(kszmap, mapsize," kSZ      ");
  if(Parameters.DoMap[TAUCODE]==1 && myid==0) ReportMapStatistics(taumap, mapsize," tau      ");
  if(Parameters.DoMap[CIBCODE]==1 && myid==0) ReportMapStatistics(cibmap, mapsize," CIB      ");
  if(Parameters.DoMap[DTBCODE]==1 && myid==0) ReportMapStatistics(dtbmap,tmapsize," 21-cm dTb");

  if (myid==0) printf("\n Report map stats complete \n");

}

int SlabShellOverlap(float **bb, float rmin, float rmax){

  for(int i1=0;i1<2;i1++){
  for(int j1=0;j1<2;j1++){
  for(int k1=0;k1<2;k1++){

    float rc1 = pow((pow(bb[i1][0],2)+pow(bb[j1][1],2)+pow(bb[k1][2],2)),0.5);

    if(rc1>rmin && rc1<rmax) return 1;

    float dr1ci = rc1 - rmin;
    float dr1co = rc1 - rmax;

    for(int i2=0;i2<2;i2++){
    for(int j2=0;j2<2;j2++){
    for(int k2=0;k2<2;k2++){

      float rc2 = pow((pow(bb[i2][0],2)+pow(bb[j2][1],2)+pow(bb[k2][2],2)),0.5);

      // inner boundary
      float dr2ci = rc2 - rmin;
      if(dr1ci*dr2ci < 0) return 1;

      // outer boundary
      float dr2co = rc2 - rmax;
      if(dr1co*dr2co < 0) return 1;

    }    
    }    
    }    
  }
  }
  }

  return 0;
  
}

float Nu2Redshift(float nu){
  float nu0 = 1.4e3;
  return nu0/nu - 1;
}

float Redshift2Nu(float z){
  float nu0 = 1.4e3;
  return nu0/(1+z);
}


void ReportMapStatistics(float *map, int mapsize, char *variable){
  
  int i;
  double  mean, var, rms, sum;
  
  // mean
  sum=0; for(i=0;i<mapsize;i++) sum += map[i];
  mean = sum / mapsize;

  // rms
  sum=0; for(i=0;i<mapsize;i++) sum += pow((map[i]-mean),2);
  var  = sum / mapsize;
  rms  = pow(var,0.5);

  // report
  printf("\n %s mean, rms = %e, %e\n",variable,mean,rms);

}



