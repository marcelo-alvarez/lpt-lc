#include <math.h>
#include "proto.h"
#include "allvars.h"
#include <unistd.h>

void usage(){

  if(myid==0){
    printf("\n usage: asm <parameterfile> [options]\n");
    
    printf("\n OPTIONS:\n");
    printf("\n   -h show this message");
    printf("\n   -v verbose            [default = OFF]"); 
    printf("\n   -D density file       [default = './delta'");
    printf("\n   -R redshift mask file [default = 'NULL']");
    printf("\n   -H halo file          [default = 'NULL']");
    printf("\n   -F flux(z,M) file     [default = 'NULL']");
    printf("\n   -o Output base name   [default = 'map']");
    printf("\n   -N grid dimension     [default = 512]");
    printf("\n   -C Number of chunks   [default = 1]");
    printf("\n   -B box size in Mpc    [default = 100 Mpc]");
    printf("\n   -p periodicity in Mpc [default = 100 Mpc]");
    printf("\n   -x center x in Mpc    [default = 100 Mpc]");
    printf("\n   -y center y in Mpc    [default = 100 Mpc]");
    printf("\n   -z center z in Mpc    [default = 100 Mpc]");
    printf("\n   -r initial redshift   [default = 0]");
    printf("\n   -l LPT code, [l]LPT   [default = 1 --> 1LPT]");
    printf("\n   -k lensing source z   [default = 1100]");
    printf("\n   -e evolution code     [default = 1 --> evolution");
    printf("\n   -b binary map only    [default = 0 --> fits and binary]");
    printf("\n   -m binary map code (e.g. kap=0 --> no ; kap=1 --> yes)");
    printf("\n       = kap * 1  (CMB lensing)");
    printf("\n       + ksz * 2  (kinetic SZ) ");
    printf("\n       + tau * 4  (tau_es)     ");
    printf("\n       + dtb * 8  (21cm)       ");
    printf("\n       + cib * 16 (CIB)        ");
    printf("\n\n");
  }

  MPI_Finalize();
  exit(0);

}

void CommandLine(int argc, char *argv[])
{

  int c;

  opterr=0;
  
  sprintf(clParameters.DeltaFile,"./input");
  sprintf(clParameters.RedshiftFile,"NULL");
  sprintf(clParameters.HaloFile,"NULL");
  sprintf(clParameters.FluxFile,"NULL");
  sprintf(clParameters.BaseOut,"./output");

  if(myid==0) system("[ -d maps ] || mkdir maps");
  
  clParameters.verbose      = 0;
  clParameters.deltain      = 1;
  clParameters.zmask        = 0;
  clParameters.halomask     = 0;  
  clParameters.BoxSize      = 100;
  clParameters.Nchunk       = 1;
  clParameters.Periodicity  = 100;
  clParameters.BoxCenter[0] = 0;
  clParameters.BoxCenter[1] = 0;
  clParameters.BoxCenter[2] = 0;
  clParameters.N            = 512;
  clParameters.zInit        = 0.;
  clParameters.lptcode      = 1; 
  clParameters.mapcode      = 0; 
  clParameters.zKappa       = 1100;
  clParameters.evolve       = 1;
  clParameters.binary_only  = 0;
  
  while ((c = getopt (argc, argv, "hvP:D:R:H:F:o:N:B:p:C:x:y:z:r:l:m:k:e:b:")) != -1)
    switch (c)
      {
      case 'h':
	usage();
	break;
      case 'v':
	clParameters.verbose = 1;
	break;
      case 'P':
	sprintf(clParameters.ParamFile,"%s",optarg);
	break;
      case 'D':
	sprintf(clParameters.DeltaFile,"./%s",optarg);
	break;
      case 'R':
	sprintf(clParameters.RedshiftFile,"./%s",optarg);
	clParameters.zmask = 1;
	break;
      case 'H':
	sprintf(clParameters.HaloFile,"./%s",optarg);
	clParameters.halomask = 1;
	break;
      case 'F':
	sprintf(clParameters.FluxFile,"./%s",optarg);
	break;
      case 'o':
	sprintf(clParameters.BaseOut,"./%s",optarg);
	break;
      case 'N':
	clParameters.N = atoi(optarg);
	break;
      case 'B':
	clParameters.BoxSize = atof(optarg);
	break;
      case 'C':
	clParameters.Nchunk = atoi(optarg);
	break;
      case 'p':
	clParameters.Periodicity = atof(optarg);
	break;
      case 'x':
	clParameters.BoxCenter[0] = atof(optarg);
	break;
      case 'y':
	clParameters.BoxCenter[1] = atof(optarg);
	break;
      case 'z':
	clParameters.BoxCenter[2] = atof(optarg);
	break;
      case 'r':
	clParameters.zInit = atof(optarg);
	break;
      case 'l':
	clParameters.lptcode = atoi(optarg);
	break;
      case 'm':
	clParameters.mapcode = atoi(optarg);
	break;
      case 'k':
	clParameters.zKappa = atof(optarg);
	break;
      case 'e':
	clParameters.evolve = atoi(optarg);
	break;
      case 'b':
	clParameters.binary_only = atoi(optarg);
	break;
      case '?':
	if (optopt == 'i'){
	  if(myid==0) fprintf (stderr, "\n Option -%c requires an argument.\n", optopt);
	  usage();
	}
	if (optopt == 'o'){
	  if(myid==0) fprintf (stderr, "\n Option -%c requires an argument.\n", optopt);
	  usage();
	}
	else if (isprint (optopt)){
	  if(myid==0) fprintf (stderr, "\n Unknown option `-%c'.\n", optopt);
	  usage();
	}
	else{
	  if(myid==0) fprintf (stderr,
		   "Unknown option character `\\x%x'.\n",
		   optopt);
	  usage();
	}
	return;
      default:
	usage();
      }

  if(clParameters.verbose == 0){
    // Redirect stdout to output file
    char fname[256];
    sprintf(fname,"%s.stdout",clParameters.BaseOut);
    freopen(fname,"w",stdout);
  }

  Parameters.CurrentCode = NUMMAPCODES;
  FillMapCodeArray(2*clParameters.mapcode,NUMMAPCODES);
  int *buff = new int[NUMMAPCODES];
  for (int j=0;j<NUMMAPCODES;j++) buff[j]=Parameters.DoMap[j];
  for (int j=0;j<NUMMAPCODES;j++) Parameters.DoMap[j]=buff[NUMMAPCODES-j-1];

  // Don't buffer stdout and stderr
  setvbuf(stdout, NULL, _IONBF, 0);
  setvbuf(stderr, NULL, _IONBF, 0);

  if(myid==0){
    printf("\n Command line:");
    for(int i=0;i<argc;i++) printf(" %s",argv[i]); printf("\n");  
  }

  Parameters.N = clParameters.N;
  Parameters.BoxSize = clParameters.BoxSize;
  return;

}

void FillMapCodeArray(int n, int m)
{
  if (n < 0 || n > pow(2,NUMMAPCODES+1)){
    if(myid==0) fprintf (stderr, "\n map code %d not an allowed value\n",n);
    MPI_Barrier(MPI_COMM_WORLD);
    usage();
  }
  if (n / 2 != 0) {
    FillMapCodeArray(n / 2,m-1);
  }
  Parameters.DoMap[m] = n % 2;
}
