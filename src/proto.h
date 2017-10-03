#ifndef PROTO_H
#define PROTO_H

#include <stdio.h>
#include "cosmology.h"
#include "geometry.h"
#include "arrayoperations.h"
#include "memorytracking.h"
#include "parallel_io.h"

#ifndef ALLVARS_H
#include "allvars.h"
#endif

// arrayoperations.C
void AllocateArrays(); 

// commandline.C
void usage();
void CommandLine(int, char **);
void FillMapCodeArray(int, int);

// lpt.C
void Displace_1LPT(float *, float *, float *, float *);
void Displace_2LPT(float *, float *, float *, float *, float *, float *, float *, float *);
void GenerateMomentum();
void GenerateCurlMomentum();
void GenerateCurlMomentumRealSpace();
int  ijk2index(int, int, int);

// io.C
void ReadGridFromFile(float *, char *);
void ReadHaloFile(char *);
void ReadParameterFile();
void WriteMaps();
void WriteLPT();

// makemaps.C
void MakeMaps();

#endif

