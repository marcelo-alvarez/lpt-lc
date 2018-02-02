#----------------------------------------------------------------------
# OPTIONS FOR RUNNING ON OSCAR MACHINES
#  Be sure to load module mpich/3.1.1 and gsl/2.3 before running Make!
#----------------------------------------------------------------------

OPTIMIZE =  -O4 -w 

MODFLAG = -module 
OMPLIB = -openmp

FFTW2_PATH = $(HOME)/21cmSimulation/C/fftw-2.1.5
CFITS_PATH = $(HOME)/21cmSimulation/C/cfitsio
GSL_PATH = $(HOME)/21cmSimulation/C/gsl-2.4
CURL_PATH = /gpfs/runtime/opt/curl/7.50.3/include

CC = mpicc
C++ = mpic++

OPTIONS = -w

#----------------------------------------------------------------------
# DO NOT MODIFY THIS FILE
#----------------------------------------------------------------------

FFTLIB = -L$(FFTW2_PATH)/lib -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
FFTINC = -I$(FFTW2_PATH)/include 

#FFTLIB = -L$(FFTW2_PATH)/lib -L$(FFTW3_PATH)/lib -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
#FFTINC = -I$(FFTW2_PATH)/include -I$(FFTW3_PATH)/include 

CFTLIB = -L$(CFITS_PATH)/lib -lcfitsio
CFTINC = -I$(CFITS_PATH)/include

GSLLIB = -L$(GSL_PATH)/lib -lgsl -lgslcblas
GSLINC = -I$(GSL_PATH)/include

CURLINC = -I$(CURL_PATH)/include
CURLLIB = -L$(CURL_PATH)/lib -lcurl

LIB = $FFTLIB $CFTLIB $GSLLIB $CURLLIB
INC = $FFTINC $CFTINC $GSLINC $CURLINC
 
# OBJECT FILES
srcdir  = ./src
objs = \
     $(srcdir)/allocate.o \
     $(srcdir)/allvars.o \
     $(srcdir)/commandline.o \
     $(srcdir)/lpt.o \
     $(srcdir)/geometry.o \
     $(srcdir)/main.o \
     $(srcdir)/makemaps.o \
     $(srcdir)/parallel_io.o \
     $(srcdir)/io.o \
     $(srcdir)/chealpix.o \
     $(srcdir)/arrayoperations.o \
     $(srcdir)/cosmology.o \
     $(srcdir)/math.o \
     $(srcdir)/memorytracking.o \
     $(srcdir)/tables.o 

bindir = ./bin

COMPILE_FLAGS = $(OPTIMIZE) $(FFTINC) $(GSLINC) $(CFTINC) $(OPTIONS) $(CURLINC) -DENABLE_FITSIO 
LINK_FLAGS    = $(OPTIMIZE) $(FFTLIB) $(GSLLIB) $(CFTLIB) $(CURLLIB)

EXEC = lin2map

OBJS     = $(objs)

.SUFFIXES: .o .f .f90 .F90 .c .C

$(srcdir)/%.o: $(srcdir)/%.c
	$(CC)  $(COMPILE_FLAGS) -c $< -o $@
$(srcdir)/%.o: $(srcdir)/%.C 
	$(C++) $(COMPILE_FLAGS) -c $< -o $@

$(EXEC): $(OBJS) 
	$(C++) $(OBJS) $(LINK_FLAGS) -o  $(bindir)/$(EXEC)  

clean:
	rm -f $(EXEC) $(OBJS) 


