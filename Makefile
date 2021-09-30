#FC    = f90jx
#FC    = frtpx
FC = ifort
#FC=gfortran

OBJ = mod_variables.o makevort.o
OBJM = main.o
#-O -Utss
OBJS=${OBJ} ${OBJM}
#FFLAGS= -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
#DEBUG
LFLAGS=-shared-intel -convert big_endian -mcmodel=large -mkl -check all -warn all -gen_interfaces -fpe0 -ftrapuv -traceback
#PERFORMANCE
#LFLAGS= -shared-intel -O3 -parallel -qopenmp -convert big_endian -mcmodel=large -lpthread -mkl  -I${MKLROOT}/include 
#FFLAGS =-shared-intel -convert big_endian -mcmodel=large -O3 -llapack -lblas -lgfortran -L/home/yoshimura/lib/lapack-3.8.0/lib -traceback -check all
#FFLAGS=-fconvert=big-endian -mcmodel=large -fbacktrace -fbounds-check -g 
#LFLAGS=-lmkl -llapack -lblas -lgfortran  -L/home/usr/local -L/home/yoshimura/lib/lapack-3.8.0/lib -L/home/morita/MKL
PROGRAM = a.out

all	: ${OBJS}
	 ${FC} -o ${PROGRAM} ${OBJS} ${FFLAGS} ${LFLAGS}

clean :
	rm -f *.o *~ *.mod ${PROGRAM}

%.o: %.f90
	${FC} -c ${FFLAGS} ${LFLAGS} -o $@ $<
