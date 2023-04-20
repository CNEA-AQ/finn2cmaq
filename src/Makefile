.SUFFIXES: .o .f90

FC   = gfortran
LIBS = -L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm 
INC  = -I/usr/include
FFLAGS = -O2 -ffree-line-length-none #-Wunused 
OBJS = finn2cmaq.o
EXE  = ../exe/finn2cmaq

.f90.o:
		${FC} ${FFLAGS} -c ${INC} $<

finn2cmaq: ${OBJS}
		${FC} -o ${EXE} ${OBJS} ${LIBS} 

clean:
		rm -f *.o *.mod
