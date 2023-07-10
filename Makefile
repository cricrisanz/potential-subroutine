# Linux compiler options
FC= gfortran -std=legacy
#FC= gfortran 
#FLAGS = -g

#Lapack and Blas in Bohr
#LDLIBS	= -L/home/cristina/libs/lapack-3.1.1 -llapack -lblas
#LDLIBS	= /home/cris/libs/lapack-3.1.1/liblapack.a
#LDLIBS2	= /home/cris/libs/lapack-3.1.1/libblas.a
#LDLIBS2	= /usr/lib/libblas.so.3gf
#MODS	= commons.f90 summations.f90 linear.f90 optimisation.f90

SRC	= main.f90 pot.f90

fit.exe : $(SRC) 
	$(FC) $(SRC) -o pot.exe 

cleanall:
	rm -f *.o *.mod *~ *.exe 

cleandat:
	rm s??-??.dat fort.*

clean:
	rm -f *.o *.mod *~ *.exe fort.*
