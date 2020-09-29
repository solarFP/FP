.SUFFIXES:
.SUFFIXES: .o .F .F90 .f

#SET PATH TO HDF5 directory: EDIT THESE TO POINT TO THE CORRECT PATH IN YOUR INSTALLATION OF HDF5
HDF5_LIB_PATH = ~/hdf5/lib
HDF5_INCLUDE_PATH = ~/hdf5/include

#OPTIONS FOR INTEL COMPILER
#F = mpiifort
#Use the following FOPT for maximum optimizations
#FOPTS = -O3 -xHost -fpe0 -traceback -I ${HDF5_INCLUDE_PATH} -heap-arrays 1
#Use the following FOPT for debug mode
#FOPTS = -g -fpe0 -traceback -check bounds,uninit,pointers, -heap-arrays 1 -I ${HDF5_INCLUDE_PATH}
#LOPTS = -L ${HDF5_LIB_PATH} -lhdf5_fortran

#GFORTRAN COMPILER
F = mpif90
#Use the following FOPT for maximum optimizations
FOPTS = -O3 -march=native -ffree-line-length-none -mcmodel=medium -I ${HDF5_INCLUDE_PATH}
#Use the following FOPT for debug mode
#FOPTS = -Wall -g -fcheck=all -ffree-line-length-none -ffpe-trap=invalid,zero,overflow  -fbacktrace -I ${HDF5_INCLUDE_PATH}
LOPTS = -L ${HDF5_LIB_PATH} -lhdf5_fortran

objects = const.o bessk.o helpers.o options.o derivs.o grid.o parmpi.o loop.o beam.o matrix.o pentdag.o fpinopts.o fp.o 

standalone = hdf5readwrite.o writer.o fpsa.o

fp : libfp.a $(standalone)
	${F} ${LOPTS} -o fp $(standalone) libfp.a

libfp.a : $(objects)
	ar cr libfp.a $(objects)

.F.o .f.o .F90.o :
	${F} ${FOPTS} -c -o $@ $<

.phony : clean
clean :
	rm fp $(objects) $(standalone) *.mod
