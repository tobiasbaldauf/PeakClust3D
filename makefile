SOURCE = pksampler_cross
#SOURCE = pksampler

F90 = gfortran
F90FLAGS = -O3 -mcmodel=medium -fopenmp -ffixed-line-length-none -ffree-line-length-none

all:	${SOURCE}.f90
	${F90} ${F90FLAGS} ${SOURCE}.f90 -o ${SOURCE}.exe
