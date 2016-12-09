all:
	gfortran main.f90 -llapack -lblas -Wall -ljsonfortran -I/usr/include
