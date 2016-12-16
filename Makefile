.PHONY: clean

SOURCES := $(wildcard src/*.f90)
SOURCES := $(filter-out src/main.f90, $(SOURCES))
OBJECTS := $(SOURCES:src/%.f90=%.o)
CFLAGS  = -shared -fPIC -Wall 

run: machlearn.so
	@gfortran -o $@ src/main.f90  lib/$< -I/usr/include -Iinclude -Jinclude -lgmxfort -ljsonfortran -fopenmp

machlearn.so: ${OBJECTS}
	@mkdir -p include
	@mkdir -p lib
	@gfortran -o lib/$@ src/*.o  ${CFLAGS} -llapack -lblas

%.o: src/%.f90
	@mkdir -p include
	@gfortran -c -o src/$@ $< -Jinclude ${CFLAGS} -llapack -lblas

clean:
	@rm -f *.mod
	@rm -f src/*.o
	@rm -rf include
	@rm -rf lib
	@rm -f run
