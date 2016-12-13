
SOURCES := $(wildcard src/*.f90)
SOURCES := $(filter-out src/main.f90, $(SOURCES))
OBJECTS := $(SOURCES:src/%.f90=%.o)

run: machlearn.so
	@gfortran src/main.f90 lib/$< -o $@ -ljsonfortran -I/usr/include -Jinclude

machlearn.so: ${OBJECTS}
	@mkdir -p include
	@mkdir -p lib
	@gfortran -o lib/$@ src/*.o -Wall -llapack -lblas -shared -fPIC

%.o: src/%.f90
	@mkdir -p include
	@gfortran -c -o src/$@ $< -Jinclude -Wall -llapack -lblas -shared -fPIC

clean:
	@rm *.mod
	@rm include/*.mod
	@rmdir include
	@rm src/*.o
	@rm lib/*.so
	@rmdir lib
	@rm run
