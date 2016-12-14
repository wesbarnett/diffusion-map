.PHONY: clean

SOURCES := $(wildcard src/*.f90)
SOURCES := $(filter-out src/main.f90, $(SOURCES))
OBJECTS := $(SOURCES:src/%.f90=%.o)
LDFLAGS += `pkg-config --libs lapack` 
CFLAGS += -shared -fPIC -Wall `pkg-config --cflags lapack`

run: machlearn.so
	@gfortran -o $@ lib/$< src/main.f90  -ljsonfortran -I/usr/include -Jinclude

machlearn.so: ${OBJECTS}
	@mkdir -p include
	@mkdir -p lib
	@gfortran -o lib/$@ src/*.o  ${LDFLAGS} ${CFLAGS}

%.o: src/%.f90
	@mkdir -p include
	@gfortran -c -o src/$@ $< -Jinclude ${LDFLAGS} ${CFLAGS}

clean:
	@rm -f *.mod
	@rm -f src/*.o
	@rm -rf include
	@rm -rf lib
	@rm -f run
