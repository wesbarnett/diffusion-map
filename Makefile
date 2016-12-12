
CFLAGS := -Wall
LDFLAGS := -llapack -lblas -ljsonfortran

run: 
	@mkdir -p include
	@gfortran src/main.f90 -o run ${CFLAGS} ${LDFLAGS} -I/usr/include -Jinclude

clean:
	@rm include/*.mod
	@rmdir include
	@rm run
