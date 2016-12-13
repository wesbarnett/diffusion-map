
run: 
	@mkdir -p include
	@gfortran -c -o src/diffusion_map.o src/diffusion_map.f90 -Jinclude -Wall -llapack -lblas
	@gfortran -c -o src/princ_comp.o src/princ_comp.f90 -Jinclude -Wall -llapack -lblas
	@gfortran src/main.f90 src/*.o -o run -ljsonfortran -I/usr/include -Jinclude -llapack -lblas

clean:
	@rm include/*.mod
	@rm src/*.o
	@rmdir include
	@rm run
