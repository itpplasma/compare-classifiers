all: example.x

example.x: example.f90
	gfortran -o example.x -J ../SIMPLE/build/SRC example.f90 \
		-L ../SIMPLE/build -lsimple

clean:
	rm -f example.x *.mod
