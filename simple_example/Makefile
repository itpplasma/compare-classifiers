all: example.x wout.nc

example.x: example.f90
	gfortran -o example.x -J ../SIMPLE/build/SRC example.f90 \
		-L ../SIMPLE/build -lsimple

wout.nc:
	wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc

clean:
	rm -f example.x *.mod *.nc
