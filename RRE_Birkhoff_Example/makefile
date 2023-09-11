FC=gfortran
FFLAGS = 
LIBS = -L$HOME/libf77/$ARCH -llapack -lblas

main.exe: birkhoff.o main.o
	$(FC) $(FFLAGS) birkhoff.o main.o $(LIBS) -o main.exe

birkhoff.o: birkhoff.f90
	$(FC) $(FFLAGS) -c birkhoff.f90

main.o: main.f90
	$(FC) $(FFLAGS) -c main.f90

run: main.exe
	./main.exe

clean:
	rm *.o *.exe *.mod