all: aout

aout: mod_afssh.o AFSSH.o
	ifort -o aout mod_afssh.o AFSSH.o -O2 -mkl

%.o: %.f90
	ifort -c $<

clean:
	rm *.o aout

