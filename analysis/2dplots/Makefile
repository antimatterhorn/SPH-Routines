SPHDIR=/Users/cody/sph/tree16
ARCH=i386-darwin
cc=gcc -m32 -g -I$(SPHDIR)/include \
	-L$(SPHDIR)/Objfiles/$(ARCH) -lsw
pp=g++ -m32
fort=/usr/local/bin/gfortran -m32

2dconduct : public_helm.o 2dconduct.c
	$(cc) -o 2dconduct 2dconduct.c public_helm.o -L/usr/local/lib/i386 -lgfortran


public_helm.o : public_helm.f90 *.dek
	$(fort) -c public_helm.f90

clean :
	rm 2dconduct *.o