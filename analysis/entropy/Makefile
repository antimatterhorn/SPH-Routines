TARGET=.
SPHDIR=/Users/cody/sph
ARCH=i386-darwin
cc=gcc -g -m32 -I$(SPHDIR)/include \
	-L$(SPHDIR)/Objfiles/$(ARCH) -lsw
fort=gfortran -m32

sdf_entropy : sdf_entropy.c entropy_helm.o
	$(cc) -o sdf_entropy sdf_entropy.c entropy_helm.o -lgfortran

entropy_helm.o : entropy_helm.f90 *.dek
	$(fort) -c entropy_helm.f90

clean :
	rm sdf_entropy *.o
	rm -r sdf_entropy.dSYM
