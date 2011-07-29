TARGET=.
SPHDIR=/Users/cody/sph/tree16
ARCH=i386-darwin
cc=gcc -g -m32 -I$(SPHDIR)/include \
	-L$(SPHDIR)/Objfiles/$(ARCH) -lsw
fort=gfortran -m32

sdf_match : sdf_match.c chop_aprox13.o public_helm.o
	$(cc) -o sdf_match sdf_match.c chop_aprox13.o public_helm.o -lgfortran
	mv sdf_match $(TARGET)

chop_aprox13.o : chop_aprox13.f90 *.dek
	$(fort) -c chop_aprox13.f90

public_helm.o : public_helm.f90 *.dek
	$(fort) -c public_helm.f90

clean :
	rm sdf_match *.o
