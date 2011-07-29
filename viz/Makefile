FILE=
TAR=.
SPHDIR=/Users/cody/sph/tree16
MPI=
ARCH=i386-darwin
cc=gcc -g -m32 -I$(SPHDIR)/include \
	-L$(SPHDIR)/Objfiles/$(ARCH) -lsw $(MPI)

sdf_$(FILE) : sdf_$(FILE).c
	$(cc) -o sdf_$(FILE) sdf_$(FILE).c
	mv sdf_$(FILE) $(TARGET)
