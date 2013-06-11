FILE=
TAR=.
SPHDIR=/Users/cody/sph/tree16
ARCH=i686-darwin
cc=gcc -g -m64 -I$(SPHDIR)/include \
	-L$(SPHDIR)/Objfiles/$(ARCH) -lsw


single : sdf_$(FILE).c
		$(cc) -o sdf_$(FILE) sdf_$(FILE).c
		mv sdf_$(FILE) $(TARGET)


parallel : sdf_$(FILE).c
		$(cc) -o sdf_$(FILE) sdf_$(FILE).c -fopenmp
		mv sdf_$(FILE) $(TARGET)
	