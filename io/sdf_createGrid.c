#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
//#include "singlio.h"
#include "Msgs.h"
#include <stdlib.h>

#define SPHOUTBODYDESC \
"struct {\n\
	double x, y, z;		/* position of body */\n\
	float mass;			/* mass of body */\n\
	float h;			/* smoothing length of body */\n\
}"

typedef struct {
	double x, y, z;             /* position of body */
	float mass;           /* mass of body */
	float h;              /* smoothing length of body */
} SDFgrid;

int main()
{
	int iter = 0;
	float dt = 0;
	float eps = 0;
	float Gnewt = 0;
	float tolerance = 0;
	float frac_tolerance = 0;
	int ndim = 3;
	float tpos = 0;
	float tvel = 0;
	float gamma = 0;
	double ke = 0;
	double pe = 0;
	double te = 0;

	char name [80];
	int nobj,gnobj,npts;
	int i,j,k,p;
	double gridsize,dgrid;
	
	
	printf("filename: ");
	gets (name);
	printf("gridpoints on a side: ");
	scanf("%d",&npts);
	printf("gridsize: ");
	scanf("%lf",&gridsize);	

	dgrid = gridsize*2 / npts;	
	gnobj = nobj = npts*npts*npts;
	

	SDFgrid *outArray = malloc(sizeof(SDFgrid) * nobj);
	
	for(i=0;i<npts;i++){
		for(j=0;j<npts;j++){
			for(k=0;k<npts;k++){
				p = i*npts*npts + j*npts + k;
				outArray[p].x		= -gridsize + dgrid*i + dgrid/2;
				outArray[p].y		= -gridsize + dgrid*j + dgrid/2;
				outArray[p].z		= -gridsize + dgrid*k + dgrid/2;
				outArray[p].h		= 3.0*dgrid;
				outArray[p].mass	= 1.0;
			}
		}
	}

	SDFwrite(name, gnobj, 
			 nobj, outArray, sizeof(SDFgrid),
			 SPHOUTBODYDESC,
			 "npart", SDF_INT, gnobj,
			 "iter", SDF_INT, iter,
			 "dt", SDF_FLOAT, dt,
			 "eps", SDF_FLOAT, eps,
			 "Gnewt", SDF_FLOAT, Gnewt,
			 "tolerance", SDF_FLOAT, tolerance,
			 "frac_tolerance", SDF_FLOAT, frac_tolerance,
			 "ndim", SDF_INT, ndim,
			 "tpos", SDF_FLOAT, tpos,
			 "tvel", SDF_FLOAT, tvel,
			 //"t_wind", SDF_FLOAT, twind_out,
			 //"R0", SDF_FLOAT, output_R0,
			 //"Omega0", SDF_FLOAT, cosmo.Omega0,
			 //"H0", SDF_FLOAT, cosmo.H0,
			 //"Lambda_prime", SDF_FLOAT, cosmo.Lambda,
			 //"hubble", SDF_FLOAT, output_h,
			 //"redshift", SDF_FLOAT, output_z,
			 "gamma", SDF_FLOAT, gamma,
			 //"centmass", SDF_FLOAT, centmass, 
			 "ke", SDF_DOUBLE, ke,
			 "pe", SDF_DOUBLE, pe,
			 "te", SDF_DOUBLE, te,
			 NULL);
	free(outArray);
	return 0;
}
