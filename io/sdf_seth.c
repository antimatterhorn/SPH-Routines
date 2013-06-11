#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
//#include "singlio.h"
#include "Msgs.h"
#include <stdlib.h>
#include <math.h>
#include "../SPHbody.h"

#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)

int main(int argc, char **argv[])
{
	int i,gnobj,nobj,conf,iter,ndim=3;
	float dt = 0;
	float Gnewt = 0;
	float tpos = 0;
	char outname [80];
	double xmax=0,vol;
	float h;

	
	if (argc < 2){
		printf("SDF file: ");
		gets (argv[1]);
	}
	singlPrintf("Reading \"%s\"\n", argv[1]);

	SDF *sdfp;
	SPHbody *body;

	sdfp = SDFreadf(argv[1], (void **)&body, &gnobj, &nobj, sizeof(SPHbody),
					"x", offsetof(SPHbody, x), &conf,
					"y", offsetof(SPHbody, y), &conf,
					"z", offsetof(SPHbody, z), &conf,
					"mass", offsetof(SPHbody, mass), &conf,
					"vx", offsetof(SPHbody, vx), &conf,
					"vy", offsetof(SPHbody, vy), &conf,
					"vz", offsetof(SPHbody, vz), &conf,
					"u", offsetof(SPHbody, u), &conf,
					"h", offsetof(SPHbody, h), &conf,
					"rho", offsetof(SPHbody, rho), &conf,
					"drho_dt", offsetof(SPHbody, drho_dt), &conf,
					"udot", offsetof(SPHbody, udot), &conf,
					"temp", offsetof(SPHbody, temp), &conf,
					"abar", offsetof(SPHbody, abar), &conf,
					"zbar", offsetof(SPHbody, zbar), &conf,
					"phi", offsetof(SPHbody, phi), &conf,
					"idt", offsetof(SPHbody, idt), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					NULL);
	//now to fill the misc. floats and ints//
	SDFgetfloatOrDefault(sdfp, "dt",  &dt, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "Gnewt",  &Gnewt, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);

	for(i=0;i<nobj;i++) xmax = MAX(xmax,body[i].x);
	
	printf("found xmax = %3.2e\n",xmax);
	
	h = xmax/pow((double)nobj,1.0/3.0);
	vol = (4.0/3.0*3.14159*pow(h,3.0));
	
	for(i=0;i<nobj;i++)
	{
		body[i].h = h;
		body[i].rho = body[i].mass/vol;
	}
	
	snprintf(outname, sizeof(outname), "h_%s", argv[1]);
	SDFwrite(outname, gnobj, 
			 nobj, body, sizeof(SPHbody),
			 SPHOUTBODYDESC,
			 "npart", SDF_INT, gnobj,
			 "iter", SDF_INT, iter,
			 "dt", SDF_FLOAT, dt,
			 "Gnewt", SDF_FLOAT, Gnewt,
			 "ndim", SDF_INT, 3,
			 "tpos", SDF_FLOAT, tpos,
			 NULL);
	
	
	SDFclose(sdfp);
	free(body);
	return 0;
}