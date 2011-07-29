#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
//#include "singlio.h"
#include "Msgs.h"
#include <stdlib.h>
#include <math.h>

#define SPHOUTBODYDESC \
"struct {\n\
	double x, y, z;		/* position of body */\n\
	float mass;			/* mass of body */\n\
	float vx, vy, vz;		/* velocity of body */\n\
	float u;			/* internal energy */\n\
	float h;			/* smoothing length */\n\
	float rho;			/* density */\n\
	float drho_dt;              /* time derivative of rho */\n\
	float udot;			/* time derivative of u */\n\
	float temp;           /* temperature of body */\n\
	float abar;           /* avg number of nucleons per particle of body */\n\
	float zbar;           /* avg number of protons per particle of body */\n\
	float ax, ay, az;		/* acceleration */\n\
	float lax, lay, laz;	/* acceleration at tpos-dt */\n\
	float gax, gay, gaz;  /* gravity acceleration of body */\n\
	float grav_mass;      /* gravitational mass of body */\n\
	float phi;			/* potential */\n\
	float idt;			/* timestep */\n\
	unsigned int nbrs;          /* number of neighbors */\n\
	unsigned int ident;		/* unique identifier */\n\
	unsigned int windid;        /* wind id */\n\
}"


typedef struct {
	double x, y, z;             /* position of body */
	float mass;           /* mass of body */
	float vx, vy, vz;     /* velocity of body */
	float u;              /* specific energy of body*/
	float h;              /* smoothing length of body */
	float rho;            /* density of body */
	//float pr;            /* pressure of body */
	float drho_dt;        /* drho/dt of body */
	float udot;           /* du/dt of body */
	float temp;           /* temperature of body */
	//float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */
	//float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */
	float abar;           /* avg number of nucleons per particle of body */
	float zbar;           /* avg number of protons per particle of body */
	float ax, ay, az;     /* acceleration of body */
	float lax, lay, laz;  /* last acceleration of body */
	float gax, gay, gaz;  /* gravity acceleration of body */
	float grav_mass;      /* gravitational mass of body */
	float phi;            /* potential at body location */
	//float tacc;           /* time of last acceleration update of body */
	float idt;
	unsigned int nbrs;     /* number of neighbors */
	unsigned int ident;    /* unique identifier */
	unsigned int windid;   /* wind id */
	//unsigned int useless;  /* to fill the double block (4int=1double)*/		
} SDFin;

int main(int argc, char **argv[])
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
	float cpos = pow(10.,-3.);
	float cmass = pow(10.,-6.);
	float ct = 1.;
	float crho = cmass / (cpos*cpos*cpos);
	float cu = (cpos*cpos)/(ct*ct);

	char name [80];
	
	if (argc < 2){
		printf("SDF file: ");
		gets (argv[1]);
	}
	singlPrintf("Reading \"%s\"\n", argv[1]);

	SDF *sdfp;
	SDFin *inArray;
	int gnobj, nobj;
	int conf;
	sdfp = SDFreadf(argv[1], (void **)&inArray, &gnobj, &nobj, sizeof(SDFin),
					"x", offsetof(SDFin, x), &conf,
					"y", offsetof(SDFin, y), &conf,
					"z", offsetof(SDFin, z), &conf,
					"mass", offsetof(SDFin, mass), &conf,
					"vx", offsetof(SDFin, vx), &conf,
					"vy", offsetof(SDFin, vy), &conf,
					"vz", offsetof(SDFin, vz), &conf,
					"u", offsetof(SDFin, u), &conf,
					"h", offsetof(SDFin, h), &conf,
					"rho", offsetof(SDFin, rho), &conf,
					//"pr", offsetof(SDFin, pr), &conf,
					"drho_dt", offsetof(SDFin, drho_dt), &conf,
					"udot", offsetof(SDFin, udot), &conf,
					"temp", offsetof(SDFin, temp), &conf,
					//"He4", offsetof(SDFin, He4), &conf,
					//"C12", offsetof(SDFin, C12), &conf,
					//"O16", offsetof(SDFin, O16), &conf,
//					"Ne20", offsetof(SDFin, Ne20), &conf,
//					"Mg24", offsetof(SDFin, Mg24), &conf,
//					"Si28", offsetof(SDFin, Si28), &conf,
//					"S32", offsetof(SDFin, S32), &conf,
//					"Ar36", offsetof(SDFin, Ar36), &conf,
//					"Ca40", offsetof(SDFin, Ca40), &conf,
//					"Ti44", offsetof(SDFin, Ti44), &conf,
//					"Cr48", offsetof(SDFin, Cr48), &conf,
//					"Fe52", offsetof(SDFin, Fe52), &conf,
//					"Ni56", offsetof(SDFin, Ni56), &conf,
					"abar", offsetof(SDFin, abar), &conf,
					"zbar", offsetof(SDFin, zbar), &conf,
					"ax", offsetof(SDFin, ax), &conf,
					"ay", offsetof(SDFin, ay), &conf,
					"az", offsetof(SDFin, az), &conf,
					"lax", offsetof(SDFin, lax), &conf,
					"lay", offsetof(SDFin, lay), &conf,
					"laz", offsetof(SDFin, laz), &conf,
					"gax", offsetof(SDFin, gax), &conf,
					"gay", offsetof(SDFin, gay), &conf,
					"gaz", offsetof(SDFin, gaz), &conf,
					"grav_mass", offsetof(SDFin, grav_mass), &conf,
					"phi", offsetof(SDFin, phi), &conf,
					//"tacc", offsetof(SDFin, tacc), &conf,
					"idt", offsetof(SDFin, idt), &conf,
					"nbrs", offsetof(SDFin, nbrs), &conf,
					"ident", offsetof(SDFin, ident), &conf,
					"windid", offsetof(SDFin, windid), &conf,
					//"useless", offsetof(SDFin, useless), &conf,
					NULL);
	//now to fill the misc. floats and ints//
	SDFgetfloatOrDefault(sdfp, "dt",  &dt, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "eps",  &eps, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "Gnewt",  &Gnewt, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "tolerance",  &tolerance, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "frac_tolerance",  &frac_tolerance, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "tvel",  &tvel, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "gamma",  &gamma, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "ke",  &ke, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "pe",  &pe, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "te",  &te, (float)0.0);
	SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);

	
	
	//Msgf(("Data read, SPHnobj=%d, SPHgnobj=%d\n", nobj, gnobj));
	//SDFout *outArray = malloc(sizeof(SDFout) * nobj);
	int i;
	for (i=0; i<nobj; i+=1) {
		inArray[i].x = inArray[i].x / cpos;
		inArray[i].y = inArray[i].y / cpos;
		inArray[i].z = inArray[i].z / cpos;
		inArray[i].h = inArray[i].h / cpos;
		inArray[i].mass = inArray[i].mass / cmass;
		inArray[i].rho = inArray[i].rho / crho;
		inArray[i].u = inArray[i].u / cu;
	}
	char outname [80];
	snprintf(outname, sizeof(outname), "%s.units", argv[1]);
	SDFwrite(outname, gnobj, 
			 nobj, inArray, sizeof(SDFin),
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
	//free(outArray);
	return 0;
}