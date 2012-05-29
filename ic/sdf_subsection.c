#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
#include "Msgs.h"
#include <stdlib.h>


/* This code will cut out a subsection of particles in a specified sphere	*/

#define SPHOUTBODYDESC \
"struct {\n\
	double x, y, z;		/* position of body */\n\
	float mass;			/* mass of body */\n\
	float vx, vy, vz;		/* velocity of body */\n\
	float u;			/* internal energy */\n\
	float h;			/* smoothing length */\n\
	float rho;			/* density */\n\
	float pr;			/* pressure */\n\
	float drho_dt;              /* time derivative of rho */\n\
	float udot;			/* time derivative of u */\n\
	float temp;           /* temperature of body */\n\
	float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */\n\
	float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */\n\
	float vsound;		/* sound speed */\n\
	float abar;           /* avg number of nucleons per particle of body */\n\
	float zbar;           /* avg number of protons per particle of body */\n\
	float ax, ay, az;		/* acceleration */\n\
	float lax, lay, laz;	/* acceleration at tpos-dt */\n\
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
	float pr;            /* pressure of body */
	float drho_dt;        /* drho/dt of body */
	float udot;           /* du/dt of body */
	float temp;           /* temperature of body */
	float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */
	float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */
	float vsound;			/* sound speed */
	float abar;           /* avg number of nucleons per particle of body */
	float zbar;           /* avg number of protons per particle of body */
	float ax, ay, az;     /* acceleration of body */
	float lax, lay, laz;  /* last acceleration of body */
	//float gax, gay, gaz;  /* gravity acceleration of body */
	//float grav_mass;      /* gravitational mass of body */
	float phi;            /* potential at body location */
	//float tacc;           /* time of last acceleration update of body */
	float idt;
	unsigned int nbrs;     /* number of neighbors */
	unsigned int ident;    /* unique identifier */
	unsigned int windid;   /* wind id */
	//unsigned int useless;  /* to fill the double block (4int=1double)*/		
} SPHbody;

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
double tmass;

double x1,x2,y1,y2,z1,z2;
double x,y,z;
double rbox;

int gnobj1, nobj1, gnobj2,nobj2;
int conf,i,j,k;



char outfile[80];

int usage()
{
	printf("\t Cuts out a subsample of particles in a specified sphere (rbox)\n");
	printf("\t Usage: [required] {optional}\n");
	//printf("\t sdf_composition [sdf file] [x1] [x2] [y1] [y2] [z1] [z2]\n");
	printf("\t sdf_composition [sdf file] [rbox]\n");
	return 0;
}

int inbox()
{
	//if(x<x2 && x>x1 && y<y2 && y>y1 && z<z2 && z>z1) return 1;
	if(sqrt(x*x+y*y+z*z)<rbox) return 1;
	else return 0;
}

int main(int argc, char **argv[])
{
	
	if (argc < 3){
		usage();
		return 0;
	}
	
	rbox = atof(argv[2]);
	
	/*
	x1 = atoi(argv[2]);
	x2 = atoi(argv[3]);
	y1 = atoi(argv[4]);
	y2 = atoi(argv[5]);
	z1 = atoi(argv[6]);
	z2 = atoi(argv[7]);
	 */
	
	SDF *sdfp;
	SPHbody *body1;

	sdfp = SDFreadf(argv[1], (void **)&body1, &gnobj1, &nobj1, sizeof(SPHbody),
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
					"pr", offsetof(SPHbody, pr), &conf,
					"drho_dt", offsetof(SPHbody, drho_dt), &conf,
					"udot", offsetof(SPHbody, udot), &conf,
					"temp", offsetof(SPHbody, temp), &conf,
					"He4", offsetof(SPHbody, He4), &conf,
					"C12", offsetof(SPHbody, C12), &conf,
					"O16", offsetof(SPHbody, O16), &conf,
					"Ne20", offsetof(SPHbody, Ne20), &conf,
					"Mg24", offsetof(SPHbody, Mg24), &conf,
					"Si28", offsetof(SPHbody, Si28), &conf,
					"S32", offsetof(SPHbody, S32), &conf,
					"Ar36", offsetof(SPHbody, Ar36), &conf,
					"Ca40", offsetof(SPHbody, Ca40), &conf,
					"Ti44", offsetof(SPHbody, Ti44), &conf,
					"Cr48", offsetof(SPHbody, Cr48), &conf,
					"Fe52", offsetof(SPHbody, Fe52), &conf,
					"Ni56", offsetof(SPHbody, Ni56), &conf,
					"vsound", offsetof(SPHbody, vsound), &conf,
					"abar", offsetof(SPHbody, abar), &conf,
					"zbar", offsetof(SPHbody, zbar), &conf,
					"ax", offsetof(SPHbody, ax), &conf,
					"ay", offsetof(SPHbody, ay), &conf,
					"az", offsetof(SPHbody, az), &conf,
					"lax", offsetof(SPHbody, lax), &conf,
					"lay", offsetof(SPHbody, lay), &conf,
					"laz", offsetof(SPHbody, laz), &conf,
					//"gax", offsetof(SPHbody, gax), &conf,
					//"gay", offsetof(SPHbody, gay), &conf,
					//"gaz", offsetof(SPHbody, gaz), &conf,
					//"grav_mass", offsetof(SPHbody, grav_mass), &conf,
					"phi", offsetof(SPHbody, phi), &conf,
					//"tacc", offsetof(SPHbody, tacc), &conf,
					"idt", offsetof(SPHbody, idt), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	//now to fill the misc. floats and ints//
	SDFgetfloatOrDefault(sdfp, "dt",  &dt, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "eps",  &eps, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "Gnewt",  &Gnewt, (float)0.000393935);
	SDFgetfloatOrDefault(sdfp, "tolerance",  &tolerance, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "frac_tolerance",  &frac_tolerance, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "tvel",  &tvel, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "gamma",  &gamma, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "ke",  &ke, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "pe",  &pe, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "te",  &te, (float)0.0);
	SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);


	nobj2 = 0;
	
	for(i=0;i<nobj1;i++)
	{
		x = body1[i].x;
		y = body1[i].y;
		z = body1[i].z;
		if(inbox()) nobj2++;
	} 
	printf("found %d particles in the box\n",nobj2);
	gnobj2 = nobj2;
	SPHbody *outbody = malloc(sizeof(SPHbody) * (nobj2));

	j=0;
	for(i=0;i<nobj1;i++)
	{
		x = body1[i].x;
		y = body1[i].y;
		z = body1[i].z;
		if(inbox())
		{
			outbody[j].x = body1[i].x;
			outbody[j].y = body1[i].y;
			outbody[j].z = body1[i].z;
			outbody[j].mass = body1[i].mass;
			outbody[j].vx = body1[i].vx;
			outbody[j].vy = body1[i].vy;
			outbody[j].vz = body1[i].vz;
			outbody[j].u = body1[i].u;
			outbody[j].h = body1[i].h;
			outbody[j].rho = body1[i].rho;
			outbody[j].pr = body1[i].pr;
			outbody[j].drho_dt = body1[i].drho_dt;
			outbody[j].udot = body1[i].udot;
			outbody[j].temp = body1[i].temp;
			outbody[j].He4 = body1[i].He4;
			outbody[j].C12 = body1[i].C12;
			outbody[j].O16 = body1[i].O16;
			outbody[j].Ne20 = body1[i].Ne20;
			outbody[j].Mg24 = body1[i].Mg24;
			outbody[j].Si28 = body1[i].Si28;
			outbody[j].S32 = body1[i].S32;
			outbody[j].Ar36 = body1[i].Ar36;
			outbody[j].Ca40 = body1[i].Ca40;
			outbody[j].Ti44 = body1[i].Ti44;
			outbody[j].Cr48 = body1[i].Cr48;
			outbody[j].Fe52 = body1[i].Fe52;
			outbody[j].Ni56 = body1[i].Ni56;
			outbody[j].vsound = body1[i].vsound;
			outbody[j].abar = body1[i].abar;
			outbody[j].zbar = body1[i].zbar;
			outbody[j].ax = body1[i].ax;
			outbody[j].ay = body1[i].ay;
			outbody[j].az = body1[i].az;
			outbody[j].lax = body1[i].lax;
			outbody[j].lay = body1[i].lay;
			outbody[j].laz = body1[i].laz;
			//		outbody[j].gax = body1[i].gax;
			//		outbody[j].gay = body1[i].gay;
			//		outbody[j].gaz = body1[i].gaz;
			//outbody[j].grav_mass = body1[i].grav_mass;
			outbody[j].phi = body1[i].phi;
			//outbody[j].tacc = body1[i].tacc;
			outbody[j].idt = body1[i].idt;
			outbody[j].nbrs = body1[i].nbrs;
			outbody[j].ident = body1[i].ident;
			outbody[j].windid = body1[i].windid;
			//outbody[j].useless = body1[i].useless;
			j++;
		}
	}
	//write the output file
	snprintf(outfile, sizeof(outfile), "sub_%s", argv[1]);
	//*fp = fopen(outfile, "w");
	
	SDFwrite(outfile, gnobj2, 
			 nobj2, outbody, sizeof(SPHbody),
			 SPHOUTBODYDESC,
			 "npart", SDF_INT, gnobj2,
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
	free(outbody);
	
	free(body1);

	
	return 0;
}


