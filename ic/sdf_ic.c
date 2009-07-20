#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
#include "Msgs.h"
#include <stdlib.h>

/* This code will setup the initial conditions SDF file		*/
/* for a two-body system, given positions and velocities	*/
/* of each.													*/

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
	float pr;            /* pressure of body */
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
} SPHbody;


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
	
	double pos1[3],pos2[3];
	float vel1[3],vel2[3];
	double dummy1,dummy2,dummy3;
	
	SDF *sdfp;
	SDF *sdfp2;
	SPHbody *body1;
	SPHbody *body2;
	int gnobj1, nobj1;
	int gnobj2, nobj2;
	int conf,i;
	
	char file1[80];
	char file2[80];
	char outfile[80];
	
	printf("file1: ");
	gets (file1);
	
	printf("file2: ");
	gets (file2);
	
	printf("x1 y1 z1: ");
	scanf("%lf %lf %lf", &dummy1, &dummy2, &dummy3);
	pos1[0] = dummy1;
	pos1[1] = dummy2;
	pos1[2] = dummy3;
	
	printf("vx1 vy1 vz1: ");
	scanf("%lf %lf %lf", &dummy1, &dummy2, &dummy3);
	vel1[0] = dummy1;
	vel1[1] = dummy2;
	vel1[2] = dummy3;
	
	printf("x2 y2 z2: ");
	scanf("%lf %lf %lf", &dummy1, &dummy2, &dummy3);	
	pos2[0] = dummy1;
	pos2[1] = dummy2;
	pos2[2] = dummy3;
	
	printf("vx2 vy2 vz2: ");
	scanf("%lf %lf %lf", &dummy1, &dummy2, &dummy3);	
	vel2[0] = dummy1;
	vel2[1] = dummy2;
	vel2[2] = dummy3;
	
	//read first file
	sdfp = SDFreadf(file1, (void **)&body1, &gnobj1, &nobj1, sizeof(SPHbody),
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
					//"He4", offsetof(SPHbody, He4), &conf,
					//"C12", offsetof(SPHbody, C12), &conf,
					//"O16", offsetof(SPHbody, O16), &conf,
					//					"Ne20", offsetof(SPHbody, Ne20), &conf,
					//					"Mg24", offsetof(SPHbody, Mg24), &conf,
					//					"Si28", offsetof(SPHbody, Si28), &conf,
					//					"S32", offsetof(SPHbody, S32), &conf,
					//					"Ar36", offsetof(SPHbody, Ar36), &conf,
					//					"Ca40", offsetof(SPHbody, Ca40), &conf,
					//					"Ti44", offsetof(SPHbody, Ti44), &conf,
					//					"Cr48", offsetof(SPHbody, Cr48), &conf,
					//					"Fe52", offsetof(SPHbody, Fe52), &conf,
					//					"Ni56", offsetof(SPHbody, Ni56), &conf,
					"abar", offsetof(SPHbody, abar), &conf,
					"zbar", offsetof(SPHbody, zbar), &conf,
					"ax", offsetof(SPHbody, ax), &conf,
					"ay", offsetof(SPHbody, ay), &conf,
					"az", offsetof(SPHbody, az), &conf,
					"lax", offsetof(SPHbody, lax), &conf,
					"lay", offsetof(SPHbody, lay), &conf,
					"laz", offsetof(SPHbody, laz), &conf,
					"gax", offsetof(SPHbody, gax), &conf,
					"gay", offsetof(SPHbody, gay), &conf,
					"gaz", offsetof(SPHbody, gaz), &conf,
					"grav_mass", offsetof(SPHbody, grav_mass), &conf,
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
	
	//read second file
	sdfp2 = SDFreadf(file2, (void **)&body2, &gnobj2, &nobj2, sizeof(SPHbody),
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
					//"He4", offsetof(SPHbody, He4), &conf,
					//"C12", offsetof(SPHbody, C12), &conf,
					//"O16", offsetof(SPHbody, O16), &conf,
					//					"Ne20", offsetof(SPHbody, Ne20), &conf,
					//					"Mg24", offsetof(SPHbody, Mg24), &conf,
					//					"Si28", offsetof(SPHbody, Si28), &conf,
					//					"S32", offsetof(SPHbody, S32), &conf,
					//					"Ar36", offsetof(SPHbody, Ar36), &conf,
					//					"Ca40", offsetof(SPHbody, Ca40), &conf,
					//					"Ti44", offsetof(SPHbody, Ti44), &conf,
					//					"Cr48", offsetof(SPHbody, Cr48), &conf,
					//					"Fe52", offsetof(SPHbody, Fe52), &conf,
					//					"Ni56", offsetof(SPHbody, Ni56), &conf,
					"abar", offsetof(SPHbody, abar), &conf,
					"zbar", offsetof(SPHbody, zbar), &conf,
					"ax", offsetof(SPHbody, ax), &conf,
					"ay", offsetof(SPHbody, ay), &conf,
					"az", offsetof(SPHbody, az), &conf,
					"lax", offsetof(SPHbody, lax), &conf,
					"lay", offsetof(SPHbody, lay), &conf,
					"laz", offsetof(SPHbody, laz), &conf,
					"gax", offsetof(SPHbody, gax), &conf,
					"gay", offsetof(SPHbody, gay), &conf,
					"gaz", offsetof(SPHbody, gaz), &conf,
					"grav_mass", offsetof(SPHbody, grav_mass), &conf,
					"phi", offsetof(SPHbody, phi), &conf,
					//"tacc", offsetof(SPHbody, tacc), &conf,
					"idt", offsetof(SPHbody, idt), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	
	//create and fill output structure
	SPHbody *outbody = malloc(sizeof(SPHbody) * (nobj1+nobj2));

	for (i=0; i<nobj1; i++) {
		outbody[i].x = body1[i].x + pos1[0];
		outbody[i].y = body1[i].y + pos1[1];
		outbody[i].z = body1[i].z + pos1[2];
		outbody[i].mass = body1[i].mass;
		outbody[i].vx = body1[i].vx + vel1[0];
		outbody[i].vy = body1[i].vy + vel1[1];
		outbody[i].vz = body1[i].vz + vel1[2];
		outbody[i].u = body1[i].u;
		outbody[i].h = body1[i].h;
		outbody[i].rho = body1[i].rho;
		outbody[i].pr = body1[i].pr;
		outbody[i].drho_dt = body1[i].drho_dt;
		outbody[i].udot = body1[i].udot;
		outbody[i].temp = body1[i].temp;
		//		outbody[i].He4 = body1[i].He4;
		//		outbody[i].C12 = body1[i].C12;
		//		outbody[i].O16 = body1[i].O16;
		//		outbody[i].Ne20 = body1[i].Ne20;
		//		outbody[i].Mg24 = body1[i].Mg24;
		//		outbody[i].Si28 = body1[i].Si28;
		//		outbody[i].S32 = body1[i].S32;
		//		outbody[i].Ar36 = body1[i].Ar36;
		//		outbody[i].Ca40 = body1[i].Ca40;
		//		outbody[i].Ti44 = body1[i].Ti44;
		//		outbody[i].Cr48 = body1[i].Cr48;
		//		outbody[i].Fe52 = body1[i].Fe52;
		//		outbody[i].Ni56 = body1[i].Ni56;
		outbody[i].abar = body1[i].abar;
		outbody[i].zbar = body1[i].zbar;
		outbody[i].ax = body1[i].ax;
		outbody[i].ay = body1[i].ay;
		outbody[i].az = body1[i].az;
		outbody[i].lax = body1[i].lax;
		outbody[i].lay = body1[i].lay;
		outbody[i].laz = body1[i].laz;
		outbody[i].gax = body1[i].gax;
		outbody[i].gay = body1[i].gay;
		outbody[i].gaz = body1[i].gaz;
		outbody[i].grav_mass = body1[i].grav_mass;
		outbody[i].phi = body1[i].phi;
		//outbody[i].tacc = body1[i].tacc;
		outbody[i].idt = body1[i].idt;
		outbody[i].nbrs = body1[i].nbrs;
		outbody[i].ident = i;
		outbody[i].windid = body1[i].windid;
		//outbody[i].useless = body1[i].useless;
	}
	for (i=0; i<nobj2; i++) {
		outbody[i+nobj1].x = body2[i].x + pos2[0];
		outbody[i+nobj1].y = body2[i].y + pos2[1];
		outbody[i+nobj1].z = body2[i].z + pos2[2];
		outbody[i+nobj1].mass = body2[i].mass;
		outbody[i+nobj1].vx = body2[i].vx + vel2[0];
		outbody[i+nobj1].vy = body2[i].vy + vel2[1];
		outbody[i+nobj1].vz = body2[i].vz + vel2[2];
		outbody[i+nobj1].u = body2[i].u;
		outbody[i+nobj1].h = body2[i].h;
		outbody[i+nobj1].rho = body2[i].rho;
		outbody[i+nobj1].pr = body2[i].pr;
		outbody[i+nobj1].drho_dt = body2[i].drho_dt;
		outbody[i+nobj1].udot = body2[i].udot;
		outbody[i+nobj1].temp = body2[i].temp;
		//		outbody[i+nobj1].He4 = body2[i].He4;
		//		outbody[i+nobj1].C12 = body2[i].C12;
		//		outbody[i+nobj1].O16 = body2[i].O16;
		//		outbody[i+nobj1].Ne20 = body2[i].Ne20;
		//		outbody[i+nobj1].Mg24 = body2[i].Mg24;
		//		outbody[i+nobj1].Si28 = body2[i].Si28;
		//		outbody[i+nobj1].S32 = body2[i].S32;
		//		outbody[i+nobj1].Ar36 = body2[i].Ar36;
		//		outbody[i+nobj1].Ca40 = body2[i].Ca40;
		//		outbody[i+nobj1].Ti44 = body2[i].Ti44;
		//		outbody[i+nobj1].Cr48 = body2[i].Cr48;
		//		outbody[i+nobj1].Fe52 = body2[i].Fe52;
		//		outbody[i+nobj1].Ni56 = body2[i].Ni56;
		outbody[i+nobj1].abar = body2[i].abar;
		outbody[i+nobj1].zbar = body2[i].zbar;
		outbody[i+nobj1].ax = body2[i].ax;
		outbody[i+nobj1].ay = body2[i].ay;
		outbody[i+nobj1].az = body2[i].az;
		outbody[i+nobj1].lax = body2[i].lax;
		outbody[i+nobj1].lay = body2[i].lay;
		outbody[i+nobj1].laz = body2[i].laz;
		outbody[i+nobj1].gax = body2[i].gax;
		outbody[i+nobj1].gay = body2[i].gay;
		outbody[i+nobj1].gaz = body2[i].gaz;
		outbody[i+nobj1].grav_mass = body2[i].grav_mass;
		outbody[i+nobj1].phi = body2[i].phi;
		//outbody[i+nobj1].tacc = body2[i].tacc;
		outbody[i+nobj1].idt = body2[i].idt;
		outbody[i+nobj1].nbrs = body2[i].nbrs;
		outbody[i+nobj1].ident = i+nobj1;
		outbody[i+nobj1].windid = body2[i].windid;
		//outbody[i].useless = body1[i].useless;
	}
	
	//write the output file
	snprintf(outfile, sizeof(outfile), "ic_%s", file1);
	SDFwrite(outfile, gnobj1+gnobj2, 
			 nobj1+nobj2, outbody, sizeof(SPHbody),
			 SPHOUTBODYDESC,
			 "npart", SDF_INT, gnobj1+gnobj2,
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
	
	return 0;
}


