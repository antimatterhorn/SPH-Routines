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
	//float gax, gay, gaz;  /* gravity acceleration of body */
	//float grav_mass;      /* gravitational mass of body */
	float phi;            /* potential at body location */
	float tacc;           /* time of last acceleration update of body */
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
	
	int i,j;

	char name [80];
	printf("filename: ");
	gets (name);
	singlPrintf("Reading \"%s\"\n", name);

	SDF *sdfp;
	SPHbody *inArray;
	int gnobj, nobj;
	int conf;
	sdfp = SDFreadf(name, (void **)&inArray, &gnobj, &nobj, sizeof(SPHbody),
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
					//"pr", offsetof(SPHbody, pr), &conf,
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
//					"gax", offsetof(SPHbody, gax), &conf,
//					"gay", offsetof(SPHbody, gay), &conf,
//					"gaz", offsetof(SPHbody, gaz), &conf,
//					"grav_mass", offsetof(SPHbody, grav_mass), &conf,
					"phi", offsetof(SPHbody, phi), &conf,
					"tacc", offsetof(SPHbody, tacc), &conf,
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

	
	
	//Msgf(("Data read, SPHnobj=%d, SPHgnobj=%d\n", nobj, gnobj));
	//SDFout *outArray = malloc(sizeof(SDFout) * nobj);

	printf("display how many particles:");
	scanf("%d",&j);
	
	for (i=0; i<j; i+=1) {
		singlPrintf("x = %f\n", inArray[i].x);
		singlPrintf("y = %f\n", inArray[i].y);
		singlPrintf("z = %f\n", inArray[i].z);
		singlPrintf("mass = %f\n", inArray[i].mass);
		singlPrintf("vx = %f\n", inArray[i].vx);
		singlPrintf("vy = %f\n", inArray[i].vy);
		singlPrintf("vz = %f\n", inArray[i].vz);
		singlPrintf("u = %f\n", inArray[i].u);
		singlPrintf("h = %f\n", inArray[i].h);
		singlPrintf("rho = %f\n", inArray[i].rho);
//		singlPrintf("pr = %f\n", inArray[i].pr);
		singlPrintf("drho_dt = %f\n", inArray[i].drho_dt);
		singlPrintf("udot = %f\n", inArray[i].udot);
		singlPrintf("temp = %f\n", inArray[i].temp);
		singlPrintf("abar = %f\n", inArray[i].abar);
		singlPrintf("zbar = %f\n", inArray[i].zbar);
		singlPrintf("ax = %f\n", inArray[i].ax);
		singlPrintf("ay = %f\n", inArray[i].ay);
		singlPrintf("az = %f\n", inArray[i].az);
		singlPrintf("lax = %f\n", inArray[i].lax);
		singlPrintf("lay = %f\n", inArray[i].lay);
		singlPrintf("laz = %f\n", inArray[i].laz);
//		singlPrintf("gax = %f\n", inArray[i].gax);
//		singlPrintf("gay = %f\n", inArray[i].gay);
//		singlPrintf("gaz = %f\n", inArray[i].gaz);
//		singlPrintf("grav_mass = %f\n", inArray[i].grav_mass);
		singlPrintf("phi = %f\n", inArray[i].phi);
		singlPrintf("tacc = %f\n", inArray[i].tacc);
		singlPrintf("idt = %f\n", inArray[i].idt);	
	}

	return 0;
}