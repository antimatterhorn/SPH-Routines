
#include <dirent.h> 
#include <stdio.h> 
#include <glob.h>
//#include "SPHbody.h"
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>

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
float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */
float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */
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
	float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */
	float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */
	float abar;           /* avg number of nucleons per particle of body */
	float zbar;           /* avg number of protons per particle of body */
	float ax, ay, az;     /* acceleration of body */
	float lax, lay, laz;  /* last acceleration of body */
	float gax, gay, gaz;  /* gravity acceleration of body */
	float grav_mass;      /* gravitational mass of body */
	float phi;            /* potential at body location */
	float tacc;           /* time of last acceleration update of body */
	float idt;
	unsigned int nbrs;     /* number of neighbors */
	unsigned int ident;    /* unique identifier */
	unsigned int windid;   /* wind id */
	//unsigned int useless;  /* to fill the double block (4int=1double)*/		
} SPHbody;

int gnobj, nobj;
int conf, i;
char filename[80];

float interpvalue(float val0, float val1, float t0, float t1, float ti)
{
	return (val1-val0)/(t1-t0)*(ti-t0) + val0;
}

int sdf_interpolate(float tint, SPHbody *body, int iter0, int iter1, char root)
{
	int iter = 0;
	float dt = 0;
	float eps = 0;
	float Gnewt = 0;
	float tolerance = 0;
	float frac_tolerance = 0;
	int ndim = 3;
	float tpos0 = 0, tpos1 = 0;
	float tvel = 0;
	float gamma = 0;
	double ke = 0;
	double pe = 0;
	double te = 0;
	
	
	
	SDF *sdfp0;
	SDF *sdfp1;
	SPHbody *body0;
	SPHbody *body1;
	snprintf(filename, sizeof(filename), "%s." + iter0, root);
	
	sdfp0 = SDFreadf(filename, (void **)&body0, &gnobj, &nobj, sizeof(SPHbody),
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
					"tacc", offsetof(SPHbody, tacc), &conf,
					"idt", offsetof(SPHbody, idt), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	SDFgetfloatOrDefault(sdfp0, "dt",  &dt, (float)0.0);
	SDFgetfloatOrDefault(sdfp0, "eps",  &eps, (float)0.0);
	SDFgetfloatOrDefault(sdfp0, "Gnewt",  &Gnewt, (float)0.0);
	SDFgetfloatOrDefault(sdfp0, "tolerance",  &tolerance, (float)0.0);
	SDFgetfloatOrDefault(sdfp0, "frac_tolerance",  &frac_tolerance, (float)0.0);
	SDFgetfloatOrDefault(sdfp0, "tpos",  &tpos0, (float)0.0);
	SDFgetfloatOrDefault(sdfp0, "tvel",  &tvel, (float)0.0);
	SDFgetfloatOrDefault(sdfp0, "gamma",  &gamma, (float)0.0);
	SDFgetdoubleOrDefault(sdfp0, "ke",  &ke, (float)0.0);
	SDFgetdoubleOrDefault(sdfp0, "pe",  &pe, (float)0.0);
	SDFgetdoubleOrDefault(sdfp0, "te",  &te, (float)0.0);
	SDFgetintOrDefault  (sdfp0, "iter",  &iter, 0);
	
	snprintf(filename, sizeof(filename), "%s." + iter1, root);
	
	sdfp1 = SDFreadf(filename, (void **)&body1, &gnobj, &nobj, sizeof(SPHbody),
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
					"tacc", offsetof(SPHbody, tacc), &conf,
					"idt", offsetof(SPHbody, idt), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	SDFgetfloatOrDefault(sdfp0, "tpos",  &tpos1, (float)0.0);
	
	body = malloc(sizeof(SPHbody) * nobj);
	
	for (i=0; i<nobj; i++) {
		body[i].x = interpvalue(body0[i].x,body1[i].x,tpos0,tpos1,tint);
		body[i].y = interpvalue(body0[i].y,body1[i].y,tpos0,tpos1,tint);
		body[i].z = interpvalue(body0[i].z,body1[i].z,tpos0,tpos1,tint);
		body[i].mass = interpvalue(body0[i].mass,body1[i].mass,tpos0,tpos1,tint);
		body[i].vx = interpvalue(body0[i].vx,body1[i].vx,tpos0,tpos1,tint);
		body[i].vy = interpvalue(body0[i].vy,body1[i].vy,tpos0,tpos1,tint);
		body[i].vz = interpvalue(body0[i].vz,body1[i].vz,tpos0,tpos1,tint);
		body[i].u = interpvalue(body0[i].u,body1[i].u,tpos0,tpos1,tint);
		body[i].h = interpvalue(body0[i].h,body1[i].h,tpos0,tpos1,tint);
		body[i].rho = interpvalue(body0[i].rho,body1[i].rho,tpos0,tpos1,tint);
		body[i].pr = interpvalue(body0[i].pr,body1[i].pr,tpos0,tpos1,tint);
		body[i].drho_dt = interpvalue(body0[i].drho_dt,body1[i].drho_dt,tpos0,tpos1,tint);
		body[i].udot = interpvalue(body0[i].udot,body1[i].udot,tpos0,tpos1,tint);
		body[i].temp = interpvalue(body0[i].temp,body1[i].temp,tpos0,tpos1,tint);
		body[i].He4 = interpvalue(body0[i].He4,body1[i].He4,tpos0,tpos1,tint);
		body[i].C12 = interpvalue(body0[i].C12,body1[i].C12,tpos0,tpos1,tint);
		body[i].O16 = interpvalue(body0[i].O16,body1[i].O16,tpos0,tpos1,tint);
		body[i].Ne20 = interpvalue(body0[i].Ne20,body1[i].Ne20,tpos0,tpos1,tint);
		body[i].Mg24 = interpvalue(body0[i].Mg24,body1[i].Mg24,tpos0,tpos1,tint);
		body[i].Si28 = interpvalue(body0[i].Si28,body1[i].Si28,tpos0,tpos1,tint);
		body[i].S32 = interpvalue(body0[i].S32,body1[i].S32,tpos0,tpos1,tint);
		body[i].Ar36 = interpvalue(body0[i].Ar36,body1[i].Ar36,tpos0,tpos1,tint);
		body[i].Ca40 = interpvalue(body0[i].Ca40,body1[i].Ca40,tpos0,tpos1,tint);
		body[i].Ti44 = interpvalue(body0[i].Ti44,body1[i].Ti44,tpos0,tpos1,tint);
		body[i].Cr48 = interpvalue(body0[i].Cr48,body1[i].Cr48,tpos0,tpos1,tint);
		body[i].Fe52 = interpvalue(body0[i].Fe52,body1[i].Fe52,tpos0,tpos1,tint);
		body[i].Ni56 = interpvalue(body0[i].Ni56,body1[i].Ni56,tpos0,tpos1,tint);
		body[i].abar = interpvalue(body0[i].abar,body1[i].abar,tpos0,tpos1,tint);
		body[i].zbar = interpvalue(body0[i].zbar,body1[i].zbar,tpos0,tpos1,tint);
		body[i].ax = interpvalue(body0[i].ax,body1[i].ax,tpos0,tpos1,tint);
		body[i].ay = interpvalue(body0[i].ay,body1[i].ay,tpos0,tpos1,tint);
		body[i].az = interpvalue(body0[i].az,body1[i].az,tpos0,tpos1,tint);
		body[i].lax = interpvalue(body0[i].lax,body1[i].lax,tpos0,tpos1,tint);
		body[i].lay = interpvalue(body0[i].lay,body1[i].lay,tpos0,tpos1,tint);
		body[i].laz = interpvalue(body0[i].laz,body1[i].laz,tpos0,tpos1,tint);
		body[i].gax = interpvalue(body0[i].gax,body1[i].gax,tpos0,tpos1,tint);
		body[i].gay = interpvalue(body0[i].gay,body1[i].gay,tpos0,tpos1,tint);
		body[i].gaz = interpvalue(body0[i].gaz,body1[i].gaz,tpos0,tpos1,tint);
		body[i].grav_mass = interpvalue(body0[i].grav_mass,body1[i].grav_mass,tpos0,tpos1,tint);
		body[i].phi = interpvalue(body0[i].phi,body1[i].phi,tpos0,tpos1,tint);
		body[i].idt = interpvalue(body0[i].idt,body1[i].idt,tpos0,tpos1,tint);
		body[i].nbrs = interpvalue(body0[i].nbrs,body1[i].nbrs,tpos0,tpos1,tint);
		body[i].ident = body0[i].ident;
		body[i].windid = body0[i].windid;
		//body[i].useless = body0[i].useless;
	}
	
	return 0;
}

int main(void)
{
	
	return 0;
}
