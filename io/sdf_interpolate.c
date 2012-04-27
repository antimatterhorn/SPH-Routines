
#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
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
	float bdot;			/* energy generation rate */\n\
	unsigned int nbrs;          /* number of neighbors */\n\
	unsigned int ident;		/* unique identifier */\n\
	unsigned int windid;        /* wind id */\n\
	unsigned int useless;		/* to fill double block */\n\
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
	float bdot;
	unsigned int nbrs;     /* number of neighbors */
	unsigned int ident;    /* unique identifier */
	unsigned int windid;   /* wind id */
	unsigned int useless;  /* to fill the double block (4int=1double)*/		
} SPHbody;

int gnobj, nobj;
int conf, i,j,k;
char filename[80];
char outfile[80];
float dt = 0;
float eps = 0;
float Gnewt = 0;
float tolerance = 0;
float frac_tolerance = 0;
int ndim = 3;
float tvel = 0;
float gamma = 0;
double ke = 0;
double pe = 0;
double te = 0;

float tstart,tend,deltat,tpos;
int iter,diter,steps,citer;	

float tpos0,tpos,tpos1;


float interpfloat(float val0, float val1, float t0, float t1, float ti)
{
	return (val1-val0)/(t1-t0)*(ti-t0) + val0;
}

double interpdouble(double val0, double val1, float t0, float t1, float ti)
{
	return (val1-val0)/((double)t1-(double)t0)*((double)ti-(double)t0) + val0;
}

int Error()
{
	printf("File open failed.\nRequested %d, got %d.\n",iter,citer);
}

int main()
{

	
	char rootin[80] = "fltmass_sph"; //can be replaced with a gets later
	char rootout[80] = "interp_sph";
	
	printf("tstart tend dt:");
	scanf("%f %f %f",&tstart,&tend,&deltat);
	
	steps = (tend-tstart)/deltat;
	printf("%d frames.\n",steps);
	
	printf("iter0:");
	scanf("%d",&iter);
	printf("diter:");
	scanf("%d",&diter);
	
	int thread;
	
 
	for(i=0;i<steps;i++)
	{

		SDF *sdfp0;
		SDF *sdfp1;
		
		tpos0 = 0;
		tpos = (float)tstart + (float)i*(float)deltat;
		tpos1 = 0;
		iter = iter - diter;
		j = i+1000;
		
		do {
			iter += diter;
			
			printf("opening %s.%04d.\n",rootin,iter);
			
			snprintf(filename, sizeof(filename), "%s.%04d", rootin, iter);
			
			SPHbody *body1;
			
			sdfp1 = SDFreadf(filename, (void **)&body1, &gnobj, &nobj, sizeof(SPHbody),NULL);
			SDFgetfloatOrDefault(sdfp1, "tpos",  &tpos1, (float)0.0);
			SDFgetintOrDefault(sdfp1, "iter",  &citer, 0);
			
			if(citer!=iter)
			{
				Error();
				return 0;
			}
			
			printf("requesting t = %3.3f, found t = %3.2f\n",tpos,tpos1);
			
			free(body1);
			SDFclose(sdfp1);
			
		} while (tpos1 < tpos);
		
		printf("Requested t=%2.3f, found t=%2.3f in %s.%04d.\n",tpos,tpos1,rootin,iter);
		
		SPHbody *body1;
		SPHbody *body0;
		SPHbody *body = malloc(sizeof(SPHbody) * nobj);
		
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
						 "bdot", offsetof(SPHbody, bdot), &conf,
						 "nbrs", offsetof(SPHbody, nbrs), &conf,
						 "ident", offsetof(SPHbody, ident), &conf,
						 "windid", offsetof(SPHbody, windid), &conf,
						 "useless", offsetof(SPHbody, useless), &conf,
						 NULL);
		
		snprintf(filename, sizeof(filename), "%s.%04d", rootin, iter-diter);
				
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
						 "bdot", offsetof(SPHbody, bdot), &conf,
						 "nbrs", offsetof(SPHbody, nbrs), &conf,
						 "ident", offsetof(SPHbody, ident), &conf,
						 "windid", offsetof(SPHbody, windid), &conf,
						 "useless", offsetof(SPHbody, useless), &conf,						
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

		
		for (k=0; k<nobj; k++) {
			body[k].x		= interpdouble(body0[k].x,body1[k].x,tpos0,tpos1,tpos);
			body[k].y		= interpdouble(body0[k].y,body1[k].y,tpos0,tpos1,tpos);
			body[k].z		= interpdouble(body0[k].z,body1[k].z,tpos0,tpos1,tpos);
			body[k].mass	= interpfloat(body0[k].mass,body1[k].mass,tpos0,tpos1,tpos);
			body[k].vx		= interpfloat(body0[k].vx,body1[k].vx,tpos0,tpos1,tpos);
			body[k].vy		= interpfloat(body0[k].vy,body1[k].vy,tpos0,tpos1,tpos);
			body[k].vz		= interpfloat(body0[k].vz,body1[k].vz,tpos0,tpos1,tpos);
			body[k].u		= interpfloat(body0[k].u,body1[k].u,tpos0,tpos1,tpos);
			body[k].h		= interpfloat(body0[k].h,body1[k].h,tpos0,tpos1,tpos);
			body[k].rho		= interpfloat(body0[k].rho,body1[k].rho,tpos0,tpos1,tpos);
			body[k].pr		= interpfloat(body0[k].pr,body1[k].pr,tpos0,tpos1,tpos);
			body[k].drho_dt = interpfloat(body0[k].drho_dt,body1[k].drho_dt,tpos0,tpos1,tpos);
			body[k].udot	= interpfloat(body0[k].udot,body1[k].udot,tpos0,tpos1,tpos);
			body[k].temp	= interpfloat(body0[k].temp,body1[k].temp,tpos0,tpos1,tpos);
			body[k].He4		= interpfloat(body0[k].He4,body1[k].He4,tpos0,tpos1,tpos);
			body[k].C12		= interpfloat(body0[k].C12,body1[k].C12,tpos0,tpos1,tpos);
			body[k].O16		= interpfloat(body0[k].O16,body1[k].O16,tpos0,tpos1,tpos);
			body[k].Ne20	= interpfloat(body0[k].Ne20,body1[k].Ne20,tpos0,tpos1,tpos);
			body[k].Mg24	= interpfloat(body0[k].Mg24,body1[k].Mg24,tpos0,tpos1,tpos);
			body[k].Si28	= interpfloat(body0[k].Si28,body1[k].Si28,tpos0,tpos1,tpos);
			body[k].S32		= interpfloat(body0[k].S32,body1[k].S32,tpos0,tpos1,tpos);
			body[k].Ar36	= interpfloat(body0[k].Ar36,body1[k].Ar36,tpos0,tpos1,tpos);
			body[k].Ca40	= interpfloat(body0[k].Ca40,body1[k].Ca40,tpos0,tpos1,tpos);
			body[k].Ti44	= interpfloat(body0[k].Ti44,body1[k].Ti44,tpos0,tpos1,tpos);
			body[k].Cr48	= interpfloat(body0[k].Cr48,body1[k].Cr48,tpos0,tpos1,tpos);
			body[k].Fe52	= interpfloat(body0[k].Fe52,body1[k].Fe52,tpos0,tpos1,tpos);
			body[k].Ni56	= interpfloat(body0[k].Ni56,body1[k].Ni56,tpos0,tpos1,tpos);
			body[k].vsound	= interpfloat(body0[k].vsound,body1[k].vsound,tpos0,tpos1,tpos);
			body[k].abar	= interpfloat(body0[k].abar,body1[k].abar,tpos0,tpos1,tpos);
			body[k].zbar	= interpfloat(body0[k].zbar,body1[k].zbar,tpos0,tpos1,tpos);
			body[k].ax		= interpfloat(body0[k].ax,body1[k].ax,tpos0,tpos1,tpos);
			body[k].ay		= interpfloat(body0[k].ay,body1[k].ay,tpos0,tpos1,tpos);
			body[k].az		= interpfloat(body0[k].az,body1[k].az,tpos0,tpos1,tpos);
			body[k].lax		= interpfloat(body0[k].lax,body1[k].lax,tpos0,tpos1,tpos);
			body[k].lay		= interpfloat(body0[k].lay,body1[k].lay,tpos0,tpos1,tpos);
			body[k].laz		= interpfloat(body0[k].laz,body1[k].laz,tpos0,tpos1,tpos);
//			body[k].gax		= interpfloat(body0[k].gax,body1[k].gax,tpos0,tpos1,tpos);
//			body[k].gay		= interpfloat(body0[k].gay,body1[k].gay,tpos0,tpos1,tpos);
//			body[k].gaz		= interpfloat(body0[k].gaz,body1[k].gaz,tpos0,tpos1,tpos);
//			body[k].grav_mass = interpfloat(body0[k].grav_mass,body1[k].grav_mass,tpos0,tpos1,tpos);
			body[k].phi		= interpfloat(body0[k].phi,body1[k].phi,tpos0,tpos1,tpos);
			body[k].idt		= interpfloat(body0[k].idt,body1[k].idt,tpos0,tpos1,tpos);
			body[k].nbrs	= body0[k].nbrs;
			body[k].ident	= body0[k].ident;
			body[k].windid	= body0[k].windid;
			body[k].bdot	= interpfloat(body0[k].bdot,body1[k].bdot,tpos0,tpos1,tpos);
		}
								 
		snprintf(outfile, sizeof(outfile), "%s.%04d", rootout, j);
		
		SDFwrite(outfile, gnobj, 
				 nobj, body, sizeof(SPHbody),
				 SPHOUTBODYDESC,
				 "npart", SDF_INT, gnobj,
				 "iter", SDF_INT, j,
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
						 
		
		free(body);
		free(body0);
		free(body1);
		SDFclose(sdfp1);
		SDFclose(sdfp0);
		
				
	}
	
	return 0;
}
