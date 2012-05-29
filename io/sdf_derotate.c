
#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
#include "Msgs.h"
#include <stdlib.h>
#include <omp.h>


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
char infile[80];
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

int iter,iter0,iterN,diter,citer;	

float tposold,tpos,deltat;

float omega,thetaold,theta,tmass,oi,x,y,z,vx,vy,vz;

float cm[3],cv[3];


int Error()
{
	
}

int main()
{

	
	char rootin[80] = "interp_sph"; //can be replaced with a gets later
	char rootout[80] = "derot_sph";
	
	printf("iter0 iterN diter:");
	scanf("%d %d %d",&iter0,&iterN,&diter);
	
	thetaold	= 0;
	tposold		= 0;
 
	for(iter=iter0;iter<iterN+diter;iter+=diter)
	{
		SDF *sdfp;
		SPHbody *body1;
		
		snprintf(infile, sizeof(infile), "%s.%04d", rootin,iter);
		
		printf("Opening %s\n",infile);
		
		sdfp = SDFreadf(infile, (void **)&body1, &gnobj, &nobj, sizeof(SPHbody),
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
						//"kappa", offsetof(SPHbody, kappa), &conf,
						//"sigma", offsetof(SPHbody, sigma), &conf,
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
		
		omega = 0.0;
		tmass = 0.0;
		
		printf("Recentering.\n");
		
		for (i=0; i<nobj; i++) 
		{
			cm[0] += body1[i].mass * body1[i].x;
			cm[1] += body1[i].mass * body1[i].y;
			cm[2] += body1[i].mass * body1[i].z;
			cv[0] += body1[i].mass * body1[i].vx;
			cv[1] += body1[i].mass * body1[i].vy;
			cv[2] += body1[i].mass * body1[i].vz;
			tmass += body1[i].mass;
		}
		
		for (i=0;i<3;i++) 
		{
			cm[i] = cm[i]/tmass;
			cv[i] = cv[i]/tmass;
		}
		
		
#pragma omp parallel for default(none) \
	shared(body1,cm,cv,nobj) private(i)
		for (i=0; i<nobj; i++) 
		{
			body1[i].x -= cm[0];
			body1[i].y -= cm[1];
			body1[i].z -= cm[2];
			body1[i].vx -= cv[0];
			body1[i].vy -= cv[1];
			body1[i].vz -= cv[2];		
		}
		
		printf("Derotating.\n");
		
		for (i=0; i<nobj; i++) 
		{
			x = body1[i].x;
			y = body1[i].y;
			vx = body1[i].vx;
			vy = body1[i].vy;
			oi = (vx*y - vy*x)/(x*x + y*y);
			omega += body1[i].mass * oi;
		}

		omega = omega/tmass;
		theta = thetaold + omega*(tpos-tposold);

#pragma omp parallel for default(none) \
	shared(body1,theta,nobj) private(i,x,y)
		for (i=0; i<nobj; i++) 
		{
			x = body1[i].x;
			y = body1[i].y;
			body1[i].x = x*cos(theta) - y*sin(theta);
			body1[i].y = x*sin(theta) + y*cos(theta);
		}
		
		thetaold	= theta;
		tposold		= tpos;
								 
		snprintf(outfile, sizeof(outfile), "%s.%04d", rootout,iter);
		
		printf("Writing %s.\n",outfile);
		
		SDFwrite(outfile, gnobj, 
				 nobj, body1, sizeof(SPHbody),
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
		

		free(body1);
		SDFclose(sdfp);
		
				
	}
	
	return 0;
}
