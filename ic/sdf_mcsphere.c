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
	
	double dist_in_cm;
	double mass_in_g;
	double dens_in_gccm;
	double energy_in_erg;
	double time_in_s;
	double specenergy_in_ergperg;
	double pressure_in_ergperccm;
	
	time_in_s				= 1;
	dist_in_cm				= 6.955e7;
	mass_in_g				= 1.989e27;
	dens_in_gccm			= mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;

	int gnobj, nobj, nobjprev;
	int conf,i,j,k,n;
	int place,fail;
	
	float x,y,z;
	
	float r,rr,rm,rho,m,h,u,abar,zbar;
	rm = 10;
	u = 5;

	char outfile[80];
	
	n = atoi(argv[1]);
	m = 1.0e6/n;
	printf("%d particles of mass %3.2e\n",n,m);	
	h = rm/pow(n,0.3333);
	rho = m/(rm*rm*rm);
	nobj = 1;
	abar = 13.71;
	zbar = 6.857;
	u = 17;
	
	printf("h = %3.2e\n",h);
	
	SPHbody *body = malloc(sizeof(SPHbody) * (n+1));
	
	srand(time(0));
	
	body[0].x		= 0.1;
	body[0].y		= 0.1;
	body[0].z		= 0.1;
	body[0].h		= h/2;
	body[0].mass	= m;
	body[0].rho		= rho;
	body[0].u		= u+3.0;
	body[0].temp	= 1e7;
	body[0].ident	= 1;
	body[0].udot	= 0;
	body[0].abar	= abar;
	body[0].zbar	= zbar;
	body[0].pr		= 0;
	body[0].nbrs	= 40;
	
	
	
	for(i=0;i<n;i++)
	{
		if (nobjprev < nobj && nobj % 1000 == 0) printf("nobj = %d\n",nobj);
		nobjprev = nobj;	
		place = 0;
		for (j=0;j<n*2;j++)
		{
			fail = 0;
			x = rand()/(double)RAND_MAX * 2 * rm - rm;
			y = rand()/(double)RAND_MAX * 2 * rm - rm;
			z = rand()/(double)RAND_MAX * 2 * rm - rm;
			//printf("%f %f %f\n",x,y,z);
			r = sqrt(x*x + y*y + z*z);
			if (r < rm && r>0){
				
				for (k=0;k<nobj;k++)
				{
					rr = pow(fabs(body[k].x - x),2.0) + pow(fabs(body[k].y - y),2.0) + pow(fabs(body[k].z - z),2.0);
					if (rr < h*h) {
						fail = 1;
						break;
					}
				}
				if (fail == 0)
				{
					nobj+=1;
					body[i+1].x		= x;
					body[i+1].y		= y;
					body[i+1].z		= z;
					body[i+1].h		= h/2;
					body[i+1].mass	= m;
					body[i+1].rho		= rho;					
					body[i+1].u		= u;
					if (r < 2.0) body[i+1].u += 3.0;
					body[i+1].temp	= 1e7;
					body[i+1].ident	= nobj;
					body[i+1].udot	= 0;
					body[i+1].abar	= abar;
					body[i+1].zbar	= zbar;
					body[i+1].pr		= 0;
					body[i+1].nbrs	= 40;
				}
			} else fail = 1;
			if (fail == 0) break;
		}
		
	}
	
	gnobj = nobj;
	
	SPHbody *outbody = malloc(sizeof(SPHbody) * nobj);
	
	for(i=0;i<nobj;i++)
	{
		outbody[i].x		= body[i].x;
		outbody[i].y		= body[i].y;
		outbody[i].z		= body[i].z;
		outbody[i].h		= h;
		outbody[i].mass		= m;
		outbody[i].rho		= rho;
		outbody[i].u		= body[i].u;
		outbody[i].temp		= 1e7;
		outbody[i].ident	= body[i].ident;
		outbody[i].udot		= 0;
		outbody[i].abar		= abar;
		outbody[i].zbar		= zbar;
		outbody[i].pr		= 0;
		outbody[i].nbrs		= 20;
	}
	
//	for (i=0; i<nobj; i++) {
//		outbody[i].x = body1[i].x + pos1[0];
//		outbody[i].y = body1[i].y + pos1[1];
//		outbody[i].z = body1[i].z + pos1[2];
//		outbody[i].mass = body1[i].mass;
//		outbody[i].vx = body1[i].vx + vel1[0];
//		outbody[i].vy = body1[i].vy + vel1[1];
//		outbody[i].vz = body1[i].vz + vel1[2];
//		outbody[i].u = body1[i].u;
//		outbody[i].h = body1[i].h;
//	}
	
	
	//write the output file
	snprintf(outfile, sizeof(outfile), "mcsphere.sdf");
	
	SDFwrite(outfile, gnobj, 
			 nobj, outbody, sizeof(SPHbody),
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
	free(body);
	
	return 0;
}


