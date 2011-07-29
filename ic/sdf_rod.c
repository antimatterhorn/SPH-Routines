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
	int conf,i,j,k,n,id;
	int nx,ny,nz;
	
	float x,y,z;
	
	float l,w,v,rho,m,h,u,abar,zbar;
	l = 20;
	w = 5;
	v = w*w*l;
	rho = 100;
	

	char outfile[80];
	
	n = atoi(argv[1]);
	h = pow(v/(4/3*3.14159*n),0.3333);
	m = rho * (4/3)*(3.14159)*h*h*h;
	
	nz = w/(2*h);
	ny = nz;
	nx = l/(2*h);
	
	printf("h = %3.2e, total mass = %3.2e\n",h,m*n);
	printf("%d x %d x %d = %d\n",nz,ny,nx,nx*ny*nz);
	nobj = nx*ny*nz;
		
	abar = 13.71;
	zbar = 6.857;
	u = 4.8;
		
	SPHbody *body = malloc(sizeof(SPHbody) * (nobj));
	
	for(k=0;k<nz;k++)
	{
		z = k*2*h-w/2;
		for(j=0;j<ny;j++)
		{
			y = j*2*h-w/2;
			for(i=0;i<nx;i++)
			{
				x = i*2*h-h;
				id = i + j*nx + k*nx*ny;
				
				body[id].ident	= id+1;
				body[id].x		= x;
				if(i % 2 == 0){
					body[id].y	= y;
					body[id].z	= z;
				} else {
					body[id].y	= y-h;
					body[id].z	= z-h;
				}
				
				
				body[id].h		= h;
				body[id].mass	= m;
				body[id].rho	= rho;
				body[id].u		= u;
				if(x < 1.0) body[id].u = 5.5;
				body[id].abar	= abar;
				body[id].zbar	= zbar;
				body[id].pr		= 0;
				body[id].nbrs	= 40;
				body[id].temp	= 1e7;
			}
		}
	}
	
	
	gnobj = nobj;
	
	
	
	//write the output file
	snprintf(outfile, sizeof(outfile), "rod.sdf");
	
	SDFwrite(outfile, gnobj, 
			 nobj, body, sizeof(SPHbody),
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


