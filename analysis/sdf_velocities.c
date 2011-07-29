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

double rv(double vx, double vy, double vz)
{
	return sqrt(vx*vx + vy*vy + vz*vz);
}

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
	
	double cm[3],mr[3],vc[3],vr[3],tmass;
	double maxv = 0;
	int bins = 100;
	double iso[bins][14]; //13 isotopes plus 1 velocity column
	double bin,v=0;
	
	double dist_in_cm	= 6.955e7;
	double mass_in_g	= 1.989e27;
	
	SDF *sdfp;
	SPHbody *body;
	int gnobj1, nobj1;
	int conf,i,j;
	
	char file1[80];
	char outfile[80];
	
	if (argc < 2){
		printf("SDF file: ");
		gets (argv[1]);
	}
	
	//read first file
	sdfp = SDFreadf(argv[1], (void **)&body, &gnobj1, &nobj1, sizeof(SPHbody),
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
	SDFgetfloatOrDefault(sdfp, "Gnewt",  &Gnewt, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "tolerance",  &tolerance, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "frac_tolerance",  &frac_tolerance, (float)0.0);
	//SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "tvel",  &tvel, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "gamma",  &gamma, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "ke",  &ke, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "pe",  &pe, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "te",  &te, (float)0.0);
	//SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);
	
	for (i=0;i<3;i++) mr[i] = vr[i] = cm[i] = vc[i] = 0;
	
	for (i=0; i<nobj1; i++) {
		mr[0] += body[i].mass * body[i].x;
		mr[1] += body[i].mass * body[i].y;
		mr[2] += body[i].mass * body[i].z;
		vr[0] += body[i].mass * body[i].vx;
		vr[1] += body[i].mass * body[i].vy;
		vr[2] += body[i].mass * body[i].vz;
		tmass += body[i].mass;
	}
	
	for (i=0;i<3;i++) {
		cm[i] = mr[i]/tmass;
		vc[i] = vr[i]/tmass;
	}

	for (i=0; i<nobj1; i++) {
		body[i].x -= cm[0];
		body[i].y -= cm[1];
		body[i].z -= cm[2];
		body[i].vx -= vc[0];
		body[i].vy -= vc[1];
		body[i].vz -= vc[2];
		if (rv(body[i].vx,body[i].vy,body[i].vz) > maxv) maxv = rv(body[i].vx,body[i].vy,body[i].vz);
	}
	
	maxv = log10(maxv*dist_in_cm);
	
	bin = (double)maxv / (double)bins;
	
	for (i=0;i<bins;i++) {
		v = i*bin;
		iso[i][0] = v;
		for (j=1;j<14;j++) iso[i][j] = 0;
	}
	
	for (i=0; i<nobj1; i++) {
		v = log10(rv(body[i].vx,body[i].vy,body[i].vz)*dist_in_cm);
		for (j=bins-1; j>-1;j--){
			if (v >= iso[j][0]) {
				iso[j][1] += body[i].mass*body[i].He4*mass_in_g/(1.99e33);
				iso[j][2] += body[i].mass*body[i].C12*mass_in_g/(1.99e33);
				iso[j][3] += body[i].mass*body[i].O16*mass_in_g/(1.99e33);
				iso[j][4] += body[i].mass*body[i].Ne20*mass_in_g/(1.99e33);
				iso[j][5] += body[i].mass*body[i].Mg24*mass_in_g/(1.99e33);
				iso[j][6] += body[i].mass*body[i].Si28*mass_in_g/(1.99e33);
				iso[j][7] += body[i].mass*body[i].S32*mass_in_g/(1.99e33);
				iso[j][8] += body[i].mass*body[i].Ar36*mass_in_g/(1.99e33);
				iso[j][9] += body[i].mass*body[i].Ca40*mass_in_g/(1.99e33);
				iso[j][10] += body[i].mass*body[i].Ti44*mass_in_g/(1.99e33);
				iso[j][11] += body[i].mass*body[i].Cr48*mass_in_g/(1.99e33);
				iso[j][12] += body[i].mass*body[i].Fe52*mass_in_g/(1.99e33);
				iso[j][13] += body[i].mass*body[i].Ni56*mass_in_g/(1.99e33);				
				break;
			}
		}
	}
	
	FILE *stream, *fopen();
	/* declare a stream and prototype fopen */ 
	
	stream = fopen("velocities.csv","w");	
	
	fprintf(stream,"logv,He4,C12,O16,Ne20,Mg24,Si28,S32,Ar36,Ca40,Ti44,Cr48,Fe52,Ni56\n");
	
	for (i=0;i<bins;i++) {
		for (j=0;j<13;j++) fprintf(stream,"%f,",iso[i][j]);
		fprintf(stream,"%f\n",iso[i][13]);
	}
	fclose(stream);
	
	free(body);
	
	return 0;
}


