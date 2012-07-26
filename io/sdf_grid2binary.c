/*
 *  sdf2csv.c
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
#include "Msgs.h"
#include <stdlib.h>
#include <math.h>

int gnobj, nobj;
int conf;
double xmax,zmax,dist,x,y,z;
double h,rho,temp,mass,kern,r,xc,zc;
int hc,ic,jc;
double time;
float tpos,norm;
int choice,pixels,threads;
int i,j,k,b;
char sdffile[80];
char asciifile[80];

double dist_in_cm;
double mass_in_g;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;
double maxrho,maxtemp;

typedef struct {
	double x, y, z;             /* position of body */							//3
	float mass;           /* mass of body */
	float vx, vy, vz;     /* velocity of body */								//4
	float u;              /* specific energy of body*/
	float h;              /* smoothing length of body */						//6
	float rho;            /* density of body */
	float pr;            /* pressure of body */
	float temp;           /* temperature of body */								//11
	float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */		//18
	float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */			//24
	float vsound;
	float abar;           /* avg number of nucleons per particle of body */
	float zbar;           /* avg number of protons per particle of body */
	float ax, ay, az;     /* acceleration of body */							//30
	float bdot;
	float kappa;			/* opacity */
	float sigma;			/* conductivity */									//38
} SPHbody;

typedef struct {
	double x, z;        /* position of cell */							
	float vx, vy;       /* velocity of cell */								
	float rho;          /* density of cell */
	float temp;         /* temperature of cell */	
    float u;            /* spec energy of cell */
	float He4;
    float C12;
    float O16;
    float Ne20;
    float Mg24;
    float Si28;
    float S32;
    float Ar36;
    float Ca40;
    float Ti44;
    float Cr48;
    float Fe52;
    float Ni56;
} grid;

int usage()
{
	printf("\t Creates a binary output of gridded paticles.\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_grid2binary [sdf file]\n");
	return 0;
}

int main(int argc, char **argv[])
{
	
	time_in_s				= 1;
	dist_in_cm				= 6.955e7;
	mass_in_g				= 1.989e27;
	dens_in_gccm			= mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	
	SDF *sdfp;
	SPHbody *body;
	
	if (argc < 2){
		usage();
		return 0;
	}
		
	sdfp = SDFreadf(argv[1], (void **)&body, &gnobj, &nobj, sizeof(SPHbody),
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
					NULL);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
	
	snprintf(asciifile, sizeof(asciifile), "%s.dat", argv[1]);
	FILE *stream, *fopen();
	stream = fopen(asciifile,"wb");
	fwrite(body,sizeof(SPHbody),nobj,stream);
	
	fclose(stream);
	free(body);

	return 0;
}