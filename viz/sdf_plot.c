/*
 *  sdf2csv.c
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "sdf_plot.h"
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
float rho,temp,mass,h;
int i;
double r;
char sdffile[80];
char csvfile[80];

int main()
{
	SDF *sdfp;
	SPHbody *body;
	
	printf("SDF file: ");
	gets (sdffile);
		
	sdfp = SDFreadf(sdffile, (void **)&body, &gnobj, &nobj, sizeof(SPHbody),
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
					//"grav_mass", offsetof(SPHbody, grav_mass), &conf,
					"phi", offsetof(SPHbody, phi), &conf,
					"tacc", offsetof(SPHbody, tacc), &conf,
					"idt", offsetof(SPHbody, idt), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	
	singlPrintf("%s has %d particles.\n", sdffile, gnobj);
	
	double **outArray;
	
	outArray = (double **) malloc(nobj*sizeof(double *));
	for(i=0;i<nobj;i++){
		outArray[i]  = (double *) malloc(5*sizeof(double));
	}
	
	
	for(i = 0; i < nobj; i++){
		r = sqrt(pow(body[i].x,2.0)+pow(body[i].y,2.0)+pow(body[i].z,2.0));
		outArray[i][0] = r;
		outArray[i][1] = body[i].rho;
		outArray[i][2] = body[i].temp;
		outArray[i][3] = body[i].mass;
		outArray[i][4] = body[i].u;
	}
	
	//open the stream file
	snprintf(csvfile, sizeof(csvfile), "%splot.csv", sdffile);
	FILE *stream, *fopen();
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");

	fprintf(stream,"r,rho,temp,mass,u\n");
	
	for(i = 0; i < nobj; i++)
	{
		fprintf(stream,"%f,%f,%f,%f,%f\n", outArray[i][0],outArray[i][1],outArray[i][2],
										outArray[i][3],outArray[i][4]);
	}
	fclose(stream);
	
	return 0;
}