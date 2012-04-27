/*
 *  sdf2csv.c
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "SPHbody.h"
#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
#include "Msgs.h"
#include <stdlib.h>
#include <math.h>
#include <bigmalloc.h>

int gnobj, nobj;
int conf;
//int i=0;
//float He4,C12,O16,Ne20,Mg24,Si28,S32;
//float Ar36,Ca40,Ti44,Cr48,Fe52,Ni56;
float tpos;
double x,y,z,rho,pr,temp,u,vsound;
double dist_in_cm;
double mass_in_g;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;
char sdffile[80];
char number[9];
char period[1];
int i,j,start,end,inc;

int main(int argc, char **argv[])
{

	time_in_s				= 1;
	dist_in_cm				= 6.955e7;
	mass_in_g				= 1.989e27;
	dens_in_gccm			= mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	
	
	if (argc < 2){
		printf("on root: ");
		gets (argv[1]);
	}
	
	if (argc < 3){
		printf("Id: ");
		scanf(argv[2]);
	}
	
	strcpy(period,".");
	printf("start end increment: ");
	scanf("%d %d %d",&start,&end,&inc);
	
	FILE *stream, *fopen();
	/* declare a stream and prototype fopen */ 
	
	stream = fopen("out.csv","w");	
	
	fprintf(stream,"x,y,z,t,rho,Pr,T,u,cs,abar,zbar,He4,C12,O16,Ne20,Mg24,Si28,S32,Ar36,Ca40,Ti44,Cr48,Fe52,Ni56,h\n");
	
	for(i=start;i<(end+1);i+=inc)
	{
		strcpy(sdffile,argv[1]);
		strcat(sdffile,period);
		sprintf(number,"%04d",i);
		strcat(sdffile,number);
		
		SDF *sdfp;
		SPHbody *body;
		
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
						"vsound", offsetof(SPHbody, vsound), &conf,
						"abar", offsetof(SPHbody, abar), &conf,
						"zbar", offsetof(SPHbody, zbar), &conf,
						"ax", offsetof(SPHbody, ax), &conf,
						"ay", offsetof(SPHbody, ay), &conf,
						"az", offsetof(SPHbody, az), &conf,
						"lax", offsetof(SPHbody, lax), &conf,
						"lay", offsetof(SPHbody, lay), &conf,
						"laz", offsetof(SPHbody, laz), &conf,
						"phi", offsetof(SPHbody, phi), &conf,
						"sigma", offsetof(SPHbody, sigma), &conf,
						"kappa", offsetof(SPHbody, kappa), &conf,
						"nbrs", offsetof(SPHbody, nbrs), &conf,
						"ident", offsetof(SPHbody, ident), &conf,
						"windid", offsetof(SPHbody, windid), &conf,
						//"useless", offsetof(SPHbody, useless), &conf,
						NULL);
		SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
		singlPrintf("%s has %d particles.\n", sdffile, gnobj);
		SDFclose(sdfp);
		
		for(j=0;j < nobj; j++)
		{
			if (body[j].ident == atoi(argv[2])) {	
				x		= body[j].x * dist_in_cm;
				y		= body[j].y * dist_in_cm;
				z		= body[j].z * dist_in_cm;
				rho		= body[j].rho * dens_in_gccm;
				pr		= body[j].pr * pressure_in_ergperccm;
				temp	= body[j].temp;
				u		= body[j].u * energy_in_erg/mass_in_g;
				vsound	= body[j].vsound * dist_in_cm;
				
				fprintf(stream,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",x,y,z,tpos,rho,pr,temp,u,vsound,body[j].abar,body[j].zbar);
				fprintf(stream,",%f",body[j].He4);
				fprintf(stream,",%f",body[j].C12);
				fprintf(stream,",%f",body[j].O16);
				fprintf(stream,",%f",body[j].Ne20);
				fprintf(stream,",%f",body[j].Mg24);
				fprintf(stream,",%f",body[j].Si28);
				fprintf(stream,",%f",body[j].S32);
				fprintf(stream,",%f",body[j].Ar36);
				fprintf(stream,",%f",body[j].Ca40);
				fprintf(stream,",%f",body[j].Ti44);
				fprintf(stream,",%f",body[j].Cr48);
				fprintf(stream,",%f",body[j].Fe52);
				fprintf(stream,",%f",body[j].Ni56);
				fprintf(stream,",%f\n",body[j].h*dist_in_cm);
			}
		}
		Free(body);
	}
		
	//close the stream file
	fclose(stream);
	
	return 0;
}