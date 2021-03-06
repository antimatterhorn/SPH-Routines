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
int i,j,k,start,end,inc;
double iso[13],mass;

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
	
	strcpy(period,".");
	printf("start end increment: ");
	scanf("%d %d %d",&start,&end,&inc);
	
	FILE *stream, *fopen();
	/* declare a stream and prototype fopen */ 
	
	stream = fopen("out.csv","w");	
	
	fprintf(stream,"t,He4,C12,O16,Ne20,Mg24,Si28,S32,Ar36,Ca40,Ti44,Cr48,Fe52,Ni56\n");
	
	for(i=start;i<(end+1);i+=inc)
	{
		for(j=0;j<13;j++) iso[j]=0;
		mass=0;
		
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
						NULL);
		SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
		singlPrintf("%s has %d particles.\n", sdffile, gnobj);
		SDFclose(sdfp);
		
		for(k=0;k < nobj; k++)
		{
			iso[0] += body[k].He4*body[k].mass;
			iso[1] += body[k].C12*body[k].mass;
			iso[2] += body[k].O16*body[k].mass;
			iso[3] += body[k].Ne20*body[k].mass;
			iso[4] += body[k].Mg24*body[k].mass;
			iso[5] += body[k].Si28*body[k].mass;
			iso[6] += body[k].S32*body[k].mass;
			iso[7] += body[k].Ar36*body[k].mass;
			iso[8] += body[k].Ca40*body[k].mass;
			iso[9] += body[k].Ti44*body[k].mass;
			iso[10] += body[k].Cr48*body[k].mass;
			iso[11] += body[k].Fe52*body[k].mass;
			iso[12] += body[k].Ni56*body[k].mass;
			
			mass+=body[k].mass;
		}
		
		fprintf(stream,"%f",tpos);
		
		for(j=0;j<13;j++) fprintf(stream,",%f",iso[j]*mass_in_g/2e33);
		fprintf(stream,"\n");
		
		Free(body);
	}
		
	//close the stream file
	fclose(stream);
	
	return 0;
}