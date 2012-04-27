/*
 *  sdf2csv.c
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "../SPHbody.h"
#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
#include "Msgs.h"
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int gnobj, nobj;
int conf;
float x,y,rho,temp,mass,pressure,h,cs,vx,tpos,u,omega;
float bx,by,bbx,bby;
float pot;
int i,j,bi,bj;
int boolMade;
int nx;
float dx,xmax;

double dist_in_cm;
double mass_in_g;
double mass_in_msun;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;
double Gnewt;

char sdffile[80];
char csv1[80];
char csv2[80];
char csv3[80];
char vszfile[80];
char pdffile[80];
char cmd[80];

int usage()
{
	printf("\t Creates a map of grav-potential in x-y\n");
	printf("\t Usage: [required] {optional,default}\n");
	printf("\t sdf_potential [sdf file] [xmax(code units)] [omega(s-1)]\n");
	return 0;
}

int main(int argc, char **argv[])
{
	if (argc < 3)
	{
		usage();
		return 0;
	}
	xmax = atoi(argv[2]);
	omega = atof(argv[3]);
	
	char *cwd = getcwd(NULL, 0);
	nx = 100;
	dx = 2*xmax/(double)nx;
	
	time_in_s				= 1;
	dist_in_cm				= 6.955e7;
	mass_in_g				= 1.989e27;
	mass_in_msun			= mass_in_g / 1.99e33;
	dens_in_gccm			= mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	Gnewt					= 0.000393935;
	
	SDF *sdfp;
	SPHbody *body;
	
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
					"drho_dt", offsetof(SPHbody, drho_dt), &conf,
					"udot", offsetof(SPHbody, udot), &conf,
					"temp", offsetof(SPHbody, temp), &conf,
					"phi", offsetof(SPHbody,phi), &conf,
					NULL);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	
	printf("tpos = %3.5f\n",tpos);
	
	//first create potential map
	
	double **map;
	
	map = (double **) malloc(nx*sizeof(double *));
	for(i=0;i<nx;i++){
		map[i]  = (double *) malloc(nx*sizeof(double));
	}
	
	for(i=0;i<nx;i++){
		for(j=0;j<nx;j++){
			map[i][j] = 0;
		}
	}
	
	double **plot; //potx,pot
	
	plot = (double **) malloc(nx*sizeof(double *));
	for(i=0;i<nx;i++){
		plot[i]  = (double *) malloc(2*sizeof(double));
	}
	
	double **part; //x,y,phi,(h-z)
	
	part = (double **) malloc(nobj*sizeof(double *));
	for(i=0;i<nobj;i++){
		part[i]  = (double *) malloc(4*sizeof(double));
	}
	
	for(bi=0;bi<nobj;bi++){
		bx = body[bi].x;
		by = body[bi].y;					
				
		if(boolMade==0){
			part[bi][0] = bx*dist_in_cm;
			part[bi][1] = by*dist_in_cm;
			
			//for(bj=0;bj<nobj;bj++){
			//						if(bj!=bi){
			//							bbx = body[bj].x;
			//							bby = body[bj].y; 
			//							part[bi][2] -= Gnewt*body[bj].mass/sqrt((bx-bbx)*(bx-bbx)+(by-bby)*(by-bby));
			//						}
			//						//printf("%d-%d\n",bi,bj);
			//					}
			part[bi][2] = body[bi].phi;
			part[bi][2] -= 0.5*omega*omega*(bx*bx+by*by);
			part[bi][3] = 0.5*body[bi].h - body[bi].z;
			
			if(bi%1000 == 0) printf("%d\n",bi);
		}
	}
	
#pragma omp parallel for private(i,j,bi,boolMade,x,y,pot,bx,by)\
	default(none) shared(body,nx,xmax,dx,dist_in_cm,map,plot,Gnewt,omega,nobj)
	
		for(j=0;j<nx;j++){
			y = (-xmax + j*dx);
			
			for(i=0;i<nx;i++){
				x = (-xmax + i*dx);
				pot = 0;
				plot[i][0] = x*dist_in_cm;
				
				for(bi=0;bi<nobj;bi++){
					bx = body[bi].x;
					by = body[bi].y;					
					pot -= Gnewt*body[bi].mass/sqrt((x-bx)*(x-bx)+(y-by)*(y-by));
				}
				pot -= 0.5*omega*omega*(x*x+y*y);
				map[j][i] = pot;
				if(j==nx>>1) plot[i][1] = pot;
				
				
			}
			printf("<%d> : %d : %d\n",omp_get_thread_num(),(i),(nx));
		}
		
	
		
	snprintf(csv1, sizeof(csv1), "%s_potmap.csv", argv[1]);
	snprintf(csv2, sizeof(csv2), "%s_map.csv", argv[1]);
	snprintf(csv3, sizeof(csv3), "%s_pot.csv", argv[1]);
	FILE *stream, *fopen(), *vsz;
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csv1,"w");
	
	//	double maxv = 0.0;
	//	double minv = pow(10.0,10.0);
	for(i = 0; i < nx; i++)
	{
		for(j = 0; j < nx; j++)
		{
			fprintf(stream,"%f ", map[i][j]);
		}
		fprintf(stream,"\n");
	}
	
	//close the stream file
	fclose(stream);
	
	stream = fopen(csv2,"w");
	fprintf(stream,"x,y,phi\n");
	for(i=0;i<nobj;i++){
		if(i%8==0) fprintf(stream,"%f,%f,%f\n",part[i][0],part[i][1],part[i][2]);
	}
	fclose(stream);
	
	stream = fopen(csv3,"w");
	fprintf(stream,"potx,pot\n");
	for(i=0;i<nx;i++){
		fprintf(stream,"%f,%f\n",plot[i][0],plot[i][1]);
	}
	fclose(stream);
	
		
	return 0;
}