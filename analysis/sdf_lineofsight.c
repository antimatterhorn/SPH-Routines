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
double pe,rangle,tpos,xmin,xmax,dx,mass;
int i,j,k,l,angle,dangle,bins,na,np;
double mass_in_g = 1.989E+27;        
double dist_in_cm = 6.955e7;     
double time_in_s = 1;
double energy_in_erg,dens_in_gccm,pressure_in_ergperccm,specenergy_in_ergperg;

char sdffile[80];
char csvfile[80];

double w(double hi, double ri)
{
	double v = fabs(ri)/hi;
	double coef = pow(pi,-1.0)*pow(hi,-3.0);
	if (v <= 1 && v >= 0)
	{
		return coef*(1.0-1.5*v*v+0.75*v*v*v);
	}
	else if (v > 1 && v <= 2)
	{
		return coef*0.25*pow(2.0-v,3.0);
	}
	else
	{
		return 0;
	}
}

int main(int argc, char **argv[])
{
	energy_in_erg = mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
	dens_in_gccm = mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
	pressure_in_ergperccm = energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	specenergy_in_ergperg = energy_in_erg/mass_in_g;
	
	SDF *sdfp;
	SPHbody *body;
	
	if (argc < 2){
		printf("SDF file: ");
		gets (argv[1]);
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
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);

	int dv = 20;
	np = nobj/dv;
	double **inArray;
	double **outArray;
	double tempArray[na];
	
	
	inArray = (double **) malloc(np*sizeof(double *));
	for(i=0;i<np;i++){
		inArray[i]  = (double *) malloc(2*sizeof(double));
	}
	
	for(i=0;i<np;i++){
		for(j=0;j<2;j++){
			inArray[i][j]  = 0.0;
		}
		
	}
	
	j=l=0;
	for(i = 0; i < nobj; i++)
	{
		if (i % dv == 0 && j < np){
			mass = body[i].mass; 
			inArray[j][0] = sqrt(pow(body[i].x,2.0)+pow(body[i].y,2.0)+pow(body[i].z,2.0));
			inArray[j][1] = mass*dv; 
			j++;
		}
	}
	
	printf("sorting %d particle positions...\n",np);
	
	for (i=0; i< (np -1); i++)    // element to be compared
    {
		for(j = (i+1); j < np; j++)   // rest of the elements
		{
			if (inArray[i][0] > inArray[j][0])          // descending order
			{
				for (k=0;k<2;k++) 
				{
					tempArray[k] = inArray[i][k];
					inArray[i][k] = inArray[j][k];
					inArray[j][k] = tempArray[k];
				}
			}
		}
	}
	
	mass = 0;
	for (i=0;i<np;i++)
	{
		mass	+= inArray[i][1];
		inArray[i][1] = mass;
	}	
	
	na		= 18;		//r,mass,rho,temp,vproj,abundances 1 - 13
	
	outArray = (double **) malloc(np*sizeof(double *));
#pragma omp parallel for default(none) \
	shared(outArray,na,np,mass_in_g,dist_in_cm,dx,xmin,inArray) \
	private(i,j)
	for(i=0;i<np;i++){
		outArray[i]  = (double *) malloc(na*sizeof(double));
		outArray[i][0] = inArray[i][0]*dist_in_cm;
		outArray[i][1] = inArray[i][1]*mass_in_g;
	}

	
	
	printf("making line-outs...\n");	
	
	dangle = 45;
	
	FILE *stream, *fopen();
	
	
	for(angle = 0; angle < 90+dangle; angle+=dangle)
	{
		rangle = angle*pi/180.0;
		snprintf(csvfile, sizeof(csvfile), "%s_%d.dat", argv[1],angle);
		

		
		#pragma omp parallel for default(none) \
			shared(outArray,inArray,dx,na,np,nobj,bins,xmin,dist_in_cm,dens_in_gccm,body,mass_in_g,rangle) \
			private(i,j,k,l)
		for(i=0;i<np;i++){
			for(j=2;j<na;j++){
				outArray[i][j]  = 0.0;
			}
			
			double x1,y1,z1,x2,y2,z2,r,kern,pd,h,vxp,vyp,vzp;
			
			x1 = inArray[i][0]*cos(rangle);
			y1 = inArray[i][0]*sin(rangle);
			z1 = 0.0;
			r  = inArray[i][0];
			l  = 0;
			
			for(k=0;k<nobj;k++)
			{
				x2	= body[k].x;
				y2	= body[k].y;
				z2	= body[k].z;
				pd	= sqrt(pow(x2-x1,2.0)+pow(y2-y1,2.0)+pow(z2-z1,2.0));
				h	= body[k].h;
				
				if(pd<h*3.0)
				{
					kern = w(h,pd);
					kern = kern*body[k].mass/body[k].rho;
					
					outArray[i][2] += kern*body[k].rho*dens_in_gccm;
					outArray[i][3] += kern*body[k].temp;
					
					/* projection of v onto rhat */
					vxp = body[k].vx*x1/r;
					vyp = body[k].vy*y1/r;
					vzp = body[k].vz*z1/r;
					
					outArray[i][4]	+= kern*(vxp+vyp+vzp)*dist_in_cm;	

					
					outArray[i][5]	+= kern*(body[k].He4);
					outArray[i][6]	+= kern*(body[k].C12);
					outArray[i][7]	+= kern*(body[k].O16);
					outArray[i][8]	+= kern*(body[k].Ne20);
					outArray[i][9]	+= kern*(body[k].Mg24);
					outArray[i][10] += kern*(body[k].Si28);
					outArray[i][11] += kern*(body[k].S32);
					outArray[i][12] += kern*(body[k].Ar36);
					outArray[i][13] += kern*(body[k].Ca40);
					outArray[i][14] += kern*(body[k].Ti44);
					outArray[i][15] += kern*(body[k].Cr48);
					outArray[i][16] += kern*(body[k].Fe52);
					outArray[i][17] += kern*(body[k].Ni56);
				}
			}

			

		}
		
		stream = fopen(csvfile,"w");
		
		fprintf(stream,"dist[cm]\tmass[gm]\trho[gm/cc]\ttemperature[k]\tvel_p[cm/s]\t");
		fprintf(stream,"X-He4\tX-C12\tX-O16\tX-Ne20\tX-Mg24\tX-Si28\tX-S32\tX-Ar36\tX-Ca40\tX-Ti44\tX-Cr48\tX-Fe52\tX-Ni56\n");
		
		for(i=0;i<np;i++)
		{
			for(j=0;j<na;j++)
			{
				fprintf(stream,"%3.3e",outArray[i][j]);
				if(j!=na-1) fprintf(stream,"\t");
			}
			fprintf(stream,"\n");
		}

		fclose(stream);
	
		printf("finished phi = %d\n",angle);
		
	}
	
	free(inArray);
	free(outArray);
	free(body);
	
	return 0;
}