/*
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "sdf_cube.h"
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
float xmax,ymax,zmax,x,y,z,rr;
float rho,temp,mass,h;
int choice,pixels;
int i,j,k,imi,imj,imk,ic,jc,kc,r,ind;
char sdffile[80];
char csvfile[80];

void x2i()
{
	ic = x*(pixels/(2.0*xmax)) + pixels/2.0;	
}

void y2j()
{
	jc = y*(pixels/(2.0*ymax)) + pixels/2.0;
}

void z2k()
{
	kc = z*(pixels/(2.0*zmax)) + pixels/2.0;
}

double w(double hi, double ri)
{
	double v = ri/hi;
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
//					"drho_dt", offsetof(SPHbody, drho_dt), &conf,
//					"udot", offsetof(SPHbody, udot), &conf,
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
//					"abar", offsetof(SPHbody, abar), &conf,
//					"zbar", offsetof(SPHbody, zbar), &conf,
//					"ax", offsetof(SPHbody, ax), &conf,
//					"ay", offsetof(SPHbody, ay), &conf,
//					"az", offsetof(SPHbody, az), &conf,
//					"lax", offsetof(SPHbody, lax), &conf,
//					"lay", offsetof(SPHbody, lay), &conf,
//					"laz", offsetof(SPHbody, laz), &conf,
//					"gax", offsetof(SPHbody, gax), &conf,
//					"gay", offsetof(SPHbody, gay), &conf,
//					"gaz", offsetof(SPHbody, gaz), &conf,
//					"grav_mass", offsetof(SPHbody, grav_mass), &conf,
//					"phi", offsetof(SPHbody, phi), &conf,
//					"tacc", offsetof(SPHbody, tacc), &conf,
//					"idt", offsetof(SPHbody, idt), &conf,
//					"nbrs", offsetof(SPHbody, nbrs), &conf,
//					"ident", offsetof(SPHbody, ident), &conf,
//					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	
	singlPrintf("%s has %d particles.\n", sdffile, gnobj);
	
	printf("Blocks on a side:");
	scanf("%d", &pixels);
	
	printf("xmax:");
	scanf("%f", &xmax);
	ymax = zmax = xmax;
	
	double **outArray;
	
	outArray = (double **) malloc(pixels*pixels*pixels*sizeof(double *));
	for(i=0;i<(pixels*pixels*pixels);i++){
		outArray[i]  = (double *) malloc(6*sizeof(double));
	}
	
	for (i=0; i<pixels; i++) {
		for (j=0; j<pixels; j++) {
			for (k=0; k<pixels; k++) {
				ind = i*pixels*pixels + j*pixels + k;
				
				outArray[ind][0]=i;
				outArray[ind][1]=j;
				outArray[ind][2]=k;
				outArray[ind][3]=0;
				outArray[ind][4]=0;
				outArray[ind][5]=0;
			}
		}
	}
	
	SPHbody *p;
	
	printf("Anchoring to grid...\n");
	
	for(p = body; p < body+gnobj; p++)
	{
		//singlPrintf("%d / %d\n",p->ident,gnobj-1);
		//cout << p << "/" << npart-1 << endl;
		
		x = p->x;
		y = p->y;
		z = p->z;
		h = p->h;
		rho = p->rho;
		temp = p->temp;
		mass = p->mass;
		
		x2i();
		y2j();
		z2k();

		r = h*(pixels/(2.0*xmax));
		
		if (r==0)
		{
			//printf("found that r was too small, filling only 1 block\n");
			if (ic >= 0 && ic < pixels && jc >= 0 && jc < pixels && kc >= 0 && kc < pixels)
			{
				ind = ic*pixels*pixels + jc*pixels + kc;
				outArray[ind][5] += rho;
				outArray[ind][4] += temp*rho;
				outArray[ind][3] += rho*rho;
			}
			
		}
		else
		{
			for(imi = ic-2*r; imi < ic+2*r+1; imi++)
			{
				for(imj = jc-2*r; imj < jc+2*r+1; imj++)
				{
					for (imk = kc-2*r; imk < kc+2*r+1; imk++) {
						if (imi >= 0 && imi < pixels && imj >= 0 && imj < pixels && imk >= 0 && imk < pixels) // inside the image
						{
							rr = sqrt(pow((imi-ic),2.0)+pow((imj-jc),2.0)+pow((imk-kc),2.0));
							if (rr <= 2*r) // inside the circle
							{							
								ind = imi*pixels*pixels + imj*pixels + imk;
								outArray[ind][5] += rho;
								outArray[ind][4] += temp*rho*w(h,2*rr*(double)xmax/(double)pixels)*pi*pow(h,3.0);
								outArray[ind][3] += rho*rho*w(h,2*rr*(double)xmax/(double)pixels)*pi*pow(h,3.0);
							}
						}
					}
				}
			}	
		}
	}
	
	printf("Writing to file...\n");
	
	//open the stream file
	snprintf(csvfile, sizeof(csvfile), "%s.cube", sdffile);
	FILE *stream, *fopen();
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");
	
	for(i=0;i<(pixels*pixels*pixels);i++){
		if (outArray[i][5] != 0) {
			fprintf(stream,"%f %f %f %g %g\n",outArray[i][0],outArray[i][1],outArray[i][2],
										(double)outArray[i][3]/(double)outArray[i][5],
										(double)outArray[i][4]/(double)outArray[i][5]);
		}
		else {
			fprintf(stream,"%f %f %f %d %d\n",outArray[i][0],outArray[i][1],outArray[i][2],0,0);
		}
	}

	//close the stream file
	fclose(stream);
	
	return 0;
}