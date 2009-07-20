/*
 *  sdf2csv.c
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "sdf_2grid.h"
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
float xmax,ymax,x,y,z,rr;
float rho,temp,mass,h;
int choice,pixels;
int i,j,imi,imj,ic,jc,r;
char sdffile[80];
char csvfile[80];

void x2i()
{
	ic = x*(pixels/(2.0*xmax)) + pixels/2.0;	
}

void y2j()
{
	jc = -y*(pixels/(2.0*ymax)) + pixels/2.0;
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
					"grav_mass", offsetof(SPHbody, grav_mass), &conf,
					"phi", offsetof(SPHbody, phi), &conf,
					"tacc", offsetof(SPHbody, tacc), &conf,
					"idt", offsetof(SPHbody, idt), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	
	singlPrintf("%s has %d particles.\n", sdffile, gnobj);
	
	printf("\nOutput Parameters \n-----------------\nChoose a Variable:\n");
	printf("1) Density \n2) Temperature\n\nChoice:");
	scanf("%d", &choice);
	
	printf("Pixels on a side:");
	scanf("%d", &pixels);
	
	printf("xmax:");
	scanf("%f", &xmax);
	ymax = xmax;
	
	double **dens;
	
	dens = (double **) malloc(pixels*sizeof(double *));
	for(i=0;i<pixels;i++){
		dens[i]  = (double *) malloc(pixels*sizeof(double));
	}
	
	double **img;
	
	img = (double **) malloc(pixels*sizeof(double *));
	for(i=0;i<pixels;i++){
		img[i]  = (double *) malloc(pixels*sizeof(double));
	}
	
	//fill arrays with 0s?
	for(i = 0; i < pixels; i++)
	{
		for(j = 0; j < pixels; j++)
		{
			img[i][j]=dens[i][j]=0;
		}
	}
	
	SPHbody *p;
	
	for(p = body; p < body+gnobj; p++)
	{
		singlPrintf("%d / %d\n",p->ident,gnobj-1);
		//cout << p << "/" << npart-1 << endl;
		
		if (fabs(p->z) < p->h) 
		{
			
			x = p->x;
			y = p->y;
			z = p->z;
			h = p->h;
			rho = p->rho;
			temp = p->temp;
			mass = p->mass;
			
			h = sqrt(pow(h,2.0)-pow(z,2.0));
			
			x2i();
			y2j();
			r = h*(pixels/(2.0*xmax));
			
			if (r==0)
			{
				printf("found that r was too small, filling only 1 pixel\n");
				if (ic >= 0 && ic < pixels && jc >= 0 && jc < pixels)
				{
					dens[ic][jc] += rho;
					switch( choice )
					{
						case 1 :
							img[ic][jc] += p->rho*rho;
						case 2 :
							img[ic][jc] += p->temp*rho;
					}
				}
				
			}
			else
			{
				for(imi = ic-2*r; imi < ic+2*r+1; imi++)
				{
					for(imj = jc-2*r; imj < jc+2*r+1; imj++)
					{
						if (imi >= 0 && imi < pixels && imj >= 0 && imj < pixels) // inside the image
						{
							rr = sqrt(pow((imi-ic),2.0)+pow((imj-jc),2.0));
							//cout << rr << "," << r << endl;
							if (rr <= 2*r) // inside the circle
							{							
								dens[imi][imj] += rho;		
								switch( choice )
								{
									case 1 :
										img[imi][imj] += p->rho*rho*w(h,2*rr*(double)xmax/(double)pixels)*pi*pow(h,3.0);
										break;
									case 2 :
										img[imi][imj] += p->temp*rho*w(h,2*rr*(double)xmax/(double)pixels)*pi*pow(h,3.0);
										break;
								}
							}
						}
						
					}
				}	
			}
		}
	}
	
	//open the stream file
	snprintf(csvfile, sizeof(csvfile), "%s.csv", sdffile);
	FILE *stream, *fopen();
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");
	
//	double maxv = 0.0;
//	double minv = pow(10.0,10.0);
	for(i = 0; i < pixels; i++)
	{
		for(j = 0; j < pixels; j++)
		{
			if (dens[i][j] != 0)
			{
				fprintf(stream,"%f ", (double)img[i][j]/(double)dens[i][j]);
//				fout << img[i][j]/dens[i][j] << " ";
//				if (img[i][j]/dens[i][j] < minv) {minv = img[i][j]/dens[i][j];}
//				if (img[i][j]/dens[i][j] > maxv) {maxv = img[i][j]/dens[i][j];}
			}
			else
			{fprintf(stream,"%f ", 0.0);}
			
		}
		fprintf(stream,"\n");
	}
	
	//close the stream file
	fclose(stream);
	
	return 0;
}