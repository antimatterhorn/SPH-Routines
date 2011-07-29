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
double xmax,ymax,x,y,z,rr,xc,yc,kern,r;
double rho,temp,mass,h;
float tpos;
int choice,pixels;
int i,j,imi,imj,ic,jc,hc;
char sdffile[80];
char csvfile[80];
char vszfile[80];
char pdffile[80];

double dist_in_cm;
double mass_in_g;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;
double maxrho;

int usage()
{
	printf("\t Creates an interpolation plot of paticles on the x-y plane.\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_2grid [sdf file] [rho=1,temp=2] [pixels] [xmax(code units)] [rhomax(cgs)]\n");
	return 0;
}

void x2i()
{
	ic = x*(pixels/(2.0*xmax)) + pixels/2.0;	
}

void y2j()
{
	jc = -y*(pixels/(2.0*ymax)) + pixels/2.0;
}

void i2x()
{
	xc = (double)i/(double)pixels*(2.0*xmax) - xmax;
}

void j2y()
{
	yc = ymax - (double)j/(double)pixels*(2.0*ymax);
}

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
	
	time_in_s				= 1;
	dist_in_cm				= 6.955e7;
	mass_in_g				= 1.989e27;
	dens_in_gccm			= mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	
	//maxrho					= 1e8;
	
	SDF *sdfp;
	SPHbody *body;
	
	if (argc < 6){
		usage();
		return 0;
	}else {
		choice = atoi(argv[2]);
		pixels = atoi(argv[3]);
		xmax = atof(argv[4]);
		maxrho = atof(argv[5]);
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
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
	
	ymax = xmax;
	
//	double **dens;
//	
//	dens = (double **) malloc(pixels*sizeof(double *));
//	for(i=0;i<pixels;i++){
//		dens[i]  = (double *) malloc(pixels*sizeof(double));
//	}
	
	double **img;
	
	img = (double **) malloc(pixels*sizeof(double *));
	for(i=0;i<pixels;i++){
		img[i]  = (double *) malloc(pixels*sizeof(double));
	}
	
	SPHbody *p;
	i=0;

	for(p = body; p < body+gnobj; p++)
	{
		if (fabs(p->z) < 2.0*p->h) /* midplane */
		{
			x = p->x;
			y = p->y;
			z = p->z;
			h = p->h;
			hc = 3*h*pixels/(2*xmax);
			rho = p->rho;
			temp = p->temp;
			mass = p->mass;
			
			x2i();
			y2j();
			
			if(hc<1){ /* particle is smaller than a pixel, fill pixel instead */
				kern = w(h,0);
				switch( choice )
				{
					case 1 :
						img[ic][jc] += mass*kern*dens_in_gccm;
						break;
					case 2 :
						img[ic][jc] += mass*temp/rho*kern;
						break;
				}
				
			}else{
				for(i=ic-hc;i<ic+hc+1;i++)
				{
					for(j=jc-hc;j<jc+hc+1;j++)
					{
						if(j>-1 && i>-1 && i<pixels && j<pixels){
							i2x();
							j2y();
							r = sqrt(pow(xc-x,2.0)+pow(yc-y,2.0)+z*z);
							kern = w(h,r);
							switch( choice )
							{
								case 1 :
									img[i][j] += mass*kern*dens_in_gccm;
									break;
								case 2 :
									img[i][j] += mass*temp/rho*kern;
									break;
							}
						}
						
					}
				}
			}
			
			
		}
	}
		
	//open the stream file
	
	switch (choice){
		case 1:
			snprintf(csvfile, sizeof(csvfile), "%s_rho.csv", argv[1]);
			break;
		case 2:
			snprintf(csvfile, sizeof(csvfile), "%s_temp.csv", argv[1]);
			break;
	}
	
	FILE *stream, *fopen(), *vsz;
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");
	
//	double maxv = 0.0;
//	double minv = pow(10.0,10.0);
	for(i = 0; i < pixels; i++)
	{
		for(j = 0; j < pixels; j++)
		{
				fprintf(stream,"%f ", (double)img[i][j]);
		}
		fprintf(stream,"\n");
	}
	
	//close the stream file
	fclose(stream);
	
	char *cwd = getcwd(NULL, 0);

	
	switch (choice) {
		case 1:
			snprintf(vszfile, sizeof(vszfile), "%s_rho.vsz", argv[1]);
			snprintf(pdffile, sizeof(pdffile), "%s_rho.png", argv[1]);
			
			vsz = fopen(vszfile,"w");
			
			fprintf(vsz,"ImportFile2D(u'%s/%s', [u'ds'], invertrows=False, invertcols=False, transpose=True, linked=True)\n",cwd,csvfile);
			fprintf(vsz,"Add('page', name='page1', autoadd=False)\n");
			fprintf(vsz,"To('page1')\n");
			fprintf(vsz,"Add('graph', name='graph1', autoadd=False)\n");
			fprintf(vsz,"To('graph1')\n");
			fprintf(vsz,"Set('leftMargin', u'0cm')\n");
			fprintf(vsz,"Set('rightMargin', u'0cm')\n");
			fprintf(vsz,"Set('topMargin', u'0cm')\n");
			fprintf(vsz,"Set('bottomMargin', u'0cm')\n");
			fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
			fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
			fprintf(vsz,"To('y')\n");
			fprintf(vsz,"Set('direction', 'vertical')\n");
			fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('label', name='label1', autoadd=False)\n");
			fprintf(vsz,"To('label1')\n");
			fprintf(vsz,"Set('label', u't=%f')\n",tpos);
			fprintf(vsz,"Set('xPos', [0.050000000000000003])\n");
			fprintf(vsz,"Set('yPos', [0.05000000000000002])\n");
			fprintf(vsz,"Set('Text/size', u'18pt')\n");
			fprintf(vsz,"Set('Text/color', u'black')\n");
			fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('colorbar', name='colorbar1', autoadd=False)\n");
			fprintf(vsz,"To('colorbar1')\n");
			fprintf(vsz,"Set('image', u'image1')\n");
			fprintf(vsz,"Set('label', u'rho [g/cc]')\n");
			fprintf(vsz,"Set('autoExtend', False)\n");
			fprintf(vsz,"Set('autoExtendZero', False)\n");
			fprintf(vsz,"Set('TickLabels/color', u'black')\n");
			fprintf(vsz,"Set('vertPosn', u'top')\n");
			fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('image', name='image1', autoadd=False)\n");
			fprintf(vsz,"To('image1')\n");
			fprintf(vsz,"Set('data', u'ds')\n");
			fprintf(vsz,"Set('min', 0.01)\n");
			fprintf(vsz,"Set('max', %f)\n",maxrho);
			fprintf(vsz,"Set('colorScaling', u'log')\n");
			fprintf(vsz,"Set('colorMap', u'heat2')\n");
			fprintf(vsz,"Set('colorInvert', True)\n");
			fprintf(vsz,"Set('smooth', True)\n");
			break;
		case 2:
			snprintf(vszfile, sizeof(vszfile), "%s_temp.vsz", argv[1]);
			snprintf(pdffile, sizeof(pdffile), "%s_temp.png", argv[1]);
			
			vsz = fopen(vszfile,"w");
			
			fprintf(vsz,"ImportFile2D(u'%s/%s', [u'ds'], invertrows=False, invertcols=False, transpose=True, linked=True)\n",cwd,csvfile);
			fprintf(vsz,"Add('page', name='page1', autoadd=False)\n");
			fprintf(vsz,"To('page1')\n");
			fprintf(vsz,"Add('graph', name='graph1', autoadd=False)\n");
			fprintf(vsz,"To('graph1')\n");
			fprintf(vsz,"Set('leftMargin', u'0cm')\n");
			fprintf(vsz,"Set('rightMargin', u'0cm')\n");
			fprintf(vsz,"Set('topMargin', u'0cm')\n");
			fprintf(vsz,"Set('bottomMargin', u'0cm')\n");
			fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
			fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
			fprintf(vsz,"To('y')\n");
			fprintf(vsz,"Set('direction', 'vertical')\n");
			fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('label', name='label1', autoadd=False)\n");
			fprintf(vsz,"To('label1')\n");
			fprintf(vsz,"Set('label', u't=%f')\n",tpos);
			fprintf(vsz,"Set('xPos', [0.050000000000000003])\n");
			fprintf(vsz,"Set('yPos', [0.90000000000000002])\n");
			fprintf(vsz,"Set('Text/size', u'18pt')\n");
			fprintf(vsz,"Set('Text/color', u'white')\n");
			fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('colorbar', name='colorbar1', autoadd=False)\n");
			fprintf(vsz,"To('colorbar1')\n");
			fprintf(vsz,"Set('image', u'image1')\n");
			fprintf(vsz,"Set('autoExtend', False)\n");
			fprintf(vsz,"Set('autoExtendZero', False)\n");
			fprintf(vsz,"Set('TickLabels/color', u'white')\n");
			fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('image', name='image1', autoadd=False)\n");
			fprintf(vsz,"To('image1')\n");
			fprintf(vsz,"Set('data', u'ds')\n");
			fprintf(vsz,"Set('min', 2e8)\n");
			fprintf(vsz,"Set('max', 'Auto')\n");
			fprintf(vsz,"Set('colorScaling', u'linear')\n");
			fprintf(vsz,"Set('colorMap', u'spectrum')\n");
			fprintf(vsz,"Set('colorInvert', False)\n");
			fprintf(vsz,"Set('smooth', True)\n");
			break;
	}

	fclose(vsz);
	
	return 0;
}