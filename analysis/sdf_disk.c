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
double xmax,zmax,x,y,z,rr,xmin,rc;
double rho,temp,mass,h;
float tpos;
int choice,pixels;
int i,j,imi,imj,ic,jc,r,k;
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
	printf("\t sdf_2grid [sdf file] [rho=1,temp=2] [pixels] [xint(cm)] [rhomax(cgs)]\n");
	return 0;
}

void x2i()
{
	ic = x*(pixels/(2.0*xmax)) + pixels/2.0;	
}

void z2j()
{
	jc = -z*(pixels/(2.0*zmax)) + pixels/2.0;
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
		xmin = atof(argv[4]);
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
	
	zmax = xmax = xmin*10;
	
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
	
	double **disk;
	
	disk = (double **) malloc(pixels/2.0*sizeof(double *));
	for(i=0;i<pixels/2.0;i++){
		disk[i] = (double *) malloc(2*sizeof(double));
	}
	
	SPHbody *p;
	i=0;
	for(p = body; p < body+gnobj; p++)
	{
		x = p->x * dist_in_cm;
		y = p->y * dist_in_cm;
		z = p->z * dist_in_cm;
		h = p->h * dist_in_cm;
		
		rc = sqrt(x*x+y*y+z*z);
		
		if (fabs(p->y) < p->h && rc > xmin) 
		{
			

			rho = p->rho * dens_in_gccm;
			temp = p->temp;
			mass = p->mass * mass_in_g;
			
			h = sqrt(pow(h,2.0)-pow(y,2.0));
			
			x2i();
			z2j();
			r = 2*h*(pixels/(2.0*xmax));
		
			if (r==0)
			{
				//printf("found that r was too small, filling only 1 pixel\n");
				if (ic >= 0 && ic < pixels && jc >= 0 && jc < pixels)
				{
					dens[ic][jc] += rho;
					switch( choice )
					{
						case 1 :
							img[ic][jc] +=rho*rho;
							break;
						case 2 :
							img[ic][jc] += temp*rho;
							break;
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
							if (rr <= r) // inside the circle
							{							
								dens[imi][imj] += rho;		
								switch( choice )
								{
									case 1 :
										img[imi][imj] += rho*rho*w(h,2*rr*(double)xmax/(double)pixels)*pi*pow(h,3.0);
										break;
									case 2 :
										img[imi][imj] += temp*rho;
										break;
								}
							}
						}
						
					}
				}	
			}
			

		}
		i++;
	}
	
	printf("image created.\n");
	
	//open the stream file
	
	switch (choice){
		case 1:
			snprintf(csvfile, sizeof(csvfile), "%s_drho.csv", argv[1]);
			break;
		case 2:
			snprintf(csvfile, sizeof(csvfile), "%s_dtemp.csv", argv[1]);
			break;
	}
	
	printf("saving image.\n");
	
	FILE *stream, *fopen(), *vsz;
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");
	
	for(i = 0; i < pixels; i++)
	{
		for(j = 0; j < pixels; j++)
		{
			if (dens[i][j] != 0)
			{
				fprintf(stream,"%f ", (double)img[i][j]/(double)dens[i][j]);
			}
			else
			{fprintf(stream,"%f ", 0.0);}
			
		}
		fprintf(stream,"\n");
	}
	
	//close the stream file
	fclose(stream);
	
	printf("creating plotfile.\n");
	
	char *cwd = getcwd(NULL, 0);

	
	switch (choice) {
		case 1:
			snprintf(vszfile, sizeof(vszfile), "%s_drho.vsz", argv[1]);
			snprintf(pdffile, sizeof(pdffile), "%s_drho.png", argv[1]);
			
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
			snprintf(vszfile, sizeof(vszfile), "%s_dtemp.vsz", argv[1]);
			snprintf(pdffile, sizeof(pdffile), "%s_dtemp.png", argv[1]);
			
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
	
	printf("calculating disk properties.\n");
	
	x = xmin*1.1;
	x2i();
	for(j=0;j<pixels/2.0;j++)
	{
		k = j+pixels/2.0;
		disk[j][0] = k*xmax/(pixels/2.0)/(6.955e10);
		disk[j][1] = (double)img[ic][k]/(double)dens[ic][k];
	}
	
	snprintf(csvfile, sizeof(csvfile), "%s_disk.csv", argv[1]);
	
	stream = fopen(csvfile,"w");
	
	fprintf(stream,"z,rho\n");
	
	for(i = 0; i < pixels/2.0; i++)
	{

		fprintf(stream,"%f,%3.2e\n",disk[i][0],disk[i][1]);
		
	}
	
	//close the stream file
	fclose(stream);
	
	return 0;
}