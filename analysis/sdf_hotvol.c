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
int conf,id,pixels;
double xmax,ymax,x,y,z,rr,volume;
double rho,temp,mass,h,mintemp,minrho,xperpix;
double averho,avetemp,ytot,mtot,etot;
float tpos;
int i,j,imi,imj,ic,jc,r;
int imin,jmin,imax,jmax,pcount;
char sdffile[80];
char csvfile[80];
char csvfile2[80];
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
	printf("\t Usage: [required] {optional,default}\n");
	printf("\t sdf_matchhead [sdf file] [id] {min temp(K),5e8} {min rho(cgs),5e6} {pixels,1000} {xmax(code units),15}\n");
	return 0;
}

int checknbrs(int ipos, int jpos, int **array, int size_el)
{
	pcount = 0;
	
	imin = ipos-1;
	jmin = jpos-1;
	jmax = jpos+1;
	imax = ipos+1;
	
	if (jmin < 0) jmin = 0;
	if (jmax >= size_el) jmax = size_el - 1;
	if (imin < 0) imin = 0;
	if (imax >= size_el) imax = size_el - 1;
	
	for (ic = imin; ic <= imax; ic++)
	{
		for (jc = jmin; jc <= jmax; jc++)
		{
			pcount += array[ic][jc];
		}
	}
	if (pcount <= 1) array[ipos][jpos] = 0;
	return pcount;
}

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

int main(int argc, char **argv[])
{
	
	time_in_s				= 1;
	dist_in_cm				= 6.955e7;
	mass_in_g				= 1.989e27;
	dens_in_gccm			= mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	
	maxrho					= 1e9;
	
	SDF *sdfp;
	SPHbody *body;
	
	if (argc < 3){
		usage();
		return 0;
	}else {
		id = atoi(argv[2]);
	}
	
	if (argc < 4){
		mintemp = 5e8;
	}else {
		mintemp = atof(argv[3]);
	}
	
	if (argc < 5){
		minrho = 5e6;
	}else {
		minrho = atof(argv[4]);
	}
	
	if (argc < 6){
		pixels = 1000;
	}else {
		pixels = atoi(argv[5]);
	}
	
	if (argc < 7){
		xmax = 15;
	}else {
		xmax = atof(argv[6]);
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
	
	ymax = xmax = xmax * dist_in_cm;
	
	xperpix = 2*xmax/(double)pixels;
	
	double **dens;
	
	dens = (double **) malloc(pixels*sizeof(double *));
	for(i=0;i<pixels;i++){
		dens[i]  = (double *) malloc(pixels*sizeof(double));
	}
	
	double **img;
	int **img2;
	
	img = (double **) malloc(pixels*sizeof(double *));
	for(i=0;i<pixels;i++){
		img[i]  = (double *) malloc(pixels*sizeof(double));
	}
	
	img2 = (int **) malloc(pixels*sizeof(int *));
	for(i=0;i<pixels;i++){
		img2[i]  = (int *) malloc(pixels*sizeof(int));
	}
	
	//fill arrays with 0s?
	for(i = 0; i < pixels; i++)
	{
		for(j = 0; j < pixels; j++)
		{
			img[i][j]=dens[i][j]=img2[i][j]=0;
		}
	}
	
	
	averho = avetemp = mtot = etot = 0;
	
	SPHbody *p;
	i=0;
	for(p = body; p < body+gnobj; p++)
	{		
		if (fabs(p->z) < p->h) 
		{
			
			x = p->x * dist_in_cm;
			y = p->y * dist_in_cm;
			z = p->z * dist_in_cm;
			h = p->h * dist_in_cm;
			rho = p->rho * dens_in_gccm;
			temp = p->temp;
			mass = p->mass * mass_in_g;
			
			h = sqrt(pow(h,2.0)-pow(z,2.0));
			
			x2i();
			y2j();
			r = h*(pixels/(2.0*xmax));
			
			if (temp > mintemp & rho > minrho)
			{
				averho += rho*fabs(y);
				avetemp += temp*fabs(y);
				ytot += fabs(y);
				mtot += mass;
				etot += p->u * specenergy_in_ergperg * mass;
			}
			
			if (r==0)
			{
				if (ic >= 0 && ic < pixels && jc >= 0 && jc < pixels)
				{
					dens[ic][jc] += rho;
					img[ic][jc] +=rho*rho;
					img2[ic][jc] = 1;
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
							if (rr <= 2*r) // inside the circle
							{							
								dens[imi][imj] += rho;		
								img[imi][imj] += rho*rho*w(h,2*rr*(double)xmax/(double)pixels)*pi*pow(h,3.0);
								if (temp > mintemp & rho > minrho) img2[imi][imj] = 1;
							}
						}
						
					}
				}	
			}
		}
		i++;
	}
	
	
	volume = 0;
	for(i = 0; i < pixels; i++)
	{
		for(j = 0; j < pixels; j++)
		{
			if (checknbrs(i,j,img2,pixels) > 1)
			{
				r = fabs(ymax-(double)i*xperpix);
				volume += img2[i][j]*pow(xperpix,2.0)*(3.14159*r);
			}
		}
	}
	
	avetemp = avetemp/ytot;
	averho = averho/ytot;
	printf("rhobar\t= %e g/cc\n",averho);
	printf("tempbar\t= %e K\n",avetemp);
	printf("volume\t= %e cc\n",volume);
	printf("mass\t= %e g\n",mtot);
	printf("energy\t= %e erg\n",etot);
	
	//open the stream file
	snprintf(csvfile, sizeof(csvfile), "%s.csv", argv[1]);
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
	
	snprintf(csvfile2, sizeof(csvfile2), "%s_crit.csv",argv[1]);
	
	stream = fopen(csvfile2,"w");
		
	for(i = 0; i < pixels; i++)
	{
		for(j = 0; j < pixels; j++)
		{
				fprintf(stream,"%d ", img2[i][j]);
		}
		fprintf(stream,"\n");
	}
	
	//close the stream file
	fclose(stream);
	
	char *cwd = getcwd(NULL, 0);


	snprintf(vszfile, sizeof(vszfile), "%s_rho.vsz", argv[1]);
	snprintf(pdffile, sizeof(pdffile), "%s_rho.png", argv[1]);
	
	vsz = fopen(vszfile,"w");
	
	fprintf(vsz,"AddImportPath('%s')\n",cwd);
	fprintf(vsz,"ImportFile2D(u'%s', [u'ds'], invertrows=False, invertcols=False, transpose=True, linked=True)\n",csvfile);
	fprintf(vsz,"ImportFile2D(u'%s', [u'dsc'], invertrows=False, invertcols=False, transpose=True, linked=True)\n",csvfile2);
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
	fprintf(vsz,"Set('TickLabels/color', u'black')\n");
	fprintf(vsz,"Set('vertPosn', u'top')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('image', name='image2', autoadd=False)\n");
	fprintf(vsz,"To('image2')\n");
	fprintf(vsz,"Set('data', u'dsc')\n");
	fprintf(vsz,"Set('colorMap', u'transblack')\n");
	fprintf(vsz,"Set('colorInvert', True)\n");
	fprintf(vsz,"Set('transparency', 80)\n");
	fprintf(vsz,"Set('smooth', True)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('image', name='image1', autoadd=False)\n");
	fprintf(vsz,"To('image1')\n");
	fprintf(vsz,"Set('data', u'ds')\n");
	fprintf(vsz,"Set('min', 1.0)\n");
	fprintf(vsz,"Set('max', %f)\n",maxrho);
	fprintf(vsz,"Set('colorScaling', u'log')\n");
	fprintf(vsz,"Set('colorMap', u'heat')\n");
	fprintf(vsz,"Set('colorInvert', True)\n");

	fclose(vsz);
	
	return 0;
}