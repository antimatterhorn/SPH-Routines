/*
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *  
 *	Usage: sdf_opacity filename.sdf blocks xmax(code_units) rhomax(cgs) theta phi alpha
 *			alpha<1 is the factor by which the farthest length is scaled
 *			e.g. at the back edge of the cube, a distance x at the center appears
 *			as alpha*x
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
float xmax,ymax,zmax,x,y,z,zzm,rr,rhomax,pixscale,scalepix,xc,yc,zc,r;
float rho,temp,mass,h,tpos,scale,theta,phi,alpha;
int choice,pixels;
int i,j,k,imi,imj,imk,ic,jc,kc,hc,ind;
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
	printf("\t Creates a 3d plot of particles interpolated on a uniform grid.\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_opacity [sdf file] [blocks<500] [xmax(code units)] [rhomax(cgs)] {theta} {phi} {alpha<1}\n");
	return 0;
}

void x2i()
{
	ic = x*pixscale + (double)(pixels>>1);	
}

void y2j()
{
	jc = -y*pixscale + (double)(pixels>>1);
}

void z2k()
{
	kc = -z*pixscale + (double)(pixels>>1);
}

void i2x()
{
	xc = (double)i/(double)pixels*(2.0*xmax) - xmax;
}

void j2y()
{
	yc = xmax-(double)j/(double)pixels*(2.0*xmax);
}

void k2z()
{
	zc = xmax-(double)k/(double)pixels*(2.0*xmax);
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
	
	maxrho					= 1e6;	
	
	SDF *sdfp;
	SPHbody *body;
	
	if (argc < 5){
		usage();
		return 0;
	}
	else {
		pixels = atoi(argv[2]);
		xmax = atof(argv[3]);
		rhomax = atof(argv[4]);
	}

	//optional args below
	
	if (argc < 6){
		theta = 0;
	}
	else {
		theta = atoi(argv[5]);
	}
	
	if (argc < 7){
		phi = 0;
	}
	else {
		phi = atoi(argv[6]);
	}
	
	if (argc < 8){
		alpha = 1;
	}
	else {
		alpha = atof(argv[7]);
	}
	
	
	pixscale = (double)pixels/(double)xmax/2.0;
	scalepix = 1.0/pixscale;
	
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
	
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
	
	
	ymax = zmax = xmax;
	
	theta = theta * 3.14159/180;
	phi = phi * 3.14159/180;
	
	double ****outArray;
	
	outArray = (double ****) malloc(pixels*pixels*pixels*sizeof(double ***));
	for(i=0;i<pixels;i++){
		outArray[i] = (double ***) malloc(pixels*sizeof(double **));
		for (j=0; j<pixels; j++) {
			outArray[i][j] = (double **) malloc(pixels*sizeof(double *));
			for (k=0; k<pixels; k++) {
				outArray[i][j][k] = (double *) malloc(3*sizeof(double));
			}
		}
	}	

	for (i=0; i<pixels; i++) {
		for (j=0; j<pixels; j++) {
			for (k=0; k<pixels; k++) {
				outArray[i][j][k][0]=0;
				outArray[i][j][k][1]=0;
				outArray[i][j][k][2]=0;
			}
		}
	}
	
	double **imgArray;
	
	imgArray = (double **) malloc(pixels*sizeof(double *));
	for(i=0;i<pixels;i++) imgArray[i] = (double *) malloc(pixels*sizeof(double));
	
	for (i=0; i<pixels; i++) {
		for (j=0; j<pixels; j++) {
			imgArray[i][j] = 0;
		}
	}
	
	SPHbody *p;
	
	printf("Anchoring to grid...\n");
	
	for(p = body; p < body+gnobj; p++)
	{
		//singlPrintf("%d / %d\n",p->ident,gnobj-1);
		//cout << p << "/" << npart-1 << endl;
		
		ind++;
		if (ind % 10000 == 0) printf("%d/%d\t %3.3f\n",ind,nobj,(double)ind/(double)nobj);
		
		x = p->x;
		y = p->y;
		z = p->z;
		h = p->h;
		
		x = x*cos(theta) + z*sin(theta);
		y = y;
		z = -x*sin(theta) + z*cos(theta);
		
		x = x;
		y = y*cos(phi) - z*sin(phi);
		z = y*sin(phi) + z*cos(phi);
		
		zzm = z/zmax;
		
		x = x*(1-(1-alpha)*(1+zzm));
		y = y*(1-(1-alpha)*(1+zzm));
		h = h*(1-(1-alpha)*(1+zzm));
		
		hc = 3*h*pixscale*(1-(1-alpha)*(1+zzm));
		
		rho = p->rho;
		temp = p->temp;
		mass = p->mass;
		
		x2i();
		y2j();
		z2k();

		if(hc<1)
		{
			outArray[ic][jc][kc][2] += rho*dens_in_gccm;
			outArray[ic][jc][kc][1] += temp*mass/rho*w(h,0);
			outArray[ic][jc][kc][0] += mass*w(h,0)*dens_in_gccm;
		}
		else
		{
			for(i=ic-hc;i<ic+hc+1;i++)
			{
				for(j=jc-hc;j<jc+hc+1;j++)
				{
					for(k=kc-hc;k<kc+hc+1;k++)
					{
						if(j>-1 && i>-1 && i<pixels && j<pixels && k>-1 && k<pixels)
						{
							i2x();
							j2y();
							k2z();
							r = sqrt(pow(xc-x,2.0)+pow(yc-y,2.0)+pow(zc-z,2.0));
							outArray[i][j][k][2] += rho*dens_in_gccm;
							outArray[i][j][k][1] += temp*mass/rho*w(h,r);
							outArray[i][j][k][0] += mass*w(h,r)*dens_in_gccm;
						}
					}
				}
			}
		}
	}

	
	for (i=0; i<pixels; i++) {
		for (j=0; j<pixels; j++) {
			for (k=0; k<pixels; k++) {
				scale = outArray[i][j][k][2]/rhomax;
				imgArray[i][j] += (double)outArray[i][j][k][0];
				if (scale > 1) break;
			}
		}
	}
	
	
	printf("Writing to file...\n");
	
	snprintf(csvfile, sizeof(csvfile), "%s_%s.csv", argv[1],argv[5]);
	FILE *stream, *fopen(), *vsz;
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");
	
	//	double maxv = 0.0;
	//	double minv = pow(10.0,10.0);
	for(i = 0; i < pixels; i++)
	{
		for(j = 0; j < pixels; j++)
		{
				fprintf(stream,"%f ", (double)imgArray[i][j]);			
		}
		fprintf(stream,"\n");
	}
	
	//close the stream file
	fclose(stream);
	
	char *cwd = getcwd(NULL, 0);
	
	
	snprintf(vszfile, sizeof(vszfile), "%s_%s.vsz", argv[1],argv[5]);
	snprintf(pdffile, sizeof(pdffile), "%s_%s.png", argv[1],argv[5]);
	
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
			fprintf(vsz,"Set('Text/color', u'white')\n");
			fprintf(vsz,"To('..')\n");
			//fprintf(vsz,"Add('colorbar', name='colorbar1', autoadd=False)\n");
			//fprintf(vsz,"To('colorbar1')\n");
			//fprintf(vsz,"Set('image', u'image1')\n");
			//fprintf(vsz,"Set('label', u'rho [g/cc]')\n");
			//fprintf(vsz,"Set('TickLabels/color', u'black')\n");
			//fprintf(vsz,"Set('vertPosn', u'top')\n");
			//fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('image', name='image1', autoadd=False)\n");
			fprintf(vsz,"To('image1')\n");
			fprintf(vsz,"Set('data', u'ds')\n");
			fprintf(vsz,"Set('min', 'Auto')\n");
			fprintf(vsz,"Set('max', 'Auto')\n");
			fprintf(vsz,"Set('colorScaling', u'log')\n");
			fprintf(vsz,"Set('colorMap', u'heat')\n");
			fprintf(vsz,"Set('colorInvert', False)\n");
			fprintf(vsz,"Set('smooth', True)\n");

	fclose(vsz);
	
	free(imgArray);
	free(outArray);
	
	return 0;
}