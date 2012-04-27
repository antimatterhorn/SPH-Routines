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
#include <omp.h>

int gnobj, nobj;
int conf;
double xmax,ymax;
double time;
float tpos;
int choice,pixels,threads;
int i,j,k;
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
int openout, openup, opendown;
int open;

int usage()
{
	printf("\t Creates an interpolation plot of paticles on the x-y plane.\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_2grid [sdf file] [rho=1,temp=2] [pixels] [xmax(code units)] [rhomax(cgs)]\n");
	return 0;
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
					"durad", offsetof(SPHbody, durad), &conf,
					"openout", offsetof(SPHbody, openout), &conf,
					"openup", offsetof(SPHbody, openup), &conf,
					"opendown", offsetof(SPHbody, opendown), &conf,
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
	
	
	threads = omp_get_max_threads();
	
	double **img, **(my_img[threads]);
	
	img = (double **) malloc(pixels*sizeof(double *));
	for(i=0;i<pixels;i++){
		img[i]  = (double *) malloc(pixels*sizeof(double));
		for(j=0;j<pixels;j++)
		{
			img[i][j] = 0;
		}
	}
	
	printf("computing %d threads...\t",threads);
	time = omp_get_wtime();
	
#pragma omp parallel private(k,i,j,openout,openup,opendown,open) \
	shared(pixels,xmax,ymax,body,choice,my_img,nobj,dens_in_gccm,threads) default(none)
	{
		int idx = omp_get_thread_num();
		my_img[idx] = (double **) malloc(pixels*sizeof(double *));
		for(i=0;i<pixels;i++){
			my_img[idx][i] = (double *) malloc(pixels*sizeof(double));
			for(j=0;j<pixels;j++)
			{
				my_img[idx][i][j] = 0;
			}
		}
		
		double x,y,z,h,rho,temp,mass,kern,r,xc,yc;
		int hc,ic,jc;
		for(k=(double)idx/(double)threads*nobj;k<(double)(idx+1)/(double)threads*nobj;k++)
		{
			if (fabs(body[k].y) < 2.0*body[k].h) /* midplane */
			{
				x = body[k].x;
				y = body[k].z;
				z = body[k].y;
				h = body[k].h;
				hc = 3*h*pixels/(2*xmax);
				rho = body[k].rho;
				temp = body[k].temp;
				mass = body[k].mass;
				openout = body[k].openout;
				openup = body[k].openup;
				opendown = body[k].opendown;
				
				open = fmax(fmax(openup,opendown),openout);
				
				ic = x*(pixels/(2.0*xmax)) + pixels/2.0;
				jc = -y*(pixels/(2.0*ymax)) + pixels/2.0;
				
				if(hc<1){ /* particle is smaller than a pixel, fill pixel instead */
					kern = w(h,0);
					switch( choice )
					{
						case 1 :
							my_img[idx][ic][jc] += mass*kern*dens_in_gccm;
							break;
						case 2 :
							my_img[idx][ic][jc] += mass*temp/rho*kern;
							break;
					}
					
				}else{
					for(i=ic-hc;i<ic+hc+1;i++)
					{
						for(j=jc-hc;j<jc+hc+1;j++)
						{
							if(j>-1 && i>-1 && i<pixels && j<pixels){
								xc = (double)i/(double)pixels*(2.0*xmax) - xmax;
								yc = ymax - (double)j/(double)pixels*(2.0*ymax);
								
								r = sqrt(pow(xc-x,2.0)+pow(yc-y,2.0)+z*z);
								kern = w(h,r);
								switch( choice )
								{
									case 1 :
										my_img[idx][i][j] += mass*kern*dens_in_gccm;
										break;
									case 2 :
										my_img[idx][i][j] += mass*temp/rho*kern;
										break;
								}
							}
							
						}
					}
				}	
			}
		}	
		//printf("thread <%d> has finished.\n",idx);
	}
	
	time = omp_get_wtime() - time;
	printf("%3.2fs\ncollecting work...\t",time);
	time = omp_get_wtime();
	
	int idx;
	for(i=0;i<pixels;i++)
	{
		for(j=0;j<pixels;j++)
		{
			for(idx=0;idx<threads;idx++) img[i][j] += my_img[idx][i][j];
		}
	}
		
	for(idx=0;idx<threads;idx++) free(my_img[idx]);
	
	time = omp_get_wtime() - time;	
	printf("%3.2fs\nwriting the files...\t",time);
	time = omp_get_wtime();
	
	switch (choice){
		case 1:
			snprintf(csvfile, sizeof(csvfile), "%s_rho.csv", argv[1]);
			break;
		case 2:
			snprintf(csvfile, sizeof(csvfile), "%s_temp.csv", argv[1]);
			break;
	}
	
#pragma omp parallel num_threads(2) \
	default(none) shared(img,choice,csvfile,vszfile,pdffile,argv,pixels,tpos,maxrho) \
	private(i,j)
	{
	#pragma omp sections nowait
		{
		
		#pragma omp section
			{
				FILE *stream, *fopen();
				
				stream = fopen(csvfile,"w");
				
				for(i = 0; i < pixels; i++)
				{
					for(j = 0; j < pixels; j++)
					{
						fprintf(stream,"%f ", (double)img[i][j]);
					}
					fprintf(stream,"\n");
				}
				
				fclose(stream);
			}
			
		#pragma omp section
			{
				
				FILE *fopen(), *vsz;
				char *cwd = getcwd(NULL, 0);
				
				switch (choice) 
				{
					case 1:
						snprintf(vszfile, sizeof(vszfile), "%s_rho.vsz", argv[1]);
						snprintf(pdffile, sizeof(pdffile), "%s_rho.png", argv[1]);
						
						vsz = fopen(vszfile,"w");
						
						fprintf(vsz,"ImportFile2D(u'%s/%s', [u'ds'], invertrows=False, invertcols=False,transpose=True, linked=True)\n",cwd,csvfile);
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
						fprintf(vsz,"Set('min', 1e7)\n");
						fprintf(vsz,"Set('max', 2e8)\n");
						fprintf(vsz,"Set('colorScaling', u'linear')\n");
						fprintf(vsz,"Set('colorMap', u'spectrum2')\n");
						fprintf(vsz,"Set('colorInvert', False)\n");
						fprintf(vsz,"Set('smooth', True)\n");
						break;
				}
				
				fclose(vsz);
			}
		}
		
	}
		
	free(body);
	free(img);
	
	time = omp_get_wtime() - time;	
	printf("%3.2fs\n",time);
	
	return 0;
}