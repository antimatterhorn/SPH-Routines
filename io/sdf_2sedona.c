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
double xmax,ymax,dist,x,z;
double time;
float tpos,norm;
int choice,pixels,threads;
int i,j,k;
char sdffile[80];
char asciifile[80];

double dist_in_cm;
double mass_in_g;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;
double maxrho,maxtemp;

int bj=128,bi=256;



typedef struct {
	double x, z;             /* position of body */							//3
	float vx, vy;     /* velocity of body */								//4
	float rho;            /* density of body */
	float temp;           /* temperature of body */								//11
	float abund[13];
} asciibody;

int usage()
{
	printf("\t Creates an interpolated ascii output of paticles on the x-z plane.\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_plot2d [sdf file] [xmax(code units)]\n");
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

double wstep(double hi, double ri)
{
	double v = fabs(ri)/hi;
	double coef = pow(pi,-1.0)*pow(hi,-3.0);
	if (v <= 2 && v >= 0)
	{
		return coef;
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
	
	if (argc < 3){
		usage();
		return 0;
	}else {
		pixels = 256;
		xmax = atof(argv[2]);
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
    dist = xmax*dist_in_cm;
	
//	double **dens;
//	
//	dens = (double **) malloc(pixels*sizeof(double *));
//	for(i=0;i<pixels;i++){
//		dens[i]  = (double *) malloc(pixels*sizeof(double));
//	}
	
	
	threads = omp_get_max_threads();
	
	double **imgrho, **(my_imgrho[threads]);
	double **imgtemp, **(my_imgtemp[threads]);
	double **imgvx, **(my_imgvx[threads]);
	double **imgvy, **(my_imgvy[threads]);
	double **imgHe4, **(my_imgHe4[threads]);
	double **imgC12, **(my_imgC12[threads]);
	double **imgO16, **(my_imgO16[threads]);
	double **imgNe20, **(my_imgNe20[threads]);
	double **imgMg24, **(my_imgMg24[threads]);
	double **imgSi28, **(my_imgSi28[threads]);
	double **imgS32, **(my_imgS32[threads]);
	double **imgAr36, **(my_imgAr36[threads]);
	double **imgCa40, **(my_imgCa40[threads]);
	double **imgTi44, **(my_imgTi44[threads]);
	double **imgCr48, **(my_imgCr48[threads]);
	double **imgFe52, **(my_imgFe52[threads]);
	double **imgNi56, **(my_imgNi56[threads]);
	
	imgrho = (double **) malloc(pixels*sizeof(double *));
	imgtemp = (double **) malloc(pixels*sizeof(double *));
	imgvx = (double **) malloc(pixels*sizeof(double *));
	imgvy = (double **) malloc(pixels*sizeof(double *));
	imgHe4 = (double **) malloc(pixels*sizeof(double *));
	imgC12 = (double **) malloc(pixels*sizeof(double *));
	imgO16 = (double **) malloc(pixels*sizeof(double *));
	imgNe20 = (double **) malloc(pixels*sizeof(double *));
	imgMg24 = (double **) malloc(pixels*sizeof(double *));
	imgSi28 = (double **) malloc(pixels*sizeof(double *));
	imgS32 = (double **) malloc(pixels*sizeof(double *));
	imgAr36 = (double **) malloc(pixels*sizeof(double *));
	imgCa40 = (double **) malloc(pixels*sizeof(double *));
	imgTi44 = (double **) malloc(pixels*sizeof(double *));
	imgCr48 = (double **) malloc(pixels*sizeof(double *));
	imgFe52 = (double **) malloc(pixels*sizeof(double *));
	imgNi56 = (double **) malloc(pixels*sizeof(double *));
	for(i=0;i<pixels;i++){
		imgrho[i] = (double *) malloc(pixels*sizeof(double));
		imgtemp[i] = (double *) malloc(pixels*sizeof(double));
		imgvx[i] = (double *) malloc(pixels*sizeof(double));
		imgvy[i] = (double *) malloc(pixels*sizeof(double));
		imgHe4[i] = (double *) malloc(pixels*sizeof(double));
		imgC12[i] = (double *) malloc(pixels*sizeof(double));
		imgO16[i] = (double *) malloc(pixels*sizeof(double));
		imgNe20[i] = (double *) malloc(pixels*sizeof(double));
		imgMg24[i] = (double *) malloc(pixels*sizeof(double));
		imgSi28[i] = (double *) malloc(pixels*sizeof(double));
		imgS32[i] = (double *) malloc(pixels*sizeof(double));
		imgAr36[i] = (double *) malloc(pixels*sizeof(double));
		imgCa40[i] = (double *) malloc(pixels*sizeof(double));
		imgTi44[i] = (double *) malloc(pixels*sizeof(double));
		imgCr48[i] = (double *) malloc(pixels*sizeof(double));
		imgFe52[i] = (double *) malloc(pixels*sizeof(double));
		imgNi56[i] = (double *) malloc(pixels*sizeof(double));
		for(j=0;j<pixels;j++)
		{
			imgrho[i][j] = 0.0;
			imgtemp[i][j] = 0.0;
			imgvx[i][j] = 0.0;
			imgvy[i][j] = 0.0;
			imgHe4[i][j] = 0.0;
			imgC12[i][j] = 0.0;
			imgO16[i][j] = 0.0;
			imgNe20[i][j] = 0.0;
			imgMg24[i][j] = 0.0;
			imgSi28[i][j] = 0.0;
			imgS32[i][j] = 0.0;
			imgAr36[i][j] = 0.0;
			imgCa40[i][j] = 0.0;
			imgTi44[i][j] = 0.0;
			imgCr48[i][j] = 0.0;
			imgFe52[i][j] = 0.0;
			imgNi56[i][j] = 0.0;
		}
	}
	
	printf("computing %d threads...\t",threads);
	time = omp_get_wtime();
	
#pragma omp parallel private(k,i,j) \
	shared(pixels,xmax,ymax,body,choice,nobj,dens_in_gccm,threads) default(none)\
    shared(my_imgrho) \
    shared(my_imgtemp) \
    shared(my_imgvx) \
    shared(my_imgvy) \
    shared(my_imgHe4) \
    shared(my_imgC12) \
    shared(my_imgO16) \
    shared(my_imgNe20) \
    shared(my_imgMg24) \
    shared(my_imgSi28) \
    shared(my_imgS32) \
    shared(my_imgAr36) \
    shared(my_imgCa40) \
    shared(my_imgTi44) \
    shared(my_imgCr48) \
    shared(my_imgFe52) \
    shared(my_imgNi56) 
	{
		int idx = omp_get_thread_num();
		my_imgrho[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgtemp[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgvx[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgvy[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgHe4[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgC12[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgO16[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgNe20[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgMg24[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgSi28[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgS32[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgAr36[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgCa40[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgTi44[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgCr48[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgFe52[idx] = (double **) malloc(pixels*sizeof(double *));
		my_imgNi56[idx] = (double **) malloc(pixels*sizeof(double *));

		for(i=0;i<pixels;i++){
			my_imgrho[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgtemp[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgvx[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgvy[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgHe4[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgC12[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgO16[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgNe20[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgMg24[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgSi28[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgS32[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgAr36[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgCa40[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgTi44[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgCr48[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgFe52[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_imgNi56[idx][i] = (double *) malloc(pixels*sizeof(double));
			for(j=0;j<pixels;j++)
			{
				my_imgrho[idx][i][j] = 0.0;
				my_imgtemp[idx][i][j] = 0.0;
				my_imgvx[idx][i][j] = 0.0;
				my_imgvy[idx][i][j] = 0.0;
				my_imgHe4[idx][i][j] = 0.0;
				my_imgC12[idx][i][j] = 0.0;
				my_imgO16[idx][i][j] = 0.0;
				my_imgNe20[idx][i][j] = 0.0;
				my_imgMg24[idx][i][j] = 0.0;
				my_imgSi28[idx][i][j] = 0.0;
				my_imgS32[idx][i][j] = 0.0;
				my_imgAr36[idx][i][j] = 0.0;
				my_imgCa40[idx][i][j] = 0.0;
				my_imgTi44[idx][i][j] = 0.0;
				my_imgCr48[idx][i][j] = 0.0;
				my_imgFe52[idx][i][j] = 0.0;
				my_imgNi56[idx][i][j] = 0.0;
			}
		}
		
		double x,y,z,h,rho,temp,mass,kern,r,xc,yc;
		int hc,ic,jc;
		for(k=(double)idx/(double)threads*nobj;k<(double)(idx+1)/(double)threads*nobj;k++)
		{

			if (fabs(body[k].y) < 2.0*body[k].h) /* x-z midplane */
			{
				x = body[k].x;
				y = body[k].z;
				z = body[k].y;
				h = body[k].h;
				
				hc = 3*h*pixels/(2*xmax);
				ic = x*(pixels/(2.0*xmax)) + pixels/2.0;
				jc = -y*(pixels/(2.0*ymax)) + pixels/2.0; //this is actually z now remember
				
				if((abs(ic)-hc) <= pixels || (abs(jc)-hc) <= pixels)
				{
					rho = body[k].rho;
					temp = body[k].temp;
					mass = body[k].mass;
					
					if(hc<1){ /* particle is smaller than a pixel, fill pixel instead */
						kern = pow(2.0*xmax/pixels,-3);
						my_imgrho[idx][ic][jc] += mass*kern*dens_in_gccm;
						my_imgtemp[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].temp;
						my_imgvx[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].vx;
						my_imgvy[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].vz;
						my_imgHe4[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].He4;
						my_imgC12[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].C12;
						my_imgO16[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].O16;
						my_imgNe20[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].Ne20;
						my_imgMg24[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].Mg24;
						my_imgSi28[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].Si28;
						my_imgS32[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].S32;
						my_imgAr36[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].Ar36;
						my_imgCa40[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].Ca40;
						my_imgTi44[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].Ti44;
						my_imgCr48[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].Cr48;
						my_imgFe52[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].Fe52;
						my_imgNi56[idx][ic][jc] += mass*kern*dens_in_gccm*body[k].Ni56;					
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
									my_imgrho[idx][i][j] += mass*kern*dens_in_gccm;
									my_imgtemp[idx][i][j] += mass*kern*dens_in_gccm*body[k].temp;
                                    my_imgvx[idx][i][j] += mass*kern*dens_in_gccm*body[k].vx;
                                    my_imgvy[idx][i][j] += mass*kern*dens_in_gccm*body[k].vz;
                                    my_imgHe4[idx][i][j] += mass*kern*dens_in_gccm*body[k].He4;
                                    my_imgC12[idx][i][j] += mass*kern*dens_in_gccm*body[k].C12;
                                    my_imgO16[idx][i][j] += mass*kern*dens_in_gccm*body[k].O16;
                                    my_imgNe20[idx][i][j] += mass*kern*dens_in_gccm*body[k].Ne20;
                                    my_imgMg24[idx][i][j] += mass*kern*dens_in_gccm*body[k].Mg24;
                                    my_imgSi28[idx][i][j] += mass*kern*dens_in_gccm*body[k].Si28;
                                    my_imgS32[idx][i][j] += mass*kern*dens_in_gccm*body[k].S32;
                                    my_imgAr36[idx][i][j] += mass*kern*dens_in_gccm*body[k].Ar36;
                                    my_imgCa40[idx][i][j] += mass*kern*dens_in_gccm*body[k].Ca40;
                                    my_imgTi44[idx][i][j] += mass*kern*dens_in_gccm*body[k].Ti44;
                                    my_imgCr48[idx][i][j] += mass*kern*dens_in_gccm*body[k].Cr48;
                                    my_imgFe52[idx][i][j] += mass*kern*dens_in_gccm*body[k].Fe52;
                                    my_imgNi56[idx][i][j] += mass*kern*dens_in_gccm*body[k].Ni56;
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
			for(idx=0;idx<threads;idx++) 
			{
				
                imgrho[i][j] += my_imgrho[idx][i][j];
                imgtemp[i][j] += my_imgtemp[idx][i][j];
                imgvx[i][j] += my_imgvx[idx][i][j];
                imgvy[i][j] += my_imgvy[idx][i][j];
                imgHe4[i][j] += my_imgHe4[idx][i][j];
                imgC12[i][j] += my_imgC12[idx][i][j];
                imgO16[i][j] += my_imgO16[idx][i][j];
                imgNe20[i][j] += my_imgNe20[idx][i][j];
                imgMg24[i][j] += my_imgMg24[idx][i][j];
                imgSi28[i][j] += my_imgSi28[idx][i][j];
                imgS32[i][j] += my_imgS32[idx][i][j];
                imgAr36[i][j] += my_imgAr36[idx][i][j];
                imgCa40[i][j] += my_imgCa40[idx][i][j];
                imgTi44[i][j] += my_imgTi44[idx][i][j];
                imgCr48[i][j] += my_imgCr48[idx][i][j];
                imgFe52[i][j] += my_imgFe52[idx][i][j];
                imgNi56[i][j] += my_imgNi56[idx][i][j];

			}
		}
	}
	
	for(i=0;i<pixels;i++)
	{
		for(j=0;j<pixels;j++)
		{
			imgtemp[i][j] = ((imgrho[i][j]>0.0) ? imgtemp[i][j]/imgrho[i][j] : 0.0);
            imgvx[i][j] = ((imgrho[i][j]>0.0) ? imgvx[i][j]/imgrho[i][j] : 0.0);
            imgvy[i][j] = ((imgrho[i][j]>0.0) ? imgvy[i][j]/imgrho[i][j] : 0.0);
            imgHe4[i][j] = ((imgrho[i][j]>0.0) ? imgHe4[i][j]/imgrho[i][j] : 0.0);
            imgC12[i][j] = ((imgrho[i][j]>0.0) ? imgC12[i][j]/imgrho[i][j] : 0.0);
            imgO16[i][j] = ((imgrho[i][j]>0.0) ? imgO16[i][j]/imgrho[i][j] : 0.0);
            imgNe20[i][j] = ((imgrho[i][j]>0.0) ? imgNe20[i][j]/imgrho[i][j] : 0.0);
            imgMg24[i][j] = ((imgrho[i][j]>0.0) ? imgMg24[i][j]/imgrho[i][j] : 0.0);
            imgSi28[i][j] = ((imgrho[i][j]>0.0) ? imgSi28[i][j]/imgrho[i][j] : 0.0);
            imgS32[i][j] = ((imgrho[i][j]>0.0) ? imgS32[i][j]/imgrho[i][j] : 0.0);
            imgAr36[i][j] = ((imgrho[i][j]>0.0) ? imgAr36[i][j]/imgrho[i][j] : 0.0);
            imgCa40[i][j] = ((imgrho[i][j]>0.0) ? imgCa40[i][j]/imgrho[i][j] : 0.0);
            imgTi44[i][j] = ((imgrho[i][j]>0.0) ? imgTi44[i][j]/imgrho[i][j] : 0.0);
            imgCr48[i][j] = ((imgrho[i][j]>0.0) ? imgCr48[i][j]/imgrho[i][j] : 0.0);
            imgFe52[i][j] = ((imgrho[i][j]>0.0) ? imgFe52[i][j]/imgrho[i][j] : 0.0);
            imgNi56[i][j] = ((imgrho[i][j]>0.0) ? imgNi56[i][j]/imgrho[i][j] : 0.0);
		}
	}
    
    for(i=0;i<pixels;i++)
	{
		for(j=0;j<pixels;j++)
		{
            norm = 0.0;
            norm += imgHe4[i][j];
            norm += imgC12[i][j];
            norm += imgO16[i][j];
            norm += imgNe20[i][j];
            norm += imgMg24[i][j];
            norm += imgSi28[i][j];
            norm += imgS32[i][j];
            norm += imgAr36[i][j];
            norm += imgCa40[i][j];
            norm += imgTi44[i][j];
            norm += imgCr48[i][j];
            norm += imgFe52[i][j];
            norm += imgNi56[i][j];
            if(norm>0)
            {
                imgHe4[i][j] = imgHe4[i][j]/norm;
                imgC12[i][j] = imgC12[i][j]/norm;
                imgO16[i][j] = imgO16[i][j]/norm;
                imgNe20[i][j] = imgNe20[i][j]/norm;
                imgMg24[i][j] = imgMg24[i][j]/norm;
                imgSi28[i][j] = imgSi28[i][j]/norm;
                imgS32[i][j] = imgS32[i][j]/norm;
                imgAr36[i][j] = imgAr36[i][j]/norm;
                imgCa40[i][j] = imgCa40[i][j]/norm;
                imgTi44[i][j] = imgTi44[i][j]/norm;
                imgCr48[i][j] = imgCr48[i][j]/norm;
                imgFe52[i][j] = imgFe52[i][j]/norm;
                imgNi56[i][j] = imgNi56[i][j]/norm;
            }
		}
	}
		
	//printf("freeing shared memory...\n");
	for(idx=0;idx<threads;idx++) 
	{
		free(my_imgrho[idx]);
        free(my_imgtemp[idx]);
        free(my_imgvx[idx]);
        free(my_imgvy[idx]);
        free(my_imgHe4[idx]);
        free(my_imgC12[idx]);
        free(my_imgO16[idx]);
        free(my_imgNe20[idx]);
        free(my_imgMg24[idx]);
        free(my_imgSi28[idx]);
        free(my_imgS32[idx]);
        free(my_imgAr36[idx]);
        free(my_imgCa40[idx]);
        free(my_imgTi44[idx]);
        free(my_imgCr48[idx]);
        free(my_imgFe52[idx]);
        free(my_imgNi56[idx]);

	}
	
	time = omp_get_wtime() - time;	
	printf("%3.2fs\nmaking files...\t\t",time);
	time = omp_get_wtime();

    snprintf(asciifile, sizeof(asciifile), "%s.ascii", argv[1]);
	FILE *stream, *fopen();
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(asciifile,"w");

     fprintf(stream,"** %s **\n",asciifile);
     fprintf(stream,"box size:\t %d x %d\n",bj,bi);
     fprintf(stream,"dimensions:\t %3.2ecm x %3.2ecm\n",xmax*dist_in_cm,xmax*2.0*dist_in_cm);
     fprintf(stream,"t-pos:\t %3.2f\n",tpos);
     fprintf(stream,"**************************\n");
     fprintf(stream,"\n");

	fprintf(stream,"x(cm) z(cm) rho(g/cc) vx(cm/s) vz(cm/s) T(K) He4 C12 O16 Ne20 Mg24 Si28 S32 Ar36 Ca40 Ti44 Cr48 Fe52 Ni56\n");

    for(i=0;i<pixels;i++)
	{
		for(j=0;j<pixels;j++)
		{
            //i spans x, j spans z
            x = (double)i/(double)pixels*2.0*xmax-xmax;
            z = ymax - (double)j/(double)pixels*2.0*ymax;
            
            if(x>=0.0)
            {
                fprintf(stream,"%03.5e %03.5e ",x*dist_in_cm,z*dist_in_cm);
                fprintf(stream,"%03.2e ",((imgrho[i][j]>1e-12) ? imgrho[i][j] : 1e-12));
                fprintf(stream,"%03.2e %03.2e ",imgvx[i][j]*dist_in_cm,imgvy[i][j]*dist_in_cm);
                fprintf(stream,"%03.2e ",imgtemp[i][j]);
                fprintf(stream,"%1.5f ",imgHe4[i][j]);
                fprintf(stream,"%1.5f ",imgC12[i][j]);
                fprintf(stream,"%1.5f ",imgO16[i][j]);
                fprintf(stream,"%1.5f ",imgNe20[i][j]);
                fprintf(stream,"%1.5f ",imgMg24[i][j]);
                fprintf(stream,"%1.5f ",imgSi28[i][j]);
                fprintf(stream,"%1.5f ",imgS32[i][j]);
                fprintf(stream,"%1.5f ",imgAr36[i][j]);
                fprintf(stream,"%1.5f ",imgCa40[i][j]);
                fprintf(stream,"%1.5f ",imgTi44[i][j]);
                fprintf(stream,"%1.5f ",imgCr48[i][j]);
                fprintf(stream,"%1.5f ",imgFe52[i][j]);
                fprintf(stream,"%1.5f\n",imgNi56[i][j]);
            }
        }
    }

	
	fclose(stream);
	
		
	time = omp_get_wtime() - time;
	printf("%3.2fs\n",time);
		
	free(body);
	free(imgrho);
    free(imgtemp);
    free(imgvx);
    free(imgvy);
    free(imgHe4);
    free(imgC12);
    free(imgO16);
    free(imgNe20);
    free(imgMg24);
    free(imgSi28);
    free(imgS32);
    free(imgAr36);
    free(imgCa40);
    free(imgTi44);
    free(imgCr48);
    free(imgFe52);
    free(imgNi56);

	
	
	return 0;
}