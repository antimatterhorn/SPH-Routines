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
double xmax,ymax,dist;
double time;
float tpos;
int choice,pixels,threads;
int i,j,k;
char sdffile[80];
char csvfile1[80];
char vszfile1[80];
char csvfile2[80];
char vszfile2[80];
char pdffile[80];

double dist_in_cm;
double mass_in_g;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;
double maxrho,maxtemp;

int usage()
{
	printf("\t Creates an interpolation plot of paticles on the x-y plane.\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_plot2d [sdf file] [enter '1'] [pixels] [xmax(code units)] [rhomax(cgs)] [tempmax(K)]\n");
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
	
	if (argc < 7){
		usage();
		return 0;
	}else {
		choice = atoi(argv[2]);
		pixels = atoi(argv[3]);
		xmax = atof(argv[4]);
		maxrho = atof(argv[5]);
        maxtemp = atof(argv[6]);
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
	
	double **img1, **img2, **(my_img1[threads]), **(my_img2[threads]);
	
	img1 = (double **) malloc(pixels*sizeof(double *));
	img2 = (double **) malloc(pixels*sizeof(double *));
	for(i=0;i<pixels;i++){
		img1[i]  = (double *) malloc(pixels*sizeof(double));
		img2[i]  = (double *) malloc(pixels*sizeof(double));
		for(j=0;j<pixels;j++)
		{
			img1[i][j] = img2[i][j] = 0;
		}
	}
	
	printf("computing %d threads...\t",threads);
	time = omp_get_wtime();
	
#pragma omp parallel private(k,i,j) \
	shared(pixels,xmax,ymax,body,choice,my_img1,my_img2,nobj,dens_in_gccm,threads) default(none)
	{
		int idx = omp_get_thread_num();
		my_img1[idx] = (double **) malloc(pixels*sizeof(double *));
		my_img2[idx] = (double **) malloc(pixels*sizeof(double *));

		for(i=0;i<pixels;i++){
			my_img1[idx][i] = (double *) malloc(pixels*sizeof(double));
			my_img2[idx][i] = (double *) malloc(pixels*sizeof(double));
			for(j=0;j<pixels;j++)
			{
				my_img1[idx][i][j] = my_img2[idx][i][j] = 0;
			}
		}
		
		double x,y,z,h,rho,temp,mass,kern,r,xc,yc;
		int hc,ic,jc;
		for(k=(double)idx/(double)threads*nobj;k<(double)(idx+1)/(double)threads*nobj;k++)
		{
			switch (choice) 
            {
                case 1:
                    if (fabs(body[k].z) < 2.0*body[k].h) /* x-y midplane */
                    {
                        x = body[k].x;
                        y = body[k].y;
                        z = body[k].z;
                        h = body[k].h;
                        
                        hc = 3*h*pixels/(2*xmax);
                        ic = x*(pixels/(2.0*xmax)) + pixels/2.0;
                        jc = -y*(pixels/(2.0*ymax)) + pixels/2.0;
                        
                        if((abs(ic)-hc) <= pixels || (abs(jc)-hc) <= pixels)
                        {
                            rho = body[k].rho;
                            temp = body[k].temp;
                            mass = body[k].mass;
                            
                            if(hc<1){ /* particle is smaller than a pixel, fill pixel instead */
                                kern = w(h,0);
                                my_img1[idx][ic][jc] += mass*kern*dens_in_gccm;
                                my_img2[idx][ic][jc] += mass*kern*dens_in_gccm*temp;					
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
                                            my_img1[idx][i][j] += mass*kern*dens_in_gccm;
                                            my_img2[idx][i][j] += mass*kern*dens_in_gccm*temp;
                                        }
                                    }
                                }
                            }  
                        }	
                    }
                    break;
                    
                case 2:
                    if (fabs(body[k].y) < 2.0*body[k].h) /* x-z midplane */
                    {
                        x = body[k].x;
                        y = body[k].z;
                        z = body[k].y;
                        h = body[k].h;
                        
                        hc = 3*h*pixels/(2*xmax);
                        ic = x*(pixels/(2.0*xmax)) + pixels/2.0;
                        jc = -y*(pixels/(2.0*ymax)) + pixels/2.0;
                        
                        if((abs(ic)-hc) <= pixels || (abs(jc)-hc) <= pixels)
                        {
                            rho = body[k].rho;
                            temp = body[k].temp;
                            mass = body[k].mass;
                            
                            if(hc<1){ /* particle is smaller than a pixel, fill pixel instead */
                                kern = w(h,0);
                                my_img1[idx][ic][jc] += mass*kern*dens_in_gccm;
                                my_img2[idx][ic][jc] += mass*kern*dens_in_gccm*temp;					
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
                                            my_img1[idx][i][j] += mass*kern*dens_in_gccm;
                                            my_img2[idx][i][j] += mass*kern*dens_in_gccm*temp;
                                        }
                                    }
                                }
                            }  
                        }	
                    }
                    break;
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
				
				
				img1[i][j] += my_img1[idx][i][j];
				img2[i][j] += my_img2[idx][i][j]; 
			}
		}
	}
	
	for(i=0;i<pixels;i++)
	{
		for(j=0;j<pixels;j++)
		{
				//printf("i=%d\tj=%d\n",i,j);
				img2[i][j] = ((img1[i][j]>0) ? img2[i][j]/img1[i][j] : 0); 

		}
	}
		
	//printf("freeing shared memory...\n");
	for(idx=0;idx<threads;idx++) 
	{
		free(my_img1[idx]);
		free(my_img2[idx]);
	}
	
	time = omp_get_wtime() - time;	
	printf("%3.2fs\nmaking files...\t\t",time);
	time = omp_get_wtime();

	snprintf(csvfile1, sizeof(csvfile1), "%s_rho.csv", argv[1]);
	snprintf(csvfile2, sizeof(csvfile2), "%s_temp.csv", argv[1]);

	FILE *stream1, *stream2, *fopen();
	
	stream1 = fopen(csvfile1,"w");
	stream2 = fopen(csvfile2,"w");
	for(i = 0; i < pixels; i++)
	{
		for(j = 0; j < pixels; j++)
		{
			fprintf(stream1,"%f ", (double)img1[i][j]);
			fprintf(stream2,"%f ", (double)img2[i][j]);
		}
		fprintf(stream1,"\n");
		fprintf(stream2,"\n");
	}
	
	fclose(stream1);
	fclose(stream2);

	
	FILE *vsz1, *vsz2;
	
	snprintf(vszfile1, sizeof(vszfile1), "%s_rho.vsz", argv[1]);						
	vsz1 = fopen(vszfile1,"w");
	
	fprintf(vsz1,"ImportFile2D(u'%s', [u'ds'], invertrows=False, invertcols=False,transpose=True, linked=True)\n",csvfile1);
	fprintf(vsz1,"Add('page', name='page1', autoadd=False)\n");
	fprintf(vsz1,"To('page1')\n");
	fprintf(vsz1,"Add('graph', name='graph1', autoadd=False)\n");
	fprintf(vsz1,"To('graph1')\n");
	fprintf(vsz1,"Set('leftMargin', u'2.3cm')\n");
	fprintf(vsz1,"Set('rightMargin', u'0.45cm')\n");
	fprintf(vsz1,"Set('topMargin', u'0.85cm')\n");
	fprintf(vsz1,"Set('bottomMargin', u'1.9cm')\n");
	fprintf(vsz1,"Add('axis', name='x', autoadd=False)\n");
    fprintf(vsz1,"To('x')\n");
    fprintf(vsz1,"Set('TickLabels/hide', True)\n");
    fprintf(vsz1,"Set('MajorTicks/hide', True)\n");
    fprintf(vsz1,"Set('MinorTicks/hide', True)\n");
    fprintf(vsz1,"To('..')\n");
	fprintf(vsz1,"Add('axis', name='y', autoadd=False)\n");
	fprintf(vsz1,"To('y')\n");
	fprintf(vsz1,"Set('direction', 'vertical')\n");
    fprintf(vsz1,"Set('TickLabels/hide', True)\n");
    fprintf(vsz1,"Set('MajorTicks/hide', True)\n");
    fprintf(vsz1,"Set('MinorTicks/hide', True)\n");
    fprintf(vsz1,"To('..')\n");
    fprintf(vsz1,"Add('axis', name='x-axis', autoadd=False)\n");
    fprintf(vsz1,"To('x-axis')\n");
    fprintf(vsz1,"Set('label', u'x [R_{\\odot}]')\n");
    fprintf(vsz1,"Set('min', %f)\n",-dist/6.96e10);
    fprintf(vsz1,"Set('max', %f)\n",dist/6.96e10);
    fprintf(vsz1,"Set('autoExtend', False)\n");
    fprintf(vsz1,"Set('Label/size', u'18pt')\n");
    fprintf(vsz1,"Set('TickLabels/size', u'18pt')\n");
    fprintf(vsz1,"To('..')\n");
    fprintf(vsz1,"Add('axis', name='y-axis', autoadd=False)\n");
    fprintf(vsz1,"To('y-axis')\n");
    fprintf(vsz1,"Set('label', u'y [R_{\\odot}]')\n");
    fprintf(vsz1,"Set('min', %f)\n",-dist/6.96e10);
    fprintf(vsz1,"Set('max', %f)\n",dist/6.96e10);
    fprintf(vsz1,"Set('autoExtend', False)\n");
    fprintf(vsz1,"Set('direction', u'vertical')\n");
    fprintf(vsz1,"Set('Label/size', u'18pt')\n");
    fprintf(vsz1,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz1,"To('..')\n");
	fprintf(vsz1,"Add('label', name='label1', autoadd=False)\n");
	fprintf(vsz1,"To('label1')\n");
	fprintf(vsz1,"Set('label', u't=%f')\n",tpos);
	fprintf(vsz1,"Set('xPos', [0.050000000000000003])\n");
	fprintf(vsz1,"Set('yPos', [0.05000000000000002])\n");
	fprintf(vsz1,"Set('Text/size', u'18pt')\n");
	fprintf(vsz1,"Set('Text/color', u'black')\n");
	fprintf(vsz1,"To('..')\n");
	fprintf(vsz1,"Add('colorbar', name='colorbar1', autoadd=False)\n");
	fprintf(vsz1,"To('colorbar1')\n");
	fprintf(vsz1,"Set('image', u'image1')\n");
	fprintf(vsz1,"Set('label', u'rho [g/cc]')\n");
	fprintf(vsz1,"Set('autoExtend', False)\n");
	fprintf(vsz1,"Set('autoExtendZero', False)\n");
	fprintf(vsz1,"Set('TickLabels/color', u'black')\n");
	fprintf(vsz1,"Set('vertPosn', u'top')\n");
	fprintf(vsz1,"To('..')\n");
	fprintf(vsz1,"Add('image', name='image1', autoadd=False)\n");
	fprintf(vsz1,"To('image1')\n");
	fprintf(vsz1,"Set('data', u'ds')\n");
	fprintf(vsz1,"Set('min', 0.01)\n");
	fprintf(vsz1,"Set('max', %f)\n",maxrho);
	fprintf(vsz1,"Set('colorScaling', u'log')\n");
	fprintf(vsz1,"Set('colorMap', u'heat')\n");
	fprintf(vsz1,"Set('colorInvert', True)\n");
	fprintf(vsz1,"Set('smooth', True)\n");
    fprintf(vsz1,"To('..')\n");
    fprintf(vsz1,"To('..')\n");
    fprintf(vsz1,"To('..')\n");

	
	fclose(vsz1);
	
	snprintf(vszfile2, sizeof(vszfile2), "%s_temp.vsz", argv[1]);						
	vsz2 = fopen(vszfile2,"w");
			
	fprintf(vsz2,"ImportFile2D(u'%s', [u'ds'], invertrows=False, invertcols=False, transpose=True, linked=True)\n",csvfile2);
	fprintf(vsz2,"Add('page', name='page1', autoadd=False)\n");
	fprintf(vsz2,"To('page1')\n");
	fprintf(vsz2,"Add('graph', name='graph1', autoadd=False)\n");
	fprintf(vsz2,"To('graph1')\n");
	fprintf(vsz2,"Set('leftMargin', u'2.3cm')\n");
	fprintf(vsz2,"Set('rightMargin', u'0.45cm')\n");
	fprintf(vsz2,"Set('topMargin', u'0.85cm')\n");
	fprintf(vsz2,"Set('bottomMargin', u'1.9cm')\n");
	fprintf(vsz2,"Add('axis', name='x', autoadd=False)\n");
    fprintf(vsz2,"To('x')\n");
    fprintf(vsz2,"Set('TickLabels/hide', True)\n");
    fprintf(vsz2,"Set('MajorTicks/hide', True)\n");
    fprintf(vsz2,"Set('MinorTicks/hide', True)\n");
    fprintf(vsz2,"To('..')\n");
	fprintf(vsz2,"Add('axis', name='y', autoadd=False)\n");
	fprintf(vsz2,"To('y')\n");
	fprintf(vsz2,"Set('direction', 'vertical')\n");
    fprintf(vsz2,"Set('TickLabels/hide', True)\n");
    fprintf(vsz2,"Set('MajorTicks/hide', True)\n");
    fprintf(vsz2,"Set('MinorTicks/hide', True)\n");
    fprintf(vsz2,"To('..')\n");
    fprintf(vsz2,"Add('axis', name='x-axis', autoadd=False)\n");
    fprintf(vsz2,"To('x-axis')\n");
    fprintf(vsz2,"Set('label', u'x [R_{\\odot}]')\n");
    fprintf(vsz2,"Set('min', %f)\n",-dist/6.96e10);
    fprintf(vsz2,"Set('max', %f)\n",dist/6.96e10);
    fprintf(vsz2,"Set('autoExtend', False)\n");
    fprintf(vsz2,"Set('Label/size', u'18pt')\n");
    fprintf(vsz2,"Set('TickLabels/size', u'18pt')\n");
    fprintf(vsz2,"To('..')\n");
    fprintf(vsz2,"Add('axis', name='y-axis', autoadd=False)\n");
    fprintf(vsz2,"To('y-axis')\n");
    fprintf(vsz2,"Set('label', u'y [R_{\\odot}]')\n");
    fprintf(vsz2,"Set('min', %f)\n",-dist/6.96e10);
    fprintf(vsz2,"Set('max', %f)\n",dist/6.96e10);
    fprintf(vsz2,"Set('autoExtend', False)\n");
    fprintf(vsz2,"Set('direction', u'vertical')\n");
    fprintf(vsz2,"Set('Label/size', u'18pt')\n");
    fprintf(vsz2,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz2,"To('..')\n");
	fprintf(vsz2,"Add('label', name='label1', autoadd=False)\n");
	fprintf(vsz2,"To('label1')\n");
	fprintf(vsz2,"Set('label', u't=%f')\n",tpos);
	fprintf(vsz2,"Set('xPos', [0.050000000000000003])\n");
	fprintf(vsz2,"Set('yPos', [0.90000000000000002])\n");
	fprintf(vsz2,"Set('Text/size', u'18pt')\n");
	fprintf(vsz2,"Set('Text/color', u'white')\n");
	fprintf(vsz2,"To('..')\n");
	fprintf(vsz2,"Add('colorbar', name='colorbar1', autoadd=False)\n");
	fprintf(vsz2,"To('colorbar1')\n");
	fprintf(vsz2,"Set('image', u'image1')\n");
	fprintf(vsz2,"Set('autoExtend', False)\n");
	fprintf(vsz2,"Set('autoExtendZero', False)\n");
	fprintf(vsz2,"Set('TickLabels/color', u'black')\n");
	fprintf(vsz2,"To('..')\n");
	fprintf(vsz2,"Add('image', name='image1', autoadd=False)\n");
	fprintf(vsz2,"To('image1')\n");
	fprintf(vsz2,"Set('data', u'ds')\n");
	fprintf(vsz2,"Set('min', 1000.0)\n");
	fprintf(vsz2,"Set('max', %f)\n",maxtemp);
	fprintf(vsz2,"Set('colorScaling', u'log')\n");
	fprintf(vsz2,"Set('colorMap', u'spectrum2')\n");
	fprintf(vsz2,"Set('colorInvert', True)\n");
	fprintf(vsz2,"Set('smooth', True)\n");
	
	
	fclose(vsz2);
		
	time = omp_get_wtime() - time;
	printf("%3.2fs\n",time);
		
	free(body);
	free(img1);
	free(img2);
	
	
	return 0;
}