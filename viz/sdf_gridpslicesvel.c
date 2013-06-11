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
int conf,idx,ctr;
double xmax,ymax,zz,dz;
double time;
float tpos;
int choice,pixels,threads;
int i,j,k,zi;
char sdffile[80];
char csvfile1[80];
char vszfile1[80];
char csvfile2[80];
char csvfile3[80];
char csvfile4[80];
char csvfile5[80];

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
	
	if (argc < 6){
		usage();
		return 0;
	}else {
		choice = atoi(argv[2]);
		pixels = atoi(argv[3]);
		xmax = atof(argv[4]);
		maxrho = atof(argv[5]);
	}
    
    zz = -xmax;
    dz = 2*xmax / (double) pixels;
		
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
	
	threads = omp_get_max_threads();
	
	double **img1, **imgvx, **imgvy, **imgvz, **imgu;
    double **(my_img1[threads]), **(my_imgvx[threads]), **(my_imgvy[threads]), **(my_imgvz[threads]), **(my_imgu[threads]);
	
	img1 = (double **) malloc(pixels*sizeof(double *));
	imgvx = (double **) malloc(pixels*sizeof(double *));
    imgvy = (double **) malloc(pixels*sizeof(double *));
    imgvz = (double **) malloc(pixels*sizeof(double *));
    imgu = (double **) malloc(pixels*sizeof(double *));
	for(i=0;i<pixels;i++){
		img1[i]  = (double *) malloc(pixels*sizeof(double));
		imgvx[i]  = (double *) malloc(pixels*sizeof(double));
        imgvy[i]  = (double *) malloc(pixels*sizeof(double));
        imgvz[i]  = (double *) malloc(pixels*sizeof(double));
        imgu[i]  = (double *) malloc(pixels*sizeof(double));
		for(j=0;j<pixels;j++)
		{
			img1[i][j] = imgvx[i][j] = imgvy[i][j] = imgvz[i][j] = imgu[i][j] = 0;
		}
	}
    
    for(idx=0;idx<threads;idx++)
    {
        my_img1[idx] = (double **) malloc(pixels*sizeof(double *));
        my_imgvx[idx] = (double **) malloc(pixels*sizeof(double *));
        my_imgvy[idx] = (double **) malloc(pixels*sizeof(double *));
        my_imgvz[idx] = (double **) malloc(pixels*sizeof(double *));
        my_imgu[idx] = (double **) malloc(pixels*sizeof(double *));
        
        for(i=0;i<pixels;i++){
            my_img1[idx][i] = (double *) malloc(pixels*sizeof(double));
            my_imgvx[idx][i] = (double *) malloc(pixels*sizeof(double));
            my_imgvy[idx][i] = (double *) malloc(pixels*sizeof(double));
            my_imgvz[idx][i] = (double *) malloc(pixels*sizeof(double));
            my_imgu[idx][i] = (double *) malloc(pixels*sizeof(double));
            for(j=0;j<pixels;j++)
            {
                my_imgu[idx][i][j] = my_img1[idx][i][j] = my_imgvx[idx][i][j] = my_imgvy[idx][i][j] = my_imgvz[idx][i][j] = 0;
            }
        }
    }

	
	printf("computing %d threads...\t",threads);
	time = omp_get_wtime();

    snprintf(csvfile1, sizeof(csvfile1), "%3.0f_rho.csv", xmax);
    snprintf(csvfile2, sizeof(csvfile2), "%3.0f_u.csv", xmax);
    snprintf(csvfile3, sizeof(csvfile3), "%3.0f_vx.csv", xmax);
    snprintf(csvfile4, sizeof(csvfile4), "%3.0f_vy.csv", xmax);
    snprintf(csvfile5, sizeof(csvfile5), "%3.0f_vz.csv", xmax);
    
    FILE *stream1, *streamu, *streamvx, *streamvy, *streamvz, *fopen();
    
    stream1 = fopen(csvfile1,"w");
    streamu = fopen(csvfile2,"w");
    streamvx = fopen(csvfile3,"w");
    streamvy = fopen(csvfile4,"w");
    streamvz = fopen(csvfile5,"w");
    
    ctr = 0;
    
    for (zi=0;zi<pixels;zi++)
    {
        zz = zi*2.0*xmax/(double)pixels-xmax;
        for(i=0;i<pixels;i++){
            for(j=0;j<pixels;j++)
            {
                img1[i][j] = imgvx[i][j] = imgvy[i][j] = imgvz[i][j] = imgu[i][j] = 0;
            }
        }
        
#pragma omp parallel private(k,i,j,idx) \
    shared(pixels,xmax,ymax,body,choice,my_img1,my_imgu,my_imgvx,my_imgvy,my_imgvz,nobj,dens_in_gccm,threads,zz,dist_in_cm,specenergy_in_ergperg) default(none)
        {
            idx = omp_get_thread_num();            
            for(i=0;i<pixels;i++)
                for(j=0;j<pixels;j++)
                    my_imgu[idx][i][j] = my_img1[idx][i][j] = my_imgvx[idx][i][j] = my_imgvy[idx][i][j] = my_imgvz[idx][i][j] = 0;
            
            double x,y,z,h,rho,temp,mass,kern,r,xc,yc,u,vx,vy,vz;
            int hc,ic,jc;
            for(k=(double)idx/(double)threads*nobj;k<(double)(idx+1)/(double)threads*nobj;k++)
            {
                if (fabs(zz-body[k].z) < 2.0*body[k].h) 
                {
                    x = body[k].x;
                    y = body[k].y;
                    z = body[k].z;
                    vx = body[k].vx;
                    vy = body[k].vy;
                    vz = body[k].vz;
                    h = body[k].h;
                    hc = 3*h*pixels/(2*xmax);
                    rho = body[k].rho;
                    temp = body[k].temp;
                    mass = body[k].mass;
                    u = body[k].u;
                    
                    ic = x*(pixels/(2.0*xmax)) + pixels/2.0;
                    jc = -y*(pixels/(2.0*ymax)) + pixels/2.0;
                    
                    if(hc<1){ /* particle is smaller than a pixel, fill pixel instead */
                        kern = w(h,0);
                        my_img1[idx][ic][jc] += mass*kern*dens_in_gccm;
                        my_imgvx[idx][ic][jc] += mass*kern*dens_in_gccm*vx*dist_in_cm;
                        my_imgvy[idx][ic][jc] += mass*kern*dens_in_gccm*vy*dist_in_cm;
                        my_imgvz[idx][ic][jc] += mass*kern*dens_in_gccm*vz*dist_in_cm;
                        my_imgu[idx][ic][jc] += mass*kern*dens_in_gccm*u*specenergy_in_ergperg;					
                    }else{
                        for(i=ic-hc;i<ic+hc+1;i++)
                        {
                            for(j=jc-hc;j<jc+hc+1;j++)
                            {
                                if(j>-1 && i>-1 && i<pixels && j<pixels){
                                    xc = (double)i/(double)pixels*(2.0*xmax) - xmax;
                                    yc = ymax - (double)j/(double)pixels*(2.0*ymax);
                                    
                                    r = sqrt(pow(xc-x,2.0)+pow(yc-y,2.0)+pow(fabs(zz-z),2.0));
                                    kern = w(h,r);
                                    my_img1[idx][i][j] += mass*kern*dens_in_gccm;
                                    my_imgvx[idx][i][j] += mass*kern*dens_in_gccm*vx*dist_in_cm;
                                    my_imgvy[idx][i][j] += mass*kern*dens_in_gccm*vy*dist_in_cm;
                                    my_imgvz[idx][i][j] += mass*kern*dens_in_gccm*vz*dist_in_cm;
                                    my_imgu[idx][i][j] += mass*kern*dens_in_gccm*u*specenergy_in_ergperg;
                                }
                            }
                        }
                    }	
                }
            }	
        }
        
        time = omp_get_wtime() - time;
        //printf("%3.2fs\ncollecting work...\t",time);
        time = omp_get_wtime();
        
        for(i=0;i<pixels;i++)
        {
            for(j=0;j<pixels;j++)
            {
                for(idx=0;idx<threads;idx++) 
                { 
                    img1[i][j] += my_img1[idx][i][j];
                    imgu[i][j] += my_imgu[idx][i][j]; 
                    imgvx[i][j] += my_imgvx[idx][i][j]; 
                    imgvy[i][j] += my_imgvy[idx][i][j]; 
                    imgvz[i][j] += my_imgvz[idx][i][j]; 
                }
            }
        }
        
        for(i=0;i<pixels;i++)
        {
            for(j=0;j<pixels;j++)
            {
				imgu[i][j] = ((img1[i][j]>0) ? imgu[i][j]/img1[i][j] : 0); 
                imgvx[i][j] = ((img1[i][j]>0) ? imgvx[i][j]/img1[i][j] : 0); 
                imgvy[i][j] = ((img1[i][j]>0) ? imgvy[i][j]/img1[i][j] : 0); 
                imgvz[i][j] = ((img1[i][j]>0) ? imgvz[i][j]/img1[i][j] : 0); 
            }
        }
		
        //printf("freeing shared memory...\n");

        
        time = omp_get_wtime() - time;	
        //printf("%3.2fs\nmaking files...\t\t",time);
        time = omp_get_wtime();
        

        for(i = 0; i < pixels; i++)
        {
            for(j = 0; j < pixels; j++)
            {
                float imv1 = img1[i][j];
                float imvu = imgu[i][j];
                float imvvx = imgvx[i][j];
                float imvvy = imgvy[i][j];
                float imvvz = imgvz[i][j];
                
                fwrite(&imv1,sizeof(float),1,stream1); 
                fwrite(&imvu,sizeof(float),1,streamu); 
                fwrite(&imvvx,sizeof(float),1,streamvx); 
                fwrite(&imvvy,sizeof(float),1,streamvy); 
                fwrite(&imvvz,sizeof(float),1,streamvz); 
                ctr++;
//                fprintf(stream1,"%f ", img1[i][j]);
//                fprintf(stream2,"%f ", img2[i][j]);
            }
//            fprintf(stream1,"\n");
//            fprintf(stream2,"\n");
        }
        
//        fprintf(stream1,"\n\n");
//        fprintf(stream2,"\n\n");
        
        printf("%d/%d\n",zi,pixels);
    } 
    
    printf("counter = %d\n",ctr);
    
    fclose(stream1); 
    fclose(streamu);
    fclose(streamvx);
    fclose(streamvy);
    fclose(streamvz);
    

    for(idx=0;idx<threads;idx++) 
    {
        free(my_img1[idx]);
        free(my_imgu[idx]);
        free(my_imgvx[idx]);
        free(my_imgvy[idx]);
        free(my_imgvz[idx]);
    }

	
			
	time = omp_get_wtime() - time;
	printf("%3.2fs\n",time);
		
	free(body);
	free(img1);
	free(imgu);
    free(imgvx);
    free(imgvy);
    free(imgvz);
    
	
	
	return 0;
}