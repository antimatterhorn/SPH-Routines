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

#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)

typedef struct
{
    double r;
    double rho;
	double temp;
	double u;
    double v;
	double He4;
	double C12;
	double O16;
	double Ne20;
	double Mg24;
	double Si28;
	double S32;
	double Ar36;
	double Ca40;
	double Ti44;
	double Cr48;
	double Fe52;
	double Ni56;
} line;

int usage()
{
	printf("\t Creates a 1D interpolation of paticles on the x-axis.\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_1D [sdf file] [angle(radians)] [zones] {rmax(code units)}\n");
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
    double dist_in_cm               = 6.955e7;
    double mass_in_g                = 1.989e27;
    double time_in_s                = 1.0;
    double dens_in_gccm             = mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    double energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    double specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    double pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
    double maxrho,tot_mass;

    int gnobj, nobj;
    int conf;
    double ang;
    double time;
    float tpos;
    int threads;
    int i,j,k,nz;
    
    double x,y,z,l,h,sx,sy,ds;
    double vdotr,vdots,rho,temp,mass,kern,r,rmax;
    int hc,ic,jc;
    
    if(argc < 4)
    {
        usage();
        return 0;
    } else
    {
        ang     = atof(argv[2]);
        nz      = atoi(argv[3]);
    }
	
	SDF *sdfp;
	SPHbody *body;
    
    line *lineOut;
    lineOut = malloc(nz*sizeof(line));
		
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
//					"pr", offsetof(SPHbody, pr), &conf,
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
    
    for(i=0;i<nobj;i++)
    {
        x   = body[i].x;
        y   = body[i].y;
        z   = body[i].z;
        rmax= MAX(rmax,sqrt(x*x+y*y+z*z));
    }
    
    
    if(argc > 4) rmax = atof(argv[4]);
    printf("rmax[code units] = %3.2e\n",rmax);
    
    for(i=0;i<nz;i++)
	{
		lineOut[i].r	= ((double)(i+1.0)/(double)nz)*rmax;
		lineOut[i].rho	= 0;
		//lineOut[i].v	= lineOut[i].r/tpos*dist_in_cm;
		lineOut[i].v	= 0;
	}
        
    
    ds = lineOut[i].r - lineOut[i-1].r;

    for(i=0;i<nobj;i++)
    {
        x = body[i].x;
        y = body[i].y;
        z = body[i].z;
        h = body[i].h;
		mass = body[i].mass;
        
        r = sqrt(x*x+y*y+z*z);
        l = r*r - pow(x*cos(ang)+y*sin(ang),2.0); // this is actually l^2
        if(l<4*h*h)   // the line intercepts particle i
        {
            vdotr = (body[i].vx*x + body[i].vy*y + body[i].vz*z)/r;
            if(2*h<ds)
            {
                // particle is smaller than zone
                for(j=0;j<nz;j++)
                    if(lineOut[j].r > sqrt(x*x+y*y)) break;
                j = ((j>0) ? j : j++);
                vdots = body[i].vx*cos(ang) + body[i].vy*sin(ang);
                kern = pow(ds,-3.0);
                lineOut[j-1].rho += mass*kern*dens_in_gccm;
                lineOut[j-1].v   += mass*kern*dens_in_gccm*vdots*dist_in_cm;
				lineOut[j-1].temp += mass*kern*body[i].temp;
            }
            else
            {
                for(j=0;j<nz;j++)
                {
                    sx  = lineOut[j].r*cos(ang);
                    sy  = lineOut[j].r*sin(ang);
                    r   = sqrt(pow(x-sx,2.0)+pow(y-sy,2.0)+pow(z,2.0));
                    if(r<h) // inside the particle
                    {
                        vdots = body[i].vx*cos(ang) + body[i].vy*sin(ang); // v * s^
                        
                        kern = w(h,r);
                        lineOut[j].rho += mass*kern*dens_in_gccm;
                        lineOut[j].v   += mass*kern*dens_in_gccm*vdots*dist_in_cm;
						lineOut[j].temp += mass*kern*body[i].temp;
                    }
                } 
            }

        }
    }
    
    for(i=0;i<nz;i++)
        lineOut[i].r = lineOut[i].r*dist_in_cm;
    
    for(i=0;i<nz;i++)
    {
        printf("%3.2e %3.2e %3.2e %3.2e\n",lineOut[i].r,lineOut[i].rho,lineOut[i].v,lineOut[i].temp);
        tot_mass += (4.0/3.0*3.14159)*((i>0) ? (pow(lineOut[i].r,3)-pow(lineOut[i-1].r,3))*lineOut[i].rho : 
                                       pow(lineOut[i].r,3)*lineOut[i].rho);
    }
        
    printf("tot mass = %3.2e\n",tot_mass);
    
    SDFclose(sdfp);
	free(body);
    free(lineOut);
	return 0;
}