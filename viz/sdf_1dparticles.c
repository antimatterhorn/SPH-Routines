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

char csvfile[80];

int usage()
{
	printf("\t Creates a 1D interpolation of paticles on the x-axis.\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_1D [sdf file] [angle(radians)]\n");
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
    double maxrho;

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
    
    if(argc < 2)
    {
        usage();
        return 0;
    } else
    {
        ang     = atof(argv[2]);
	}
	
	SDF *sdfp;
	SPHbody *body;
		
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

	FILE *stream, *fopen();
	snprintf(csvfile, sizeof(csvfile), "%s_%3.3f.csv", argv[1],ang);
	stream = fopen(csvfile,"w");
	
	fprintf(stream,"time:  %3.4f\n",tpos);
	fprintf(stream,"angle: %3.4f\n",ang);
	
	fprintf(stream,"x,y,rho,temp,u,He4,C12,O16,Ne20,Mg24,Si28,S32,Ar36,Ca40,Ti44,Cr48,Fe52,Ni56\n");
	
    for(i=0;i<nobj;i++)
    {
        x = body[i].x;
        y = body[i].y;
        z = body[i].z;
        h = body[i].h;
		mass = body[i].mass;
        
        r = sqrt(x*x+y*y+z*z);
        l = r*r - pow(x*cos(ang)+y*sin(ang),2.0); // this is actually l^2
        if(l<h*h)   // the line intercepts particle i
        {
            x = x*dist_in_cm;
			y = y*dist_in_cm;
			
			fprintf(stream,"%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e\n",x,y,
				   body[i].rho*dens_in_gccm,
				   body[i].temp,
				   body[i].u*specenergy_in_ergperg,
				   body[i].He4,
				   body[i].C12,
				   body[i].O16,
				   body[i].Ne20,
				   body[i].Mg24,
				   body[i].Si28,
				   body[i].S32,
				   body[i].Ar36,
				   body[i].Ca40,
				   body[i].Ti44,
				   body[i].Cr48,
				   body[i].Fe52,
				   body[i].Ni56);
        }
    }
    
    fclose(stream);
    SDFclose(sdfp);
	free(body);
	return 0;
}