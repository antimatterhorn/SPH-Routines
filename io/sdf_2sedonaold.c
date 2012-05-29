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
double xmax,zmax;
double time;
float tpos;
int pixels,threads;
int i,j,k,ic,jc,hc;
float xc,zc;
float r,kern,h,mass,x,z,y,rho,temp,vx,vy,He4,C12,O16,Ne20,Mg24,Si28,S32,Ar36,Ca40,Ti44,Cr48,Fe52,Ni56;
float norm;
char sdffile[80];
char asciifile[80];

double dist_in_cm;
double mass_in_g;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;

int bj=128,bi=256;



typedef struct {
	double x, z;             /* position of body */							//3
	float vx, vy;     /* velocity of body */								//4
	float rho;            /* density of body */
	float temp;           /* temperature of body */								//11
	float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */		//18
	float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */			//24
} asciibody;



int usage()
{
	printf("\t Creates...\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_2sedona [sdf file] [xmax(code units)]\n");
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
		xmax = atof(argv[2]);

	}
	
	zmax = xmax;
	
	pixels = 256;
		
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
					NULL);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
	
	asciibody ascii[bi*bi];
	
	for(i=0;i<bi;i++)
	{
		for(j=0;j<bi;j++)
		{
			ascii[i*bi+j].x = ((float)j/(float)bi)*2.0*xmax-xmax;
			ascii[i*bi+j].z = zmax-((float)i/(float)bi)*2.0*zmax;
		}
	}
	
	for(k=0;k<nobj;k++)
	{
		if (fabs(body[k].y) < 2.0*body[k].h) /* midplane */
		{
			x = body[k].x;
			y = body[k].y;
			z = body[k].z;
			h = body[k].h;
			hc = 2*h*bi/(2*zmax);
			mass = body[k].mass;
			
			jc = x*(bi/(2.0*xmax)) + bi/2.0;
			ic = -z*(bi/(2.0*zmax)) + bi/2.0;
			
			
			
			
			if(hc<1){ /* particle is smaller than a pixel, fill pixel instead */
				kern = pow(2.0*xmax/pixels,-3);
				ascii[ic*bi+jc].rho		+= mass*kern*dens_in_gccm;
				ascii[ic*bi+jc].temp	+= mass*kern*dens_in_gccm*body[k].temp;
				ascii[ic*bi+jc].vx		+= mass*kern*dens_in_gccm*body[k].vx;
				ascii[ic*bi+jc].vy		+= mass*kern*dens_in_gccm*body[k].vy;
				ascii[ic*bi+jc].He4		+= mass*kern*dens_in_gccm*body[k].He4;
				ascii[ic*bi+jc].C12		+= mass*kern*dens_in_gccm*body[k].C12;
				ascii[ic*bi+jc].O16		+= mass*kern*dens_in_gccm*body[k].O16;
				ascii[ic*bi+jc].Ne20	+= mass*kern*dens_in_gccm*body[k].Ne20;
				ascii[ic*bi+jc].Mg24	+= mass*kern*dens_in_gccm*body[k].Mg24;
				ascii[ic*bi+jc].Si28	+= mass*kern*dens_in_gccm*body[k].Si28;
				ascii[ic*bi+jc].S32		+= mass*kern*dens_in_gccm*body[k].S32;
				ascii[ic*bi+jc].Ar36	+= mass*kern*dens_in_gccm*body[k].Ar36;
				ascii[ic*bi+jc].Ca40	+= mass*kern*dens_in_gccm*body[k].Ca40;
				ascii[ic*bi+jc].Ti44	+= mass*kern*dens_in_gccm*body[k].Ti44;
				ascii[ic*bi+jc].Cr48	+= mass*kern*dens_in_gccm*body[k].Cr48;
				ascii[ic*bi+jc].Fe52	+= mass*kern*dens_in_gccm*body[k].Fe52;
				ascii[ic*bi+jc].Ni56	+= mass*kern*dens_in_gccm*body[k].Ni56;
			}else{
				for(i=ic-hc;i<ic+hc+1;i++)
				{
					for(j=jc-hc;j<jc+hc+1;j++)
					{
						if(j>-1 && i>-1 && i<bi && j<bi){
							xc = (double)j/(double)bi*(2.0*xmax) - xmax;
							zc = zmax - (double)i/(double)pixels*(2.0*zmax);
							
							r = sqrt(pow(xc-x,2.0)+pow(zc-z,2.0)+y*y);
							kern = w(h,r);
							ascii[i*bi+j].rho	+= mass*kern*dens_in_gccm;
							ascii[i*bi+j].temp	+= mass*kern*dens_in_gccm*body[k].temp;
							ascii[i*bi+j].vx	+= mass*kern*dens_in_gccm*body[k].vx;
							ascii[i*bi+j].vy	+= mass*kern*dens_in_gccm*body[k].vy;
							ascii[i*bi+j].He4	+= mass*kern*dens_in_gccm*body[k].He4;
							ascii[i*bi+j].C12	+= mass*kern*dens_in_gccm*body[k].C12;
							ascii[i*bi+j].O16	+= mass*kern*dens_in_gccm*body[k].O16;
							ascii[i*bi+j].Ne20	+= mass*kern*dens_in_gccm*body[k].Ne20;
							ascii[i*bi+j].Mg24	+= mass*kern*dens_in_gccm*body[k].Mg24;
							ascii[i*bi+j].Si28	+= mass*kern*dens_in_gccm*body[k].Si28;
							ascii[i*bi+j].S32	+= mass*kern*dens_in_gccm*body[k].S32;
							ascii[i*bi+j].Ar36	+= mass*kern*dens_in_gccm*body[k].Ar36;
							ascii[i*bi+j].Ca40	+= mass*kern*dens_in_gccm*body[k].Ca40;
							ascii[i*bi+j].Ti44	+= mass*kern*dens_in_gccm*body[k].Ti44;
							ascii[i*bi+j].Cr48	+= mass*kern*dens_in_gccm*body[k].Cr48;
							ascii[i*bi+j].Fe52	+= mass*kern*dens_in_gccm*body[k].Fe52;
							ascii[i*bi+j].Ni56	+= mass*kern*dens_in_gccm*body[k].Ni56;
						}
					}
				}
			}	
		}
	}
	
	for(i=0;i<bi;i++)
	{
		for(j=0;j<bi;j++)
		{
			norm = 0;
			
			ascii[i*bi+j].x *= dist_in_cm;
			ascii[i*bi+j].z *= dist_in_cm;
			
			ascii[i*bi+j].temp	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].temp/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].vx	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].vx/ascii[i*bi+j].rho : 0)*dist_in_cm;
			ascii[i*bi+j].vy	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].vy/ascii[i*bi+j].rho : 0)*dist_in_cm;
			ascii[i*bi+j].He4	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].He4/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].C12	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].C12/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].O16	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].O16/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].Ne20	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].Ne20/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].Mg24	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].Mg24/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].Si28	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].Si28/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].S32	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].S32/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].Ar36	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].Ar36/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].Ca40	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].Ca40/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].Ti44	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].Ti44/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].Cr48	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].Cr48/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].Fe52	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].Fe52/ascii[i*bi+j].rho : 0);
			ascii[i*bi+j].Ni56	= ((ascii[i*bi+j].rho>0) ? ascii[i*bi+j].Ni56/ascii[i*bi+j].rho : 0);
			
			norm += ascii[i*bi+j].He4;
			norm += ascii[i*bi+j].C12;
			norm += ascii[i*bi+j].O16;
			norm += ascii[i*bi+j].Ne20;
			norm += ascii[i*bi+j].Mg24;
			norm += ascii[i*bi+j].Si28;
			norm += ascii[i*bi+j].S32;
			norm += ascii[i*bi+j].Ar36;
			norm += ascii[i*bi+j].Ca40;
			norm += ascii[i*bi+j].Ti44;
			norm += ascii[i*bi+j].Cr48;
			norm += ascii[i*bi+j].Fe52;
			norm += ascii[i*bi+j].Ni56;
			
			if(norm==0) norm=1;
			
			ascii[i*bi+j].He4 = ascii[i*bi+j].He4/norm;
			ascii[i*bi+j].C12 = ascii[i*bi+j].C12/norm;
			ascii[i*bi+j].O16 = ascii[i*bi+j].O16/norm;
			ascii[i*bi+j].Ne20 = ascii[i*bi+j].Ne20/norm;
			ascii[i*bi+j].Mg24 = ascii[i*bi+j].Mg24/norm;
			ascii[i*bi+j].Si28 = ascii[i*bi+j].Si28/norm;
			ascii[i*bi+j].S32 = ascii[i*bi+j].S32/norm;
			ascii[i*bi+j].Ar36 = ascii[i*bi+j].Ar36/norm;
			ascii[i*bi+j].Ca40 = ascii[i*bi+j].Ca40/norm;
			ascii[i*bi+j].Ti44 = ascii[i*bi+j].Ti44/norm;
			ascii[i*bi+j].Cr48 = ascii[i*bi+j].Cr48/norm;
			ascii[i*bi+j].Fe52 = ascii[i*bi+j].Fe52/norm;
			ascii[i*bi+j].Ni56 = ascii[i*bi+j].Ni56/norm;
		}
	}

	snprintf(asciifile, sizeof(asciifile), "%s.ascii", argv[1]);
	FILE *stream, *fopen();
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(asciifile,"w");
	/*
	fprintf(stream,"**** header **************\n");
	fprintf(stream,"box size:\t %d x %d\n",bj,bi);
	fprintf(stream,"dimensions:\t %3.2ecm x %3.2ecm\n",xmax,xmax*2.0);
	fprintf(stream,"t-pos:\t %3.2f\n",tpos);
	fprintf(stream,"**************************\n");
	fprintf(stream,"\n");
	 */
	fprintf(stream,"x(cm) z(cm) rho(g/cc) vx(cm/s) vy(cm/s) T(K) He4 C12 O16 Ne20 Mg24 Si28 S32 Ar36 Ca40 Ti44 Cr48 Fe52 Ni56\n");
	
	
	for(i=0;i<bi;i++)
	{
		for(j=0;j<bi;j++)
		{
			if((ascii[i*bi+j].x >= 0) && (abs(ascii[i*bi+j].z) <=1.0))
			{
				fprintf(stream,"%03.2e %03.2e ",
						ascii[i*bi+j].x,ascii[i*bi+j].z);
				if(ascii[i*bi+j].rho > 1e-8)
				{
					fprintf(stream,"%03.2e %03.2e %03.2e %03.2e ",
							ascii[i*bi+j].rho,ascii[i*bi+j].vx,ascii[i*bi+j].vy,ascii[i*bi+j].temp);
					
					fprintf(stream,"%1.5f ",ascii[i*bi+j].He4);
					fprintf(stream,"%1.5f ",ascii[i*bi+j].C12);
					fprintf(stream,"%1.5f ",ascii[i*bi+j].O16);
					fprintf(stream,"%1.5f ",ascii[i*bi+j].Ne20);
					fprintf(stream,"%1.5f ",ascii[i*bi+j].Mg24);
					fprintf(stream,"%1.5f ",ascii[i*bi+j].Si28);
					fprintf(stream,"%1.5f ",ascii[i*bi+j].S32);
					fprintf(stream,"%1.5f ",ascii[i*bi+j].Ar36);
					fprintf(stream,"%1.5f ",ascii[i*bi+j].Ca40);
					fprintf(stream,"%1.5f ",ascii[i*bi+j].Ti44);
					fprintf(stream,"%1.5f ",ascii[i*bi+j].Cr48);
					fprintf(stream,"%1.5f ",ascii[i*bi+j].Fe52);
					fprintf(stream,"%1.5f\n",ascii[i*bi+j].Ni56);
				}
				else 
				{
					fprintf(stream,"%03.2e %03.2e %03.2e %03.2e ",0.0,0.0,0.0,0.0);
					fprintf(stream,"%1.5f ",0.0);
					fprintf(stream,"%1.5f ",0.0);
					fprintf(stream,"%1.5f ",0.0);
					fprintf(stream,"%1.5f ",0.0);
					fprintf(stream,"%1.5f ",0.0);
					fprintf(stream,"%1.5f ",0.0);		
					fprintf(stream,"%1.5f ",0.0);
					fprintf(stream,"%1.5f ",0.0);		
					fprintf(stream,"%1.5f ",0.0);
					fprintf(stream,"%1.5f ",0.0);		
					fprintf(stream,"%1.5f ",0.0);
					fprintf(stream,"%1.5f ",0.0);
					fprintf(stream,"%1.5f\n",0.0);
				}


				
			}
		}
	}
	
	fclose(stream);
	
	
	free(body);

	
	
	return 0;
}