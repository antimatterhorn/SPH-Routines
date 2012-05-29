/*
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


int gnobj, nobj;
int conf;
double xc,yc,zc,hc;
float xmax,ymax,zmax,x,y,z,rr,kern,norm;
float rho,temp,mass,h,tpos;
int choice,pixels;
int i,j,k,l,imi,imj,imk,ic,jc,kc,r,ind;
char sdffile[80];
char csvfile[80];

double dist_in_cm;
double mass_in_g;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;

void x2i()
{
	ic = x*(pixels/(2.0*xmax)) + pixels/2.0;	
}

void y2j()
{
	jc = y*(pixels/(2.0*ymax)) + pixels/2.0;
}

void z2k()
{
	kc = z*(pixels/(2.0*zmax)) + pixels/2.0;
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
    
    
    
    SDF *sdfp;
	SPHbody *body;
	
	if (argc < 2){
		printf("SDF file: ");
		gets (argv[1]);
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
	
    if (argc < 3){
		printf("xmax:");
		scanf("%f", &xmax);
	}
	else {
		xmax = atof(argv[2]);
	}
    
    if (argc < 4){
		pixels = 50;
	}
	else {
		pixels = atoi(argv[3]);
	}
	


	ymax = zmax = xmax;
	
    int items = 13+8;
	double ***(outArray[items]);
	
    for(ind=0;ind<items;ind++)
    {
        outArray[ind] = (double ***) malloc(pixels*sizeof(double **));
        for(i=0;i<pixels;i++)
        {
            outArray[ind][i] = (double **) malloc(pixels*sizeof(double *));
            for(j=0;j<pixels;j++)
            {
                outArray[ind][i][j] = (double *) malloc(pixels*sizeof(double));
                for(k=0;k<pixels;k++)
                {
                    outArray[ind][i][j][k] = 0;
                }
            }
        }
    }
	
    /*
     0	x
     1	y
     2	z
     3	rho
     4	temp
     5	vx
     6	vy
     7	vz
     8	He4
     9	C12
     10	O16
     11	Ne20
     12	Mg24
     13	Si28
     14	S32
     15	Ar36
     16	Ca40
     17	Ti44
     18	Cr48
     19	Fe52
     20	Ni56
    */
    
    
    
	SPHbody *p;
	
	printf("Anchoring to grid...\n");
	
	for(p = body; p < body+gnobj; p++)
	{
        //printf("%d / %d\n",p,nobj-1);
        
		x = p->x;
		y = p->y;
		z = p->z;
		h = p->h;
		rho = p->rho;
		temp = p->temp;
		mass = p->mass;
		
		x2i();
		y2j();
		z2k();

		r = h*(pixels/(2.0*xmax));
		
		if (r==0)
		{
			//printf("found that r was too small, filling only 1 block\n");
			if (ic >= 0 && ic < pixels && jc >= 0 && jc < pixels && kc >= 0 && kc < pixels)
			{
				kern = pow(2.0*xmax/pixels,-3);
                i = ic;
                j = jc;
                k = kc;
                outArray[3][i][j][k] += mass*kern*dens_in_gccm;
                outArray[4][i][j][k] += mass*kern*dens_in_gccm*p->temp;
                outArray[5][i][j][k] += mass*kern*dens_in_gccm*p->vx;
                outArray[6][i][j][k] += mass*kern*dens_in_gccm*p->vy;
                outArray[7][i][j][k] += mass*kern*dens_in_gccm*p->vz;
                outArray[8][i][j][k] += mass*kern*dens_in_gccm*p->He4;
                outArray[9][i][j][k] += mass*kern*dens_in_gccm*p->C12;
                outArray[10][i][j][k] += mass*kern*dens_in_gccm*p->O16;
                outArray[11][i][j][k] += mass*kern*dens_in_gccm*p->Ne20;
                outArray[12][i][j][k] += mass*kern*dens_in_gccm*p->Mg24;
                outArray[13][i][j][k] += mass*kern*dens_in_gccm*p->Si28;
                outArray[14][i][j][k] += mass*kern*dens_in_gccm*p->S32;
                outArray[15][i][j][k] += mass*kern*dens_in_gccm*p->Ar36;
                outArray[16][i][j][k] += mass*kern*dens_in_gccm*p->Ca40;
                outArray[17][i][j][k] += mass*kern*dens_in_gccm*p->Ti44;
                outArray[18][i][j][k] += mass*kern*dens_in_gccm*p->Cr48;
                outArray[19][i][j][k] += mass*kern*dens_in_gccm*p->Fe52;
                outArray[20][i][j][k] += mass*kern*dens_in_gccm*p->Ni56;
			}
			
		}
		else
		{
			for(i = ic-2*r; i < ic+2*r+1; i++)
			{
				for(j = jc-2*r; j < jc+2*r+1; j++)
				{
					for (k = kc-2*r; k < kc+2*r+1; k++) {
						if (i >= 0 && i < pixels && j >= 0 && j < pixels && k >= 0 && k < pixels) // inside the image
						{
							xc = (double)i/(double)pixels*(2.0*xmax) - xmax;
                            yc = (double)j/(double)pixels*(2.0*ymax) - ymax;
                            zc = (double)k/(double)pixels*(2.0*ymax) - ymax;
                            rr = sqrt(pow((xc-x),2.0)+pow((yc-y),2.0)+pow((zc-z),2.0));
							kern = w(h,rr);
                            if (rr <= 2*r) // inside the circle
							{							
								outArray[3][i][j][k] += mass*kern*dens_in_gccm;
                                outArray[4][i][j][k] += mass*kern*dens_in_gccm*p->temp;
                                outArray[5][i][j][k] += mass*kern*dens_in_gccm*p->vx;
                                outArray[6][i][j][k] += mass*kern*dens_in_gccm*p->vy;
                                outArray[7][i][j][k] += mass*kern*dens_in_gccm*p->vz;
                                outArray[8][i][j][k] += mass*kern*dens_in_gccm*p->He4;
                                outArray[9][i][j][k] += mass*kern*dens_in_gccm*p->C12;
                                outArray[10][i][j][k] += mass*kern*dens_in_gccm*p->O16;
                                outArray[11][i][j][k] += mass*kern*dens_in_gccm*p->Ne20;
                                outArray[12][i][j][k] += mass*kern*dens_in_gccm*p->Mg24;
                                outArray[13][i][j][k] += mass*kern*dens_in_gccm*p->Si28;
                                outArray[14][i][j][k] += mass*kern*dens_in_gccm*p->S32;
                                outArray[15][i][j][k] += mass*kern*dens_in_gccm*p->Ar36;
                                outArray[16][i][j][k] += mass*kern*dens_in_gccm*p->Ca40;
                                outArray[17][i][j][k] += mass*kern*dens_in_gccm*p->Ti44;
                                outArray[18][i][j][k] += mass*kern*dens_in_gccm*p->Cr48;
                                outArray[19][i][j][k] += mass*kern*dens_in_gccm*p->Fe52;
                                outArray[20][i][j][k] += mass*kern*dens_in_gccm*p->Ni56;
							}
						}
					}
				}
			}	
		}
	}
	
	printf("Writing to file...\n");
	
	//open the stream file
	snprintf(csvfile, sizeof(csvfile), "%f.cube", tpos);
	FILE *stream, *fopen();
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");
	
    fprintf(stream,"x(cm) y(cm) z(cm) rho(g/cc) temp(K) vx(cm/s) vy(cm/s) vz(cm/s) He4 C12 O16 Ne20 Mg24 Si28 S32 Ar36 Ca40 Ti44 Cr48 Fe52 Ni56\n");
    
    for(i=0;i<pixels;i++)
        for(j=0;j<pixels;j++)
            for(k=0;k<pixels;k++)
            {
                x = (double)i/(double)pixels*(2.0*xmax) - xmax;
                y = (double)j/(double)pixels*(2.0*ymax) - ymax;
                z = (double)k/(double)pixels*(2.0*ymax) - ymax;
                x = x*dist_in_cm;
                y = y*dist_in_cm;
                z = z*dist_in_cm;
                for(ind=4;ind<items;ind++)
                    outArray[ind][i][j][k] = ((outArray[3][i][j][k] > 0) ? (outArray[ind][i][j][k]/outArray[3][i][j][k]) : 1e-12);
                norm = 0.0;
                for(ind=8;ind<21;ind++)
                    norm += outArray[ind][i][j][k];
                if(norm > 1e-10)
                    for(ind=8;ind<21;ind++)
                        outArray[ind][i][j][k] = outArray[ind][i][j][k]/norm;
                fprintf(stream,"%3.3e %3.3e %3.3e %3.3e %3.3e ",x,y,z,((outArray[3][i][j][k]>0) ? outArray[3][i][j][k] : 1e-12),outArray[4][i][j][k]);
                for(ind=5;ind<8;ind++)
                    fprintf(stream,"%3.3e ",outArray[ind][i][j][k]*dist_in_cm);
                for(ind=8;ind<21;ind++)
                    fprintf(stream,"%3.4f ",outArray[ind][i][j][k]);
                fprintf(stream,"\n");
                
            }


	//close the stream file
	fclose(stream);
    for(ind=0;ind<items;ind++)
        free(outArray[ind]);
    free(body);

	
	return 0;
}