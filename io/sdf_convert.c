/*
 *  sdf2csv.c
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

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
double xmax,zmax,dist,x,y,z;
double h,rho,temp,mass,kern,r,xc,zc;
int hc,ic,jc;
double time;
float tpos,norm;
int choice,pixels,threads;
int i,j,k,b;
char sdffile[80];
char asciifile[80];

double dist_in_cm;
double mass_in_g;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;
double maxrho,maxtemp,abund;

typedef struct {
	double x, y, z;             /* position of body */							//3
	float mass;           /* mass of body */
	float vx, vy, vz;     /* velocity of body */								//4
	float u;              /* specific energy of body*/
	float h;              /* smoothing length of body */						//6
	float rho;            /* density of body */
	float pr;            /* pressure of body */
	float temp;           /* temperature of body */								//11
	float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */		//18
	float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */			//24
	float vsound;
	float abar;           /* avg number of nucleons per particle of body */
	float zbar;           /* avg number of protons per particle of body */
	float ax, ay, az;     /* acceleration of body */							//30
	float bdot;
	float kappa;			/* opacity */
	float sigma;			/* conductivity */									//38
} SPHbody;

typedef struct {
    float x, y, z;        /* position of cell */	
	float vx, vy, vz;       /* velocity of cell */
    float rho;          /* density of cell */
	float u;         /* temperature of cell */	
    float He4;
    float C12;
    float O16;
} SPHsmallbody;

int usage()
{
	printf("\t Converts an SDF file to other formats in cgi units.\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_convert [sdf file] {choice(1-3)}\n");
	return 0;
}

int struct_cmp_by_pos(const void *a, const void *b) 
{ 
    const SPHbody *ia = (const SPHbody *)a;
    const SPHbody *ib = (const SPHbody *)b;
    
	/* return 1 or -1 if z is not equal */
	if (ia->z > ib->z) return 1;
	if (ia->z < ib->z) return -1;
	
	/* return 1 or -1 if y is not equal */
	if (ia->y > ib->y) return 1;
	if (ia->y < ib->y) return -1;
	
	/* return 1 or -1 if x is not equal */
	if (ia->x > ib->x) return 1;
	if (ia->x < ib->x) return -1;
	
	/* return 0 if all are equal */
	return 0;
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
		usage();
		return 0;
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
	
	
    

	for(i=0;i<nobj;i++)
	{
		body[i].x       = body[i].x*dist_in_cm;
		body[i].y       = body[i].y*dist_in_cm;
		body[i].z       = body[i].z*dist_in_cm;
        body[i].h       = body[i].h*dist_in_cm;
		body[i].vx      = body[i].vx*dist_in_cm;
		body[i].vy      = body[i].vy*dist_in_cm;
		body[i].vz      = body[i].vz*dist_in_cm;
        body[i].mass    = body[i].mass*mass_in_g;
		body[i].rho     = body[i].rho*dens_in_gccm;
		body[i].temp    = body[i].temp;
		body[i].u       = body[i].u*specenergy_in_ergperg;
		body[i].pr      = body[i].pr*pressure_in_ergperccm;
		body[i].vsound  = body[i].vsound*dist_in_cm;
		body[i].He4     = body[i].He4;
		body[i].C12     = body[i].C12;
		body[i].O16     = body[i].O16;
		body[i].Ne20    = body[i].Ne20;
		body[i].Mg24    = body[i].Mg24;
		body[i].Si28    = body[i].Si28;
		body[i].S32     = body[i].S32;
		body[i].Ar36    = body[i].Ar36;
		body[i].Ca40    = body[i].Ca40;
		body[i].Ti44    = body[i].Ti44;
		body[i].Cr48    = body[i].Cr48;
		body[i].Fe52    = body[i].Fe52;
		body[i].Ni56    = body[i].Ni56;
	}
	
	printf("Sorting the array.\n");
	
	qsort(body,nobj,sizeof(SPHbody),struct_cmp_by_pos);
	
    if(argc < 3)
    {
        printf("Choose an output type:\n");
        printf("\t1. 3D Binary\n");
        printf("\t2. 3D Binary (reduced)\n");
        printf("\t3. 3D Ascii\n");
        printf("\t4. Quit\n");
        printf(":> ");
        scanf("%d",&choice);
    }
    else choice = atoi(argv[2]);
    
    
    FILE *stream, *fopen();
    
    switch (choice) {
        case 1:
        {
            snprintf(asciifile, sizeof(asciifile), "%s.dat", argv[1]);
            stream = fopen(asciifile,"wb");
            fwrite(body,sizeof(SPHbody),nobj,stream);
            fclose(stream);
            break;
        }
        case 2:
        {
            SPHsmallbody *smallbody;
            smallbody = malloc(nobj*sizeof(SPHsmallbody));
            for(i=0;i<nobj;i++)
            {
                smallbody[i].x   = body[i].x;
                smallbody[i].y   = body[i].y;
                smallbody[i].z   = body[i].z;
                smallbody[i].vx  = body[i].vx;
                smallbody[i].vy  = body[i].vy;
                smallbody[i].vz  = body[i].vz;
                smallbody[i].rho = ((body[i].rho>0) ? body[i].rho : 1e-10);
                smallbody[i].u   = body[i].u; 
                
                abund = body[i].He4 + body[i].C12 + body[i].O16;
                body[i].He4 = ((body[i].rho>1e-10) ? body[i].He4 / abund : 0.0);
                body[i].C12 = ((body[i].rho>1e-10) ? body[i].C12 / abund : 0.5);
                body[i].O16 = ((body[i].rho>1e-10) ? body[i].O16 / abund : 0.5);
                
                smallbody[i].He4 = body[i].He4;
                smallbody[i].C12 = body[i].C12;
                smallbody[i].O16 = body[i].O16;
            }
            snprintf(asciifile, sizeof(asciifile), "%s.dat", argv[1]);
            stream = fopen(asciifile,"wb");
            fwrite(smallbody,sizeof(SPHsmallbody),nobj,stream);
            fclose(stream);
            free(smallbody);
            break;
        }
        case 3:
        {
            snprintf(asciifile, sizeof(asciifile), "%s.ascii", argv[1]);
            stream = fopen(asciifile,"w");
            
            fprintf(stream,"x(cm) y(cm) z(cm) h(cm) mass(g) rho(g/cc) pr(erg/cc) u(erg/g) vx(cm/s) vy(cm/s) vz(cm/s) T(K) He4 C12 O16 Ne20 Mg24 Si28 S32 Ar36 Ca40 Ti44 Cr48 Fe52 Ni56\n");
            for(i=0;i<nobj;i++)
            {
                fprintf(stream,"%3.3e ",body[i].x);
                fprintf(stream,"%3.3e ",body[i].y);
                fprintf(stream,"%3.3e ",body[i].z);
                fprintf(stream,"%3.3e ",body[i].h);
                fprintf(stream,"%3.3e ",body[i].mass);
                fprintf(stream,"%3.3e ",body[i].rho);
                fprintf(stream,"%3.3e ",body[i].pr);
                fprintf(stream,"%3.3e ",body[i].u);
                fprintf(stream,"%3.3e ",body[i].vx);
                fprintf(stream,"%3.3e ",body[i].vy);
                fprintf(stream,"%3.3e ",body[i].vz);
                fprintf(stream,"%3.3e ",body[i].temp);
                fprintf(stream,"%1.5f ",body[i].He4);
                fprintf(stream,"%1.5f ",body[i].C12);
                fprintf(stream,"%1.5f ",body[i].O16);
                fprintf(stream,"%1.5f ",body[i].Ne20);
                fprintf(stream,"%1.5f ",body[i].Mg24);
                fprintf(stream,"%1.5f ",body[i].Si28);
                fprintf(stream,"%1.5f ",body[i].S32);
                fprintf(stream,"%1.5f ",body[i].Ar36);
                fprintf(stream,"%1.5f ",body[i].Ca40);
                fprintf(stream,"%1.5f ",body[i].Ti44);
                fprintf(stream,"%1.5f ",body[i].Cr48);
                fprintf(stream,"%1.5f ",body[i].Fe52);
                fprintf(stream,"%1.5f\n",body[i].Ni56);
            }
            
            fclose(stream);
            break;
        }  
        default:
            break;
    }
    
	

	free(body);
    

	return 0;
}
