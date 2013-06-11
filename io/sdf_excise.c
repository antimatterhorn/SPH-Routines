/*
 *  sdf_masscutoff.c
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  reads an sdf file and outputs a subset with a chosen mass
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


char abovefile[80];
char belowfile[80];

double dist_in_cm;
double mass_in_g;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;

int nabove,i,j,k;
double rhomax;

double ke,pe,te;
float gam,tvel,Gnewt,eps,dt,tpos,tolerance,frac_tolerance;
int ndim,iter;


int usage()
{
	printf("\t Creates subset files with particles excised above a chosen density.\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_excise [sdf file] [rhomax(cgs)]\n");
	return 0;
}

int main(int argc, char **argv[])
{
	    
    if(argc < 3)
    {
        usage();
        return 0;
    }
    
    time_in_s				= 1;
	dist_in_cm				= 6.955e7;
	mass_in_g				= 1.989e27;
	dens_in_gccm			= mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
    
    rhomax = atof(argv[2]);
    
    printf("Excising %s around %3.2e (%3.3f).\n",argv[1],rhomax,rhomax/dens_in_gccm);
    
    rhomax = rhomax / dens_in_gccm;
    
    SDF *sdfp;
	SPHbody *body, *above, *below;
		
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
					"abar", offsetof(SPHbody, abar), &conf,
					"zbar", offsetof(SPHbody, zbar), &conf,
					"ax", offsetof(SPHbody, ax), &conf,
					"ay", offsetof(SPHbody, ay), &conf,
					"az", offsetof(SPHbody, az), &conf,
					"lax", offsetof(SPHbody, lax), &conf,
					"lay", offsetof(SPHbody, lay), &conf,
					"laz", offsetof(SPHbody, laz), &conf,
//					"gax", offsetof(SPHbody, gax), &conf,
//					"gay", offsetof(SPHbody, gay), &conf,
//					"gaz", offsetof(SPHbody, gaz), &conf,
//					"grav_mass", offsetof(SPHbody, grav_mass), &conf,
					"phi", offsetof(SPHbody, phi), &conf,
//					"tacc", offsetof(SPHbody, tacc), &conf,
					"idt", offsetof(SPHbody, idt), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					NULL);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
    SDFgetfloatOrDefault(sdfp, "dt",  &dt, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "eps",  &eps, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "Gnewt",  &Gnewt, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "tolerance",  &tolerance, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "frac_tolerance",  &frac_tolerance, (float)0.0);
	//SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "tvel",  &tvel, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "gamma",  &gam, (float)0.0);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
    
    for(i=0;i<nobj;i++) if(body[i].rho >= rhomax) nabove++;
    
    printf("%d particles will go into above_%s\n",nabove,argv[1]);
    printf("%d particles will go into below_%s\n",nobj-nabove,argv[1]);
    
    snprintf(abovefile, sizeof(abovefile), "above_%s", argv[1]);
    snprintf(belowfile, sizeof(belowfile), "below_%s", argv[1]);
    
    above = malloc(nabove*sizeof(SPHbody));
    below = malloc((nobj-nabove)*sizeof(SPHbody));
    
    for(i=0;i<nobj;i++)
    {
        if(body[i].rho >= rhomax)
        {
            above[j].x      = body[i].x;
            above[j].y      = body[i].y;
            above[j].z      = body[i].z;
            above[j].mass   = body[i].mass;
            above[j].vx     = body[i].vx;
            above[j].vy     = body[i].vy;
            above[j].vz     = body[i].vz;
            above[j].u      = body[i].u;
            above[j].h      = body[i].h;
            above[j].rho    = body[i].rho;
            above[j].pr     = body[i].pr;
            above[j].drho_dt = body[i].drho_dt;
            above[j].udot   = body[i].udot;
            above[j].temp   = body[i].temp;
            above[j].He4    = body[i].He4;
            above[j].C12    = body[i].C12;
            above[j].O16    = body[i].O16;
            above[j].Ne20   = body[i].Ne20;
            above[j].Mg24   = body[i].Mg24;
            above[j].Si28   = body[i].Si28;
            above[j].S32    = body[i].S32;
            above[j].Ar36   = body[i].Ar36;
            above[j].Ca40   = body[i].Ca40;
            above[j].Ti44   = body[i].Ti44;
            above[j].Cr48   = body[i].Cr48;
            above[j].Fe52   = body[i].Fe52;
            above[j].Ni56   = body[i].Ni56;
            above[j].abar   = body[i].abar;
            above[j].zbar   = body[i].zbar;
            above[j].ax     = body[i].ax;
            above[j].ay     = body[i].ay;
            above[j].az     = body[i].az;
            above[j].lax    = body[i].lax;
            above[j].lay    = body[i].lay;
            above[j].laz    = body[i].laz;
            above[j].phi    = body[i].phi;
            above[j].idt    = body[i].idt;
            above[j].nbrs   = body[i].nbrs;
            above[j].ident  = body[i].ident;
            above[j].windid = body[i].windid;
            j++;
        }
        else
        {
            below[k].x      = body[i].x;
            below[k].y      = body[i].y;
            below[k].z      = body[i].z;
            below[k].mass   = body[i].mass;
            below[k].vx     = body[i].vx;
            below[k].vy     = body[i].vy;
            below[k].vz     = body[i].vz;
            below[k].u      = body[i].u;
            below[k].h      = body[i].h;
            below[k].rho    = body[i].rho;
            below[k].pr     = body[i].pr;
            below[k].drho_dt = body[i].drho_dt;
            below[k].udot   = body[i].udot;
            below[k].temp   = body[i].temp;
            below[k].He4    = body[i].He4;
            below[k].C12    = body[i].C12;
            below[k].O16    = body[i].O16;
            below[k].Ne20   = body[i].Ne20;
            below[k].Mg24   = body[i].Mg24;
            below[k].Si28   = body[i].Si28;
            below[k].S32    = body[i].S32;
            below[k].Ar36   = body[i].Ar36;
            below[k].Ca40   = body[i].Ca40;
            below[k].Ti44   = body[i].Ti44;
            below[k].Cr48   = body[i].Cr48;
            below[k].Fe52   = body[i].Fe52;
            below[k].Ni56   = body[i].Ni56;
            below[k].abar   = body[i].abar;
            below[k].zbar   = body[i].zbar;
            below[k].ax     = body[i].ax;
            below[k].ay     = body[i].ay;
            below[k].az     = body[i].az;
            below[k].lax    = body[i].lax;
            below[k].lay    = body[i].lay;
            below[k].laz    = body[i].laz;
            below[k].phi    = body[i].phi;
            below[k].idt    = body[i].idt;
            below[k].nbrs   = body[i].nbrs;
            below[k].ident  = body[i].ident;
            below[k].windid = body[i].windid;
            k++;
        }
    }
    
    iter=tpos=0;
    SDFwrite(abovefile, nabove, 
			 nabove, above, sizeof(SPHbody),
			 SPHOUTBODYDESC,
			 "npart", SDF_INT, nabove,
			 "iter", SDF_INT, iter,
			 "dt", SDF_FLOAT, dt,
			 "eps", SDF_FLOAT, eps,
			 "Gnewt", SDF_FLOAT, Gnewt,
			 "tolerance", SDF_FLOAT, tolerance,
			 "frac_tolerance", SDF_FLOAT, frac_tolerance,
			 "ndim", SDF_INT, ndim,
			 "tpos", SDF_FLOAT, tpos,
			 "tvel", SDF_FLOAT, tvel,
			 //"t_wind", SDF_FLOAT, twind_out,
			 //"R0", SDF_FLOAT, output_R0,
			 //"Omega0", SDF_FLOAT, cosmo.Omega0,
			 //"H0", SDF_FLOAT, cosmo.H0,
			 //"Lambda_prime", SDF_FLOAT, cosmo.Lambda,
			 //"hubble", SDF_FLOAT, output_h,
			 //"redshift", SDF_FLOAT, output_z,
			 "gamma", SDF_FLOAT, gam,
			 //"centmass", SDF_FLOAT, centmass, 
			 NULL);
    SDFwrite(belowfile, nobj-nabove, 
			 nobj-nabove, below, sizeof(SPHbody),
			 SPHOUTBODYDESC,
			 "npart", SDF_INT, nobj-nabove,
			 "iter", SDF_INT, iter,
			 "dt", SDF_FLOAT, dt,
			 "eps", SDF_FLOAT, eps,
			 "Gnewt", SDF_FLOAT, Gnewt,
			 "tolerance", SDF_FLOAT, tolerance,
			 "frac_tolerance", SDF_FLOAT, frac_tolerance,
			 "ndim", SDF_INT, ndim,
			 "tpos", SDF_FLOAT, tpos,
			 "tvel", SDF_FLOAT, tvel,
			 //"t_wind", SDF_FLOAT, twind_out,
			 //"R0", SDF_FLOAT, output_R0,
			 //"Omega0", SDF_FLOAT, cosmo.Omega0,
			 //"H0", SDF_FLOAT, cosmo.H0,
			 //"Lambda_prime", SDF_FLOAT, cosmo.Lambda,
			 //"hubble", SDF_FLOAT, output_h,
			 //"redshift", SDF_FLOAT, output_z,
			 "gamma", SDF_FLOAT, gam,
			 //"centmass", SDF_FLOAT, centmass, 
			 NULL);

        		
	free(body);
    free(above);
    free(below);
    SDFclose(sdfp);
	
	return 0;
}