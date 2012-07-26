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

int gnobj, nobj, ndim=3;
int conf;
float tpos,Gnewt,tvel,gma,tolerance,frac_tolerance,dt,eps;
double ke,pe,te;

char outfile[80];

double dist_in_cm;
double mass_in_g;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;
double mass_in_msol;

double kin,pot,mass;

int usage()
{
	printf("\t Determines how much mass is unbound from the system. \n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_masslost [sdf file]\n");
	return 0;
}

int main(int argc, char **argv[])
{
	if(argc < 2)
	{
		usage();
		return 0;
	}
	
	int i,j,k,id;    

    double cc = 2.99792e10;
	double cc2 = cc*cc;
    
	time_in_s				= 1;
	dist_in_cm				= 6.955e7;
	mass_in_g				= 1.989e27;
    mass_in_msol            = mass_in_g/1.99e33;
	dens_in_gccm			= mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
    
    double **ptcls;
    
    SDF *sdfp;
	SPHbody *body;
    Darkbody *nbody = malloc(sizeof(Darkbody)*1);;

		
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
	SDFgetfloatOrDefault(sdfp, "tvel",  &tvel, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "gamma",  &gma, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "ke",  &ke, 0.0);
	SDFgetdoubleOrDefault(sdfp, "pe",  &pe, 0.0);
	SDFgetdoubleOrDefault(sdfp, "te",  &te, 0.0);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
    
    nbody[0].ident = 0;

    j=0;
    for(i=0;i<nobj;i++)
    {
        kin = 0;
        pot = fabs(body[i].phi*body[i].mass);
		kin += pow(body[i].vx,2.0);
		kin += pow(body[i].vy,2.0);
		kin += pow(body[i].vz,2.0);
		kin = kin*0.5*body[i].mass;
        
		if(kin>pot) 
        {
            mass += body[i].mass;
            j++;
        } else 
        {
            nbody[0].mass += body[i].mass;
            nbody[0].x += body[i].mass*body[i].x;
            nbody[0].y += body[i].mass*body[i].y;
            nbody[0].z += body[i].mass*body[i].z;
            nbody[0].vx += body[i].mass*body[i].vx;
            nbody[0].vy += body[i].mass*body[i].vy;
            nbody[0].vz += body[i].mass*body[i].vz;
        }
    }
    

    nbody[0].x = nbody[0].x/nbody[0].mass;
    nbody[0].y = nbody[0].y/nbody[0].mass;
    nbody[0].z = nbody[0].z/nbody[0].mass;
    nbody[0].vx = nbody[0].vx/nbody[0].mass;
    nbody[0].vy = nbody[0].vy/nbody[0].mass;
    nbody[0].vz = nbody[0].vz/nbody[0].mass;
    
    nbody[0].lx = nbody[0].ly = nbody[0].lz = nbody[0].accmass= 0;
    

	printf("%d particles massing %3.3e g are unbound.\n",j,mass*mass_in_g);
    printf("Nbody particle %3.3e g:\ncm = [%3.3e %3.3e %3.3e]\n",nbody[0].mass*mass_in_g,nbody[0].x,nbody[0].y,nbody[0].z);
    printf("cv = [%3.3e %3.3e %3.3e]\n",nbody[0].vx,nbody[0].vy,nbody[0].vz);
    
    SPHbody *outbody = malloc(sizeof(SPHbody) * (j));
    
    
    j=0;
    for(i=0;i<nobj;i++)
    {
        kin = 0;
        pot = fabs(body[i].phi*body[i].mass);
		kin += pow(body[i].vx,2.0);
		kin += pow(body[i].vy,2.0);
		kin += pow(body[i].vz,2.0);
		kin = kin*0.5*body[i].mass;
        
		if(kin>pot) 
        {
            outbody[j].x = body[i].x;
            outbody[j].y = body[i].y;
            outbody[j].z = body[i].z;
            outbody[j].mass = body[i].mass;
            outbody[j].vx = body[i].vx;
            outbody[j].vy = body[i].vy;
            outbody[j].vz = body[i].vz;
            outbody[j].u = body[i].u;
            outbody[j].h = body[i].h;
            outbody[j].rho = body[i].rho;
            outbody[j].pr = body[i].pr;
            outbody[j].drho_dt = body[i].drho_dt;
            outbody[j].udot = body[i].udot;
            outbody[j].temp = body[i].temp;
            outbody[j].He4 = body[i].He4;
            outbody[j].C12 = body[i].C12;
            outbody[j].O16 = body[i].O16;
            outbody[j].Ne20 = body[i].Ne20;
            outbody[j].Mg24 = body[i].Mg24;
            outbody[j].Si28 = body[i].Si28;
            outbody[j].S32 = body[i].S32;
            outbody[j].Ar36 = body[i].Ar36;
            outbody[j].Ca40 = body[i].Ca40;
            outbody[j].Ti44 = body[i].Ti44;
            outbody[j].Cr48 = body[i].Cr48;
            outbody[j].Fe52 = body[i].Fe52;
            outbody[j].Ni56 = body[i].Ni56;
            outbody[j].vsound = body[i].vsound;
            outbody[j].abar = body[i].abar;
            outbody[j].zbar = body[i].zbar;
            outbody[j].ax = body[i].ax;
            outbody[j].ay = body[i].ay;
            outbody[j].az = body[i].az;
            outbody[j].lax = body[i].lax;
            outbody[j].lay = body[i].lay;
            outbody[j].laz = body[i].laz;
            outbody[j].phi = body[i].phi;
            outbody[j].idt = body[i].idt;
            outbody[j].nbrs = body[i].nbrs;
            outbody[j].ident = i;
            outbody[j].windid = body[i].windid;
            j++;
        }
    }
    
    snprintf(outfile, sizeof(outfile), "sph_%s", argv[1]);
    
    SDFwrite(outfile, j, 
			 j, outbody, sizeof(SPHbody),
			 SPHOUTBODYDESC,
			 "npart", SDF_INT, j,
			 "iter", SDF_INT, 0,
			 "dt", SDF_FLOAT, dt,
			 "eps", SDF_FLOAT, eps,
			 "Gnewt", SDF_FLOAT, Gnewt,
			 "tolerance", SDF_FLOAT, tolerance,
			 "frac_tolerance", SDF_FLOAT, frac_tolerance,
			 "ndim", SDF_INT, ndim,
			 "tpos", SDF_FLOAT, tpos,
			 "tvel", SDF_FLOAT, tvel,
			 "gamma", SDF_FLOAT, gma,
			 "ke", SDF_DOUBLE, ke,
			 "pe", SDF_DOUBLE, pe,
			 "te", SDF_DOUBLE, te,
			 NULL);
    
    snprintf(outfile, sizeof(outfile), "nbody_%s", argv[1]);
    
    pe = 0;
    ke = 0.5*nbody[0].mass*(pow(nbody[0].vx,2.0)+pow(nbody[0].vy,2.0)+pow(nbody[0].vz,2.0));
    te = ke;
    
    SDFwrite(outfile, 1, 
			 1, nbody, sizeof(Darkbody),
			 DARKOUTBODYDESC,
			 "npart", SDF_INT, 1,
			 "iter", SDF_INT, 0,
			 "dt", SDF_FLOAT, dt,
			 "eps", SDF_FLOAT, eps,
			 "Gnewt", SDF_FLOAT, Gnewt,
			 "tolerance", SDF_FLOAT, tolerance,
			 "frac_tolerance", SDF_FLOAT, frac_tolerance,
			 "ndim", SDF_INT, ndim,
			 "tpos", SDF_FLOAT, tpos,
			 "tvel", SDF_FLOAT, tvel,
			 "gamma", SDF_FLOAT, gma,
			 "ke", SDF_DOUBLE, ke,
			 "pe", SDF_DOUBLE, pe,
			 "te", SDF_DOUBLE, te,
			 NULL);
    
    free(outbody);
    free(nbody);
	free(body);
    SDFclose(sdfp);
	
	return 0;
}