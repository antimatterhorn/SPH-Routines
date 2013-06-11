/*
 *  sdf_masscutoff.c
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  reads an sdf file and outputs a subset with a chosen mass
 *
 */

#include "../SPHbody.h"
#include "sdf_detonate.h"
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


char outfile[80];
int swapped;

double dist_in_cm;
double mass_in_g;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;
double mass_in_msol;
double mass,tmass,rho,energy,rhomax,rhomin,etot;
int nexp;
double mni,mcore;

double ke,pe,te;
float gam,tvel,Gnewt,eps,dt,tpos,tolerance,frac_tolerance;
int ndim,iter;

double dettemp,cfact,f;

double abundin[14],abundout[14],a1[14],a2[14],a3[14];

int usage()
{
	printf("\t Instantly burns a star to a predetermined composition and deposits that energy back to the particles. \n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_postexp [sdf file] {min rho (g/cc)}\n");
	return 0;
}

int burn(double rhoin)
{
	int i,j;
	double tot;
	
	if(rhoin > 5.0e9) rhoin = 5.0e9;
	if(rhoin < 1.0e4) rhoin = 1.0e4;
	
	if(mcore >= 1.2)
	{
		for(i=0;i<378;i++)
		{
			if(rho1p2[i] < rhoin)
			{
				for(j=0;j<14;j++)
				{
					a2[j] = abund1p2[i][j];
					a1[j] = abund1p2[i-1][j];
					abundout[j] = (a2[j]-a1[j])/(rho1p2[i]-rho1p2[i-1])*(rhoin-rho1p2[i-1])+a1[j];
				}
				break;
			}
		}
	}
	else if(mcore >= 1.0)
	{
		for(i=0;i<378;i++)
		{
			if(rho1p2[i] < rhoin)
			{
				for(j=0;j<14;j++)
				{
					a2[j] = abund1p2[i][j];
					a1[j] = abund1p2[i-1][j];
					a3[j] = (a2[j]-a1[j])/(rho1p2[i]-rho1p2[i-1])*(rhoin-rho1p2[i-1])+a1[j];
				}
				break;
			}
		}
		for(i=0;i<378;i++)
		{
			if(rho1p0[i] < rhoin)
			{
				for(j=0;j<14;j++)
				{
					a2[j] = abund1p0[i][j];
					a1[j] = abund1p0[i-1][j];
					abundout[j] = (a2[j]-a1[j])/(rho1p0[i]-rho1p0[i-1])*(rhoin-rho1p0[i-1])+a1[j];
					abundout[j] = (a3[j]-abundout[j])/(1.2-1.0)*(mcore-1.0)+abundout[j];
				}
				break;
			}
		}
		
	}
	else if(mcore > 0.8)
	{
		for(i=0;i<378;i++)
		{
			if(rho1p0[i] < rhoin)
			{
				for(j=0;j<14;j++)
				{
					a2[j] = abund1p0[i][j];
					a1[j] = abund1p0[i-1][j];
					a3[j] = (a2[j]-a1[j])/(rho1p0[i]-rho1p0[i-1])*(rhoin-rho1p0[i-1])+a1[j];
				}
				break;
			}
		}
		for(i=0;i<378;i++)
		{
			if(rho1p0[i] < rhoin)
			{
				for(j=0;j<14;j++)
				{
					a2[j] = abund0p8[i][j];
					a1[j] = abund0p8[i-1][j];
					abundout[j] = (a2[j]-a1[j])/(rho0p8[i]-rho0p8[i-1])*(rhoin-rho0p8[i-1])+a1[j];
					abundout[j] = (a3[j]-abundout[j])/(1.0-0.8)*(mcore-0.8)+abundout[j];
				}
				break;
			}
		}
	}
	else 
	{
		for(i=0;i<378;i++)
		{
			if(rho0p8[i] < rhoin)
			{
				for(j=0;j<14;j++)
				{
					a2[j] = abund0p8[i][j];
					a1[j] = abund0p8[i-1][j];
					abundout[j] = (a2[j]-a1[j])/(rho0p8[i]-rho0p8[i-1])*(rhoin-rho0p8[i-1])+a1[j];
				}
				break;
			}
		}
	}
	
	for(i=0;i<14;i++)
	{
		tot += abundout[i];
	}
		
	for(i=0;i<14;i++)
		abundout[i] = abundout[i]/tot;
	
    return 0;
}

int main(int argc, char **argv[])
{
	if(argc < 2)
	{
		usage();
		return 0;
	}
	
	if(argc < 3)
		rhomin = 5.0e4;
	else
		rhomin = atof(argv[2]);

	
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
	SDFgetdoubleOrDefault(sdfp, "ke",  &ke, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "pe",  &pe, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "te",  &te, (float)0.0);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);

    
    rhomax = 0;
    
    for(i=0;i<nobj;i++)
		if(body[i].rho*dens_in_gccm > rhomax) 
			rhomax = body[i].rho*dens_in_gccm;
    printf("Found rhomax = %3.2e g/cc\n",rhomax);
	
	if(rhomax >= rho1p2[1]) mcore = 1.2;
	else if(rhomax >= rho1p0[1]) mcore = (1.2-1.0)/(rho1p2[1]-rho1p0[1])*(rhomax-rho1p0[1])+1.0;
	else if(rhomax >= rho0p8[1]) mcore = (1.0-0.8)/(rho1p0[1]-rho0p8[1])*(rhomax-rho0p8[1])+0.8;
	else mcore = 0.8;
	printf("Found mcore is most similar to %f msol\n",mcore);

    
    for(i=0;i<nobj;i++)
    {
        energy = 0;
        
        rho = body[i].rho*dens_in_gccm;
        if(rho > rhomin)
        {
            burn(rho);
            
            for(j=0;j<14;j++) abundin[j] = 0.;
            
            abundin[1] = 0.50;
			abundin[2] = 0.50;
			
            for(j=0;j<14;j++)
				energy += (abundin[j]-abundout[j])*ea[j]/aa[j];
            energy = energy*6.02214e23/specenergy_in_ergperg;
			
            printf("u: %3.2e",body[i].u);
            body[i].u += energy;
            printf(" + %3.2e {",energy);
            etot += energy * body[i].mass;
            tmass += body[i].mass;
            for(j=0;j<14;j++) printf("%1.3f ",abundout[j]);
            printf("}\n");
			
			/*
			 He4 = 0
			 C12 = 1
			 O16 = 2
			 Ne20= 3
			 Mg24= 4
			 Si28= 5
			 S32 = 6
			 Ar36= 7
			 Ca40= 8
			 Ti44= 9
			 Cr48= 10
			 Fe52= 11
			 Fe54= 12
			 Ni56= 13
			 */
			
            body[i].He4  = abundout[0];
            body[i].C12  = abundout[1];
            body[i].O16  = abundout[2];
            body[i].Ne20 = abundout[3];
            body[i].Mg24 = abundout[4];
            body[i].Si28 = abundout[5];
            body[i].S32  = abundout[6];
            body[i].Ar36 = abundout[7];
            body[i].Ca40 = abundout[8];
            body[i].Ti44 = abundout[9]+abundout[10];	//Cr48 -> Ti48 ignoring V48
            body[i].Cr48 = abundout[11];				//Fe52 -> Cr48
            body[i].Fe52 = abundout[12];				//Fe54 is stable
            body[i].Ni56 = abundout[13];
			k++;
        }
        else
        {
//			body[i].He4	 = 0.;
//            body[i].C12  = 0.5;
//            body[i].O16  = 0.5;
//            body[i].Ne20 = 0.;
//            body[i].Mg24 = 0.;
//            body[i].Si28 = 0.;
//            body[i].S32  = 0.;
//            body[i].Ar36 = 0.;
//            body[i].Ca40 = 0.;
//            body[i].Ti44 = 0.;
//            body[i].Cr48 = 0.; 
//            body[i].Fe52 = 0.; 
//            body[i].Ni56 = 0.;
        }

    }

	printf("Burned %3.3e g\n",tmass*mass_in_g);
    printf("~ %3.3e erg added\n",etot*energy_in_erg);
    
    snprintf(outfile, sizeof(outfile), "exp_%s", argv[1]);
    
    iter=tpos=0;
    SDFwrite(outfile, nobj, 
			 nobj, body, sizeof(SPHbody),
			 SPHOUTBODYDESC,
			 "npart", SDF_INT, nobj,
			 "iter", SDF_INT, iter,
			 "dt", SDF_FLOAT, dt,
			 "eps", SDF_FLOAT, eps,
			 "Gnewt", SDF_FLOAT, Gnewt,
			 "tolerance", SDF_FLOAT, tolerance,
			 "frac_tolerance", SDF_FLOAT, frac_tolerance,
			 "ndim", SDF_INT, ndim,
			 "tpos", SDF_FLOAT, tpos,
			 "tvel", SDF_FLOAT, tvel,
			 "gamma", SDF_FLOAT, gam,
			 "ke", SDF_DOUBLE, ke,
			 "pe", SDF_DOUBLE, pe,
			 "te", SDF_DOUBLE, te,
			 NULL);

        		
	free(body);
    SDFclose(sdfp);
	
	return 0;
}