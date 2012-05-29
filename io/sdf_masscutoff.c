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


char outfile[80];
int swapped;

double dist_in_cm;
double mass_in_g;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;
double maxmass;
double mass,tmass,r,rmax;
int nsub;

double ke,pe,te;
float gam,tvel,Gnewt,eps,dt,tpos,tolerance,frac_tolerance;
int ndim,iter;

int usage()
{
	printf("\t Creates a subset file for a chosen mass cutoff. ((Must isothermalize the result!))\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_2grid [sdf file] [rho=1,temp=2] [pixels] [xmax(code units)] [rhomax(cgs)]\n");
	return 0;
}

void quicksort(double** data,int m,int n)
{
    int key,i,j,k;
    double *temp;
    if( m < n)
    {
        key = m;
        i = m;
        j = n;
        while(i < j)
        {
            while((i < n) && (data[i][0] <= data[key][0]))
                i++;
            while(data[j][0] > data[key][0])
                j--;
            if( i < j)
            {
                temp = data[i];
                data[i] = data[j];
                data[j] = temp;
            }
                
        }
        temp = data[key];
        data[key] = data[j];
        data[j] = temp;

        quicksort(data,m,j-1);
        quicksort(data,j+1,n);
    }
}

int main(int argc, char **argv[])
{
	int i,j,k;
    
    
	time_in_s				= 1;
	dist_in_cm				= 6.955e7;
	mass_in_g				= 1.989e27;
	dens_in_gccm			= mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	
	//maxrho					= 1e8;
	
    if (argc < 3){
		usage();
		return 0;
	}else {
        
	}
    
    double **ptcls;
    
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
    
    ptcls = malloc(nobj*sizeof(double *));
	for(i=0;i<nobj;i++){
		/* r, mass */
        ptcls[i]  = malloc(2*sizeof(double));
	}
    
    for(i=0;i<nobj;i++)
    {
        r = sqrt(pow(body[i].x,2.0)+pow(body[i].y,2.0)+pow(body[i].z,2.0));
        ptcls[i][0] = r;
        ptcls[i][1] = body[i].mass;        
    }
    
//    swapped = 1;
//    for(i=1;(i<=k) && swapped;i++)
//    {
//        swapped = 0;
//        for(j=0;j<(nobj-1);j++)
//        {
//            if(ptcls[j+1][0] < ptcls[j][0])
//            {
//                mass        = ptcls[j][1];
//                r           = ptcls[j][0];
//                ptcls[j][0] = ptcls[j+1][0];
//                ptcls[j][1] = ptcls[j+1][1];
//                ptcls[j+1][0] = r;
//                ptcls[j+1][1] = mass;
//                swapped = 1;
//            }
//        }
//        if(i % 2000 == 0) printf("%3.2f\%\n",(double)i/(double)nobj*100);
//    }
    
    quicksort(ptcls,0,nobj-1);
    
    nsub=tmass=0;
    maxmass = atof(argv[2])*(1.989e33);
    for(i=0;(i<nobj) && (tmass<=maxmass);i++)
    {
        tmass += ptcls[i][1]*mass_in_g;
        nsub++;
        rmax = ptcls[i][0];
    }
    
    printf("Making subset with %3.2fMsun. %d particles\n",tmass/(1.989e33),nsub);
    
    SPHbody *outbody = malloc(sizeof(SPHbody) * (nsub));
    
    j=0;
    for(i=0;(i<nobj) && (j<nsub);i++)
    {
        r = sqrt(pow(body[i].x,2.0)+pow(body[i].y,2.0)+pow(body[i].z,2.0));
        if(r<=rmax)
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
            //		outbody[j].gax = body[i].gax;
            //		outbody[j].gay = body[i].gay;
            //		outbody[j].gaz = body[i].gaz;
            //outbody[j].grav_mass = body[i].grav_mass;
            outbody[j].phi = body[i].phi;
            //outbody[j].tacc = body[i].tacc;
            outbody[j].idt = body[i].idt;
            outbody[j].nbrs = body[i].nbrs;
            outbody[j].ident = j;
            outbody[j].windid = body[i].windid;

            j++;
        }
    }
    
    snprintf(outfile, sizeof(outfile), "sub_%s", argv[1]);
    
    iter=tpos=0;
    SDFwrite(outfile, nsub, 
			 nsub, outbody, sizeof(SPHbody),
			 SPHOUTBODYDESC,
			 "npart", SDF_INT, nsub,
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
			 "ke", SDF_DOUBLE, ke,
			 "pe", SDF_DOUBLE, pe,
			 "te", SDF_DOUBLE, te,
			 NULL);
	free(outbody);    
		
	free(body);
    free(ptcls);
    SDFclose(sdfp);
	
	return 0;
}