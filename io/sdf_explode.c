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
double mass,tmass,r,rr,rmax,rmin,e1,e2,e3,x,y,z,xp,yp,zp;
int nexp;

double ke,pe,te;
float gam,tvel,Gnewt,eps,dt,tpos,tolerance,frac_tolerance;
int ndim,iter;

double dettemp,cfact,f;

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
	int i,j,k,id;
    double v = 8e8;
    
    nexp = 30;
    dettemp = 1e9;
    
	time_in_s				= 1;
	dist_in_cm				= 6.955e7;
	mass_in_g				= 1.989e27;
	dens_in_gccm			= mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	
    v = v/dist_in_cm;
    
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
        ptcls[i]  = malloc(4*sizeof(double));
	}
    
    rmax = 0;
    
    for(i=0;i<nobj;i++)
    {
        r = sqrt(pow(body[i].x,2.0)+pow(body[i].y,2.0)+pow(body[i].z,2.0));
        ptcls[i][0] = r;
        ptcls[i][1] = body[i].ident;  
        ptcls[i][2] = body[i].h;
        ptcls[i][3] = 0;
        if(r>rmax) rmax=r;
    }
    printf("Found rmax = %3.2f\n",rmax);
    
    quicksort(ptcls,0,nobj-1);
    
    /*calc floor energy as efermi + 3/2*k*Tsmall and make sure udot 
     doesn't take u below that.
     note that this does not conserve energy since work is done 
     without removing it from the particle that did the work
     the hope is that this is so small that it hardly registers
     as an energy imbalance in the end - must check 
    efloor = (3/5)*2.50899e13*pow(p->rho*dens_in_gccm/p->abar,2/3)/p->abar;
    efloor += (3/2)*8.31447e7*small_temp/p->abar;*/
    
    cfact = (3/2)*8.31447e7/specenergy_in_ergperg;
    f = pow(2.0,-1.0/6.0);
	
    /* single thermal bomb method */
/*    
    for(i=0;i<nobj;i++)
    {
        id = 0;
        for(j=0;j<nexp;j++)
            if(ptcls[j][1]==body[i].ident) id=1;
        if(id==1)
        {
            printf("u0 = %3.2e\t",body[i].u);
            body[i].u += cfact*dettemp;
            printf("--> %3.2e\n",body[i].u);
            
            r = sqrt(pow(body[i].x,2.0)+pow(body[i].y,2.0)+pow(body[i].z,2.0));
        }
    }
 */

    /* atmosphere detonation method */

	
    rmin = rmax*0.1;
    rr = rmax*0.5;
    
    for(i=0;i<nobj;i++)
    {
        x = body[i].x;
        y = body[i].y;
        z = body[i].z;
        
        if(sqrt(pow(x-rr,2.0)+y*y+z*z)<rmin)
        {
            printf("u0 = %3.2e\t",body[i].u);
            body[i].u += cfact*dettemp;
            printf("--> %3.2e\n",body[i].u);
        }
    }

    
    /* multiple ignitions method */
    /*
    rmax = ptcls[nobj-1][0];
    rmin = 0.1*rmax;
    rr = 2.0*ptcls[0][2];

    srand((unsigned)time(NULL));
    
    for(j=0;j<10;j++)
    {
        e1 = ((double)rand()/(double)RAND_MAX)*rmin;        //r
        e2 = ((double)rand()/(double)RAND_MAX)*3.14159;     //theta
        e3 = ((double)rand()/(double)RAND_MAX)*2*3.14159;   //phi
        
        x = e1*cos(e3)*sin(e2);
        y = e1*sin(e3)*sin(e2);
        z = e1*cos(e2);
        
        for(i=0;i<nobj;i++)
        {
            xp = body[i].x;
            yp = body[i].y;
            zp = body[i].z;
            
            if((sqrt(pow(x-xp,2.0)+pow(y-yp,2.0)+pow(z-zp,2.0))<rr) && body[i].Ni56<1e-6)
            {
                body[i].u = body[i].u*2.0;
                body[i].rho = body[i].rho*2.0;
                body[i].h *= f;
                body[i].Ni56 = 1e-5;
            }
        }
    }
    

    */
    
    snprintf(outfile, sizeof(outfile), "det_%s", argv[1]);
    
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

        		
	free(body);
    free(ptcls);
    SDFclose(sdfp);
	
	return 0;
}