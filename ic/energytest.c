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
double mass_in_msol;
double maxmass;
double mass,tmass,r,rr,rmax,rmin,e1,e2,e3,x,y,z,xp,yp,zp,m,b,rho;
int nexp;
double mni,mcore;

double ke,pe,te;
float gam,tvel,Gnewt,eps,dt,tpos,tolerance,frac_tolerance;
int ndim,iter;

double dettemp,cfact,f;

double abundin[13],abundout[13];

double za[13] = {2.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28.};
double aa[13] = {4.,12.,16.,20.,24.,28.,32.,36.,40.,44.,48.,52.,56.};
double ma[13] = {4.00260325415,12.,15.99491461956,20.,24.,27.9769265325,32.,36.,40.,44.,48.,52.,55.942132};

double igex[14] = {0.,3.642e6,4.098e6,4.963e6,6.1e6,6.763e6,7.499e6,8.313e6,9.354e6,1.150e7,
                    1.294e7,1.477e7,1.898e7,1.e8};
double igey[14] = {0.,0.003,0.014,0.033,0.052,0.071,0.108,0.184,0.373,0.691,0.831,0.917,1.0,1.0};

double imex[30] = {0.,4.9e4,6.4e4,7.89e4,1.49e6,2.314e5,2.642e5,2.973e5,3.601e5,5.943e5,1.279e6,
                    1.549e6,1.795e6,2.020e6,2.306e6,2.595e6,2.877e6,3.285e6,3.696e6,4.890e6,6.191e6,
                    6.763e6,7.499e6,8.192e6,9.082e6,1.133e7,1.294e7,1.499e7,1.955e7,1.e8};
double imey[30] = {0.,0.,0.010,0.040,0.260,0.388,0.407,0.415,0.422,0.430,0.441,0.452,0.475,0.532,
                    0.660,0.834,0.944,0.993,0.997,0.970,0.955,0.933,0.895,0.815,0.706,0.331,0.169,
                    0.089,0.,0.};

double carbx[12] = {0.,3.95e4,5.2995e4,6.51e4,8.13e4,1.21e5,1.444e5,1.698e5,2.247e5,2.681e5,3.933e5,1.e8};
double carby[12] = {0.5,0.498,0.479,0.449,0.388,0.180,0.101,0.063,0.021,0.010,0.0,0.};

double oxygx[20] = {0.,5.222e4,6.419e4,1.210e5,1.423e5,1.698e5,2.214e5,2.603e5,3.107e5,5.205e5,
                    1.260e6,1.397e6,1.549e6,1.743e6,2.050e6,2.306e6,2.963e6,3.285e6,3.696e6,1.e8};
double oxygy[20] = {0.5,0.517,0.539,0.634,0.641,0.634,0.596,0.581,0.577,0.573,0.566,0.562,0.547,0.528,
                    0.468,0.350,0.059,0.010,0.,0.};

int usage()
{
	printf("\t Instantly burns a star to a predetermined composition and deposits that energy back to the particles. \n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_postexp [sdf file]\n");
	return 0;
}

int burn(double rhoin)
{
    int i;
    double atot=0;
    
    for(i=0;i<13;i++)
    {
        abundout[i] = 0.;
    }
    
    abundout[1] = abundout[2] = 0.5;
    
    //carbon
    if(rhoin > 2.e4)
    {
        for(i=0;i<12;i++)
        {
            if(carbx[i]>rhoin) break;
        }
        abundout[1] = (carby[i]-carby[i-1])/(carbx[i]-carbx[i-1])*(rhoin-carbx[i-1])+carby[i-1];
    }
    
    //oxygen
    if(rhoin > 2.e4)
    {
        for(i=0;i<20;i++)
        {
            if(oxygx[i]>rhoin) break;
        }
        abundout[2] = (oxygy[i]-oxygy[i-1])/(oxygx[i]-oxygx[i-1])*(rhoin-oxygx[i-1])+oxygy[i-1];
    }
    
    //silicon
    if(rhoin > 4.9e4)
    {
        for(i=0;i<30;i++)
        {
            if(imex[i]>rhoin) break;
        }
        abundout[5] = (imey[i]-imey[i-1])/(imex[i]-imex[i-1])*(rhoin-imex[i-1])+imey[i-1];
    }
    
    //nickel
    if(rhoin > 3.6e6)
    {
        for(i=0;i<14;i++)
        {
            if(igex[i]>rhoin) break;
        }
        abundout[12] = (igey[i]-igey[i-1])/(igex[i]-igex[i-1])*(rhoin-igex[i-1])+igey[i-1];
    }
    
    for(i=0;i<13;i++)
    {
        if(abundout[i]<0) abundout[i] = 0;
        atot += abundout[i];
    }
        
    if(atot>0) for(i=0;i<13;i++) abundout[i] = abundout[i]/atot;
    
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
	if(argc < 2)
	{
		usage();
		return 0;
	}
	
	int i,j,k,id;    

    double cc = 2.99792e10;
    
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
        ptcls[i]  = malloc(4*sizeof(double));
	}
    
    rmax = 0;
    
    for(i=0;i<nobj;i++)
    {
        r = sqrt(pow(body[i].x,2.0)+pow(body[i].y,2.0)+pow(body[i].z,2.0));
        ptcls[i][0] = r;
        ptcls[i][1] = body[i].ident;  
        ptcls[i][2] = body[i].mass;
        ptcls[i][3] = 0;
        if(r>rmax) rmax=r;
    }
    printf("Found rmax = %3.2f\n",rmax);
    
    quicksort(ptcls,0,nobj-1);

    /*
    mass = 0;
    for(i=0;i<nobj;i++)
    {
        mass += ptcls[i][2];
        r = ptcls[i][0];
        if(mass*mass_in_msol<mni) rr=r;
		if(mass*mass_in_msol>mcore) break;
    }
    */
    
    for(i=0;i<nobj;i++)
    {
        x = body[i].x;
        y = body[i].y;
        z = body[i].z;
        mass = 0;
        
        rho = body[i].rho*dens_in_gccm;
        if(rho > 1.e6)
        {
            burn(rho);
            
            for(j=0;j<13;j++) abundin[j] = 0.;
            /*
             abundin[0] = body[i].He4;
             abundin[1] = body[i].C12;
             abundin[2] = body[i].O16;
             abundin[5] = body[i].Si28;
             */
            
            abundin[1] = abundin[2] = 0.5;
            for(j=0;j<12;j++)
            {
                for(k=j;k<13;k++)
                {
                    mass += abundin[j]*abundout[k]*(1-(aa[j]/aa[k])*(ma[k]/ma[j]));
                }
            }
            mass = mass*cc*cc/specenergy_in_ergperg;
            printf("u: %3.2e",body[i].u);
            body[i].u += mass;
            printf("->%3.2e\n",body[i].u);
            body[i].He4 = abundout[0];
            body[i].C12 = abundout[1];
            body[i].O16 = abundout[2];
            body[i].Si28 = abundout[5];
            body[i].Ni56 = abundout[12];
        }

 
/*
        if(sqrt(x*x+y*y+z*z)<rr) //turn to Ni56
        {
           
            
            mass += (56./4.*mHe4 - mNi56)*body[i].He4/mHe4;
            mass += (56./12.*mC12 - mNi56)*body[i].C12/mC12;
            mass += (56./16.*mO16 - mNi56)*body[i].O16/mO16;
            
            body[i].He4     = 0.0;
            body[i].C12     = 0.0;
            body[i].O16     = 0.0;
            body[i].Ne20    = 0.0;
            body[i].Mg24    = 0.0;
            body[i].Si28    = 0.0;
            body[i].S32     = 0.0;
            body[i].Ar36    = 0.0;
            body[i].Ca40    = 0.0;
            body[i].Ti44    = 0.0;
            body[i].Cr48    = 0.0;
            body[i].Fe52    = 0.0;
            body[i].Ni56    = 1.0;
        }
        else if(sqrt(x*x+y*y+z*z)<r) //turn to Si28
        {
            mass += (28./4.*mHe4 - mSi28)*body[i].He4/mHe4;
            mass += (28./12.*mC12 - mSi28)*body[i].C12/mC12;
            mass += (28./16.*mO16 - mSi28)*body[i].O16/mO16;
            
            body[i].He4 = 0.0;
            body[i].C12 = 0.0;
            body[i].O16 = 0.0;
            body[i].Ne20 = 0.0;
            body[i].Mg24 = 0.0;
            body[i].Si28 = 1.0;
            body[i].S32 = 0.0;
            body[i].Ar36 = 0.0;
            body[i].Ca40 = 0.0;
            body[i].Ti44 = 0.0;
            body[i].Cr48 = 0.0;
            body[i].Fe52 = 0.0;
            body[i].Ni56 = 0.0;
        }
*/
        

    }

    
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