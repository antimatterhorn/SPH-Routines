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

double mass,tmass,rho,energy,rhomax,rhomin,etot;
int nexp;
double mni,mcore;

double ke,pe,te;
float gam,tvel,Gnewt,eps,dt,tpos,tolerance,frac_tolerance;
int ndim,iter;

double dettemp,cfact,f;

double abundin[14],abundout[14],a1[14],a2[14],a3[14];
double output[2][1000];

int usage()
{
	printf("\t Computes the abundance vector and energy for a range of preburn densities. \n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_\n");
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
	int i,j,k,id;    

    double cc = 2.99792e10;
	double cc2 = cc*cc;
    
    rhomax = 4e9;
	
	if(rhomax >= rho1p2[1]) mcore = 1.2;
	else if(rhomax >= rho1p0[1]) mcore = (1.2-1.0)/(rho1p2[1]-rho1p0[1])*(rhomax-rho1p0[1])+1.0;
	else if(rhomax >= rho0p8[1]) mcore = (1.0-0.8)/(rho1p0[1]-rho0p8[1])*(rhomax-rho0p8[1])+0.8;
	else mcore = 0.8;
	//printf("Found mcore is most similar to %f msol\n",mcore);
	
	
	
	
	printf("rho,e,He,C,O,Ne,Mg,Si,S,Ar,Ca,Ti,Cr,Fe,Ni\n");
    
    for(i=0;i<1000;i++)
    {
        energy = 0;
        
        rho = 1.0e3 + pow(10.0,9.6/1000*i);
		burn(rho);
		
		for(j=0;j<14;j++) abundin[j] = 0.;
		
		abundin[1] = 0.50;
		abundin[2] = 0.50;
		
		for(j=0;j<14;j++)
			energy += (abundin[j]-abundout[j])*ea[j]/aa[j];
		energy = energy*6.02214e23;

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
		
		printf("%3.4e,%3.4e,",rho,energy);
		for(j=0;j<9;j++)
			printf("%3.4e,",abundout[j]*10.0);
		printf("%3.4e,",10.0*(abundout[9]+abundout[10]));	//Cr48 -> Ti48 ignoring V48
		printf("%3.4e,",abundout[11]*10.0);				//Fe52 -> Cr48
		printf("%3.4e,",abundout[12]*10.0);				//Fe54 is stable
		printf("%3.4e\n",abundout[13]*10.0);
    }
	
	return 0;
}