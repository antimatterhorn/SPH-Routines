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

int gnobj, nobj;
int conf;
int i=0;
float He4,C12,O16,Ne20,Mg24,Si28,S32;
float Ar36,Ca40,Ti44,Cr48,Fe52,Ni56;
float mtot;
char sdffile[80];

int main(int argc, char **argv[])
{
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
					"phi", offsetof(SPHbody, phi), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
	
	
	SPHbody *p;
	
	for(p = body; p < body+gnobj; p++)
	{
		//singlPrintf("%d / %d\n",i,gnobj-1);
		He4 += p->mass * p->He4;
		C12 += p->mass * p->C12;
		O16 += p->mass * p->O16;
		Ne20 += p->mass * p->Ne20;
		Mg24 += p->mass * p->Mg24;
		Si28 += p->mass * p->Si28;
		S32 += p->mass * p->S32;
		Ar36 += p->mass * p->Ar36;
		Ca40 += p->mass * p->Ca40;
		Ti44 += p->mass * p->Ti44;
		Cr48 += p->mass * p->Cr48;
		Fe52 += p->mass * p->Fe52;
		Ni56 += p->mass * p->Ni56;	
		mtot += p->mass;
		
		i++;
	}
	
	He4 = He4 * pow(10.,-6.);
	C12 = C12 * pow(10.,-6.);
	O16 = O16 * pow(10.,-6.);
	Ne20 = Ne20 * pow(10.,-6.);
	Mg24 = Mg24 * pow(10.,-6.);
	Si28 = Si28 * pow(10.,-6.);
	S32 = S32 * pow(10.,-6.);
	Ar36 = Ar36 * pow(10.,-6.);
	Ca40 = Ca40 * pow(10.,-6.);
	Ti44 = Ti44 * pow(10.,-6.);
	Cr48 = Cr48 * pow(10.,-6.);
	Fe52 = Fe52 * pow(10.,-6.);
	Ni56 = Ni56 * pow(10.,-6.);
	mtot = mtot * pow(10.,-6.);
	
	singlPrintf("\n\nHe4\t= %f\n",He4);
	singlPrintf("C12\t= %f\n",C12);
	singlPrintf("O16\t= %f\n",O16);
	singlPrintf("Ne20\t= %f\n",Ne20);
	singlPrintf("Mg24\t= %f\n",Mg24);
	singlPrintf("Si28\t= %f\n",Si28);
	singlPrintf("S32\t= %f\n",S32);
	singlPrintf("Ar36\t= %f\n",Ar36);
	singlPrintf("Ca40\t= %f\n",Ca40);
	singlPrintf("Ti44\t= %f\n",Ti44);
	singlPrintf("Cr48\t= %f\n",Cr48);
	singlPrintf("Fe52\t= %f\n",Fe52);
	singlPrintf("Ni56\t= %f\n",Ni56);
	singlPrintf("------------------\nTotal\t= %f\n\n",mtot);
	
	return 0;
}