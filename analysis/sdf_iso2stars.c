/*
 *  sdf2csv.c
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "SPHbody.h"
#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
#include "Msgs.h"
#include <stdlib.h>
#include <math.h>

int gnobj, nobj;
int conf,id;
int i,j;
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
		
	if (argc < 3){
		printf("minid: ");
		scanf(argv[2]);
	}
	
	id = atoi(argv[2]);
		
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
					"gax", offsetof(SPHbody, gax), &conf,
					"gay", offsetof(SPHbody, gay), &conf,
					"gaz", offsetof(SPHbody, gaz), &conf,
					"grav_mass", offsetof(SPHbody, grav_mass), &conf,
					"phi", offsetof(SPHbody, phi), &conf,
					"tacc", offsetof(SPHbody, tacc), &conf,
					"idt", offsetof(SPHbody, idt), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
	
	double **outArray;
	
	outArray = (double **) malloc(2*sizeof(double *));
	for(i=0;i<2;i++){
		outArray[i]  = (double *) malloc(14*sizeof(double));
	}
	
	for(i=0;i<2;i++){
		for(j=0;j<14;j++) outArray[i][j] = 0;
	}
	
	SPHbody *p;
	
	for(p = body; p < body+gnobj; p++)
	{
		if (p->ident < id) {
			outArray[0][0] += p->mass * p->He4;
			outArray[0][1] += p->mass * p->C12;
			outArray[0][2] += p->mass * p->O16;
			outArray[0][3] += p->mass * p->Ne20;
			outArray[0][4] += p->mass * p->Mg24;
			outArray[0][5] += p->mass * p->Si28;
			outArray[0][6] += p->mass * p->S32;
			outArray[0][7] += p->mass * p->Ar36;
			outArray[0][8] += p->mass * p->Ca40;
			outArray[0][9] += p->mass * p->Ti44;
			outArray[0][10] += p->mass * p->Cr48;
			outArray[0][11] += p->mass * p->Fe52;
			outArray[0][12] += p->mass * p->Ni56;	
			outArray[0][13] += p->mass;
		} else {
			outArray[1][0] += p->mass * p->He4;
			outArray[1][1] += p->mass * p->C12;
			outArray[1][2] += p->mass * p->O16;
			outArray[1][3] += p->mass * p->Ne20;
			outArray[1][4] += p->mass * p->Mg24;
			outArray[1][5] += p->mass * p->Si28;
			outArray[1][6] += p->mass * p->S32;
			outArray[1][7] += p->mass * p->Ar36;
			outArray[1][8] += p->mass * p->Ca40;
			outArray[1][9] += p->mass * p->Ti44;
			outArray[1][10] += p->mass * p->Cr48;
			outArray[1][11] += p->mass * p->Fe52;
			outArray[1][12] += p->mass * p->Ni56;	
			outArray[1][13] += p->mass;
		}
	}
	
	for(i=0;i<2;i++){
		for(j=0;j<14;j++) outArray[i][j] = outArray[i][j] * pow(10.,-6.);
	}

	singlPrintf("\nHe4\t= %f + %f = %f\n",outArray[0][0],outArray[1][0],outArray[0][0]+outArray[1][0]);
	singlPrintf("C12\t= %f + %f = %f\n",outArray[0][1],outArray[1][1],outArray[0][1]+outArray[1][1]);
	singlPrintf("O16\t= %f + %f = %f\n",outArray[0][2],outArray[1][2],outArray[0][2]+outArray[1][2]);
	singlPrintf("Ne20\t= %f + %f = %f\n",outArray[0][3],outArray[1][3],outArray[0][3]+outArray[1][3]);
	singlPrintf("Mg24\t= %f + %f = %f\n",outArray[0][4],outArray[1][4],outArray[0][4]+outArray[1][4]);
	singlPrintf("Si28\t= %f + %f = %f\n",outArray[0][5],outArray[1][5],outArray[0][5]+outArray[1][5]);
	singlPrintf("S32\t= %f + %f = %f\n",outArray[0][6],outArray[1][6],outArray[0][6]+outArray[1][6]);
	singlPrintf("Ar36\t= %f + %f = %f\n",outArray[0][7],outArray[1][7],outArray[0][7]+outArray[1][7]);
	singlPrintf("Ca40\t= %f + %f = %f\n",outArray[0][8],outArray[1][8],outArray[0][8]+outArray[1][8]);
	singlPrintf("Ti44\t= %f + %f = %f\n",outArray[0][9],outArray[1][9],outArray[0][9]+outArray[1][9]);
	singlPrintf("Cr48\t= %f + %f = %f\n",outArray[0][10],outArray[1][10],outArray[0][10]+outArray[1][10]);
	singlPrintf("Fe52\t= %f + %f = %f\n",outArray[0][11],outArray[1][11],outArray[0][11]+outArray[1][11]);
	singlPrintf("Ni56\t= %f + %f = %f\n",outArray[0][12],outArray[1][12],outArray[0][12]+outArray[1][12]);
	singlPrintf("----------------------------------------\nTotal\t= %f + %f = %f\n\n",
				outArray[0][13],outArray[1][13],outArray[0][13]+outArray[1][13]);
	
	free(outArray);
	
	return 0;
}