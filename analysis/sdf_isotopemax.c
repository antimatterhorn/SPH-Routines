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
int conf;
int i=0;
float He4,C12,O16,Ne20,Mg24,Si28,S32;
float Ar36,Ca40,Ti44,Cr48,Fe52,Ni56;
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
	
	
	SPHbody *p;
	
	for(p = body; p < body+gnobj; p++)
	{
		if (p->C12 > 0.5) {	
			singlPrintf("He4	= %f\n",p->He4);
			singlPrintf("C12	= %f\n",p->C12);
			singlPrintf("O16	= %f\n",p->O16);
			singlPrintf("Ne20	= %f\n",p->Ne20);
			singlPrintf("Mg24	= %f\n",p->Mg24);
			singlPrintf("Si28	= %f\n",p->Si28);
			singlPrintf("S32	= %f\n",p->S32);
			singlPrintf("Ar36	= %f\n",p->Ar36);
			singlPrintf("Ca40	= %f\n",p->Ca40);
			singlPrintf("Ti44	= %f\n",p->Ti44);
			singlPrintf("Cr48	= %f\n",p->Cr48);
			singlPrintf("Fe52	= %f\n",p->Fe52);
			singlPrintf("Ni56	= %f\n",p->Ni56);
			singlPrintf("rho	= %f\n",p->rho);
			singlPrintf("temp	= %f\n",p->temp);
			singlPrintf("u	= %f\n",p->u);
			singlPrintf("h	= %f\n",p->h);
			singlPrintf("mass	= %f\n",p->mass);
			singlPrintf("xyz	= %f %f %f\n",p->x,p->y,p->z);
			singlPrintf("id	= %d\n",i);
			singlPrintf("-------------\n\n");
		}
		i++;
	}
	
	return 0;
}