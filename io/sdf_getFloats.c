#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
//#include "singlio.h"
#include "Msgs.h"
#include <stdlib.h>
#include "SPHbody.h"

int usage()
{
	printf("\t Gets thermovars for a chosen particle id\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_getfloats [sdf file] [id1] {id2} {id3} ...\n");
	return 0;
}

int main(int argc, char **argv[])
{
	double dist_in_cm;
	double mass_in_g;
	double dens_in_gccm;
	double energy_in_erg;
	double time_in_s;
	double specenergy_in_ergperg;
	double pressure_in_ergperccm;
	
	
	time_in_s				= 1;
	dist_in_cm				= 6.955e7;
	mass_in_g				= 1.989e27;
	dens_in_gccm			= mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	
	int id;
	double D;
	
	int i,j;

	if (argc < 2){
		usage();
		return 0;
	}

	
	
	SDF *sdfp;
	SPHbody *inArray;
	int gnobj, nobj;
	int conf;
	sdfp = SDFreadf(argv[1], (void **)&inArray, &gnobj, &nobj, sizeof(SPHbody),
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
					"vsound", offsetof(SPHbody, vsound), &conf,
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
					"sigma", offsetof(SPHbody, sigma), &conf,
					"kappa", offsetof(SPHbody, kappa), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	//now to fill the misc. floats and ints//


	
	
	//Msgf(("Data read, SPHnobj=%d, SPHgnobj=%d\n", nobj, gnobj));
	//SDFout *outArray = malloc(sizeof(SDFout) * nobj);

	printf("id x y z mass h u rho temp sigma D alpha abar zbar kappa\n");
	
	for (j=2; j<argc;j++){
		id = atoi(argv[j]);
		for (i=0; i<nobj; i+=1) {
			if (inArray[i].ident == id){
				printf("%06d ",id);
				singlPrintf("%3.2e %3.2e %3.2e ", inArray[i].x,inArray[i].y,inArray[i].z);
				printf("%3.2e ",inArray[i].mass);
				printf("%3.2e ",inArray[i].h);
				singlPrintf("%3.2e %3.2e ", inArray[i].u, inArray[i].rho);
				printf("%3.2e ",inArray[i].temp);
				printf("%3.2e ",inArray[i].sigma);
				D = inArray[i].sigma * inArray[i].abar/(inArray[i].rho*6.02e23*1.435e-59);
				printf("%3.2e ",D);
				printf("%3.2e ",1/D);
				//			singlPrintf("drho_dt = %3.2e\n", inArray[i].drho_dt);
				//			singlPrintf("udot = %3.2e\n", inArray[i].udot);
				//			singlPrintf("temp = %3.2e\n", inArray[i].temp);
				singlPrintf("%3.2e ", inArray[i].abar);
				singlPrintf("%3.2e ", inArray[i].zbar);
				printf("%3.2e ",inArray[i].kappa);
				//			singlPrintf("ax = %3.2e\n", inArray[i].ax);
				//			singlPrintf("ay = %3.2e\n", inArray[i].ay);
				//			singlPrintf("az = %3.2e\n", inArray[i].az);
				//			singlPrintf("lax = %3.2e\n", inArray[i].lax);
				//			singlPrintf("lay = %3.2e\n", inArray[i].lay);
				//			singlPrintf("laz = %3.2e\n", inArray[i].laz);
				//			singlPrintf("phi = %3.2e\n", inArray[i].phi);
				//			singlPrintf("tacc = %3.2e\n", inArray[i].tacc);
				//			singlPrintf("idt = %3.2e\n", inArray[i].idt);
				//			singlPrintf("bdot = %3.2e\n\n\n", inArray[i].bdot);
			}
		}
		printf("\n");
	}


	return 0;
}