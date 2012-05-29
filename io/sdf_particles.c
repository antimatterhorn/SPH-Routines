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
float x,rho,temp,mass,pressure,h,cs,vx,tpos,u;
int i,j=0;
double g = 4.0/3.0;
double mass_in_g = 1.989E+27;        
double dist_in_cm = 6.955e7;     
double time_in_s = 1;
double energy_in_erg,dens_in_gccm,pressure_in_ergperccm,specenergy_in_ergperg;

char sdffile[80];
char asciifile[80];

int main(int argc, char **argv[])
{
	energy_in_erg = mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
	dens_in_gccm = mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
	pressure_in_ergperccm = energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	specenergy_in_ergperg = energy_in_erg/mass_in_g;
	
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
					//"grav_mass", offsetof(SPHbody, grav_mass), &conf,
					"phi", offsetof(SPHbody, phi), &conf,
					"tacc", offsetof(SPHbody, tacc), &conf,
					"idt", offsetof(SPHbody, idt), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	
	
	snprintf(asciifile, sizeof(asciifile), "%s.ascii", argv[1]);
	FILE *stream, *fopen();
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(asciifile,"w");
	
	fprintf(stream,"**** header **************\n");
	fprintf(stream,"N-part\t= %d\n",nobj);
	fprintf(stream,"t-pos\t= %f\n",tpos);
	fprintf(stream,"**************************\n");
	fprintf(stream,"\n");
	fprintf(stream,"x y z mass vx vy vz u h rho pr temp cs He4 C12 O16 Ne20 Mg24 Si28 S32 Ar36 Ca40 Ti44 Cr48 Fe52 Ni56\n");
	
	for(i = 0; i < nobj; i++)
	{
		fprintf(stream,"%e ",body[i].x * dist_in_cm);
		fprintf(stream,"%e ",body[i].y * dist_in_cm);
		fprintf(stream,"%e ",body[i].z * dist_in_cm);
		fprintf(stream,"%e ",body[i].mass * mass_in_g);
		fprintf(stream,"%e ",body[i].vx * dist_in_cm);
		fprintf(stream,"%e ",body[i].vy * dist_in_cm);
		fprintf(stream,"%e ",body[i].vz * dist_in_cm);
		fprintf(stream,"%e ",body[i].u * specenergy_in_ergperg);
		fprintf(stream,"%e ",body[i].h * dist_in_cm);
		fprintf(stream,"%e ",body[i].rho * dens_in_gccm);
		fprintf(stream,"%e ",body[i].pr * pressure_in_ergperccm);
		fprintf(stream,"%e ",body[i].temp);
		fprintf(stream,"%e ",body[i].vsound * dist_in_cm);
		fprintf(stream,"%e ",body[i].He4);
		fprintf(stream,"%e ",body[i].C12);
		fprintf(stream,"%e ",body[i].O16);
		fprintf(stream,"%e ",body[i].Ne20);
		fprintf(stream,"%e ",body[i].Mg24);
		fprintf(stream,"%e ",body[i].Si28);
		fprintf(stream,"%e ",body[i].S32);
		fprintf(stream,"%e ",body[i].Ar36);
		fprintf(stream,"%e ",body[i].Ca40);
		fprintf(stream,"%e ",body[i].Ti44);
		fprintf(stream,"%e ",body[i].Cr48);
		fprintf(stream,"%e ",body[i].Fe52);
		fprintf(stream,"%e\n",body[i].Ni56);
		
	}
	fclose(stream);
	

	return 0;
}