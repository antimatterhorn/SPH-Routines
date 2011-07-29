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

#define Fortran(x) x

#define Fortran2(x) x##_

void Fortran2(init_helm_table)();


void Fortran2(helmeos)(double *temperature,double *den_row,double *etot_row, double *abar_row, 
							   double *zbar_row);

void Fortran2(wrapper_helmeos)(int *npart, double *den_row, 
							   double *etot_row, double *abar_row, 
							   double *zbar_row, double *temperature,
							   double *pressure, double *entropy);

void Fortran2(azbar)(double *xmass, double *aion, double *zion, int *ionmax, 
					 double *ymass, double *abar, double *zbar);


int gnobj, nobj;
int conf;
double mass_in_g = 1.989E+27;        
double dist_in_cm = 6.955e7;     
double time_in_s = 1;
double energy_in_erg,dens_in_gccm,pressure_in_ergperccm,specenergy_in_ergperg;

double x,y,z,h,rho,temp,mass,pressure,cs,u,s,abar,zbar,dt;
int npart=1,i,j,k;
double small_temp=1e3;

char sdffile[80];
char csvfile[80];
char vszfile[80];


int main(int argc, char **argv[])
{
	energy_in_erg = mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
	dens_in_gccm = mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
	pressure_in_ergperccm = energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	specenergy_in_ergperg = energy_in_erg/mass_in_g;
	
	double abund[13], abund2[13], abar, zbar;
	double zarray[13]	={ 2,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28};
	double aarray[13]	={ 4, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56};
	double molarabund[13];
	int nion=13;
	int npart = 1;
		
	singlPrintf("Reading in Helm-table... \n"); 
	Fortran2(init_helm_table)();
	
	SDF *sdfp1;
	SPHbody *body1;
	
	if (argc < 2){
		printf("SDF file: ");
		gets (argv[1]);
	}
	
	sdfp1 = SDFreadf(argv[1], (void **)&body1, &gnobj, &nobj, sizeof(SPHbody),
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
	SDFgetdoubleOrDefault(sdfp1, "dt", &dt, (double)0.0);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
	
	k = nobj/1000;
	j=0;
	
	for(i = 0; i < nobj; i++){
		if(i%k==0) j++;
	}	
	
	double **outArray;
	
	outArray = (double **) malloc(j*sizeof(double *));
	for(i=0;i<j;i++){
		outArray[i]  = (double *) malloc(11*sizeof(double));
	}
	
	j=0;
	for(i = 0; i < nobj; i++){
		if(i%k==0){
			x			= body1[i].x * dist_in_cm;
			y			= body1[i].y * dist_in_cm;
			z			= body1[i].z * dist_in_cm;
			h			= body1[i].h * dist_in_cm;
			rho			= body1[i].rho * dens_in_gccm;
			temp		= body1[i].temp;
			mass		= body1[i].mass * mass_in_g;
			pressure	= body1[i].pr * pressure_in_ergperccm;
			cs			= body1[i].vsound * dist_in_cm;
			u			= body1[i].u * specenergy_in_ergperg;
			abar		= body1[i].abar;
			zbar		= body1[i].zbar;
			s			= 0;
			
			abund[0] = body1[i].He4;
			abund[1] = body1[i].C12;
			abund[2] = body1[i].O16;
			abund[3] = body1[i].Ne20;
			abund[4] = body1[i].Mg24;
			abund[5] = body1[i].Si28;
			abund[6] = body1[i].S32;
			abund[7] = body1[i].Ar36;
			abund[8] = body1[i].Ca40;
			abund[9] = body1[i].Ti44;
			abund[10] = body1[i].Cr48;
			abund[11] = body1[i].Fe52;
			abund[12] = body1[i].Ni56;	
			
			Fortran2(wrapper_helmeos)(&npart, &rho, &u, &abar, &zbar, &temp, &pressure, &s);
			
			outArray[j][0] = x;
			outArray[j][1] = y;
			outArray[j][2] = z;
			outArray[j][3] = h;
			outArray[j][4] = rho;
			outArray[j][5] = temp;
			outArray[j][6] = mass;
			outArray[j][7] = pressure;
			outArray[j][8] = cs;
			outArray[j][9] = u;
			outArray[j][10] = s;
			j++;
		}				
	}
	

	
	
	
	
	
	
	
	
	
	printf("printing files...\n");
	
	//open the stream file
	snprintf(csvfile, sizeof(csvfile), "%s_ent.csv", argv[1]);
	FILE *stream, *fopen(), *vsz;
		
	stream = fopen(csvfile,"w");
	
	fprintf(stream,"x,y,z,h,rho,temp,mass,pressure,cs,u,s\n");
	
	for(i = 0; i < j; i++)
	{
		fprintf(stream,"%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e\n", 
				outArray[i][0],outArray[i][1],outArray[i][2],
				outArray[i][3],outArray[i][4],outArray[i][5],
				outArray[i][6],outArray[i][7],outArray[i][8],
				outArray[i][9],outArray[i][10]);
	}
	fclose(stream);
/*	
	char *cwd = getcwd(NULL, 0);
	snprintf(vszfile, sizeof(vszfile), "%s_edot.vsz", argv[1]);
	
	vsz = fopen(vszfile,"w");
	
	fprintf(vsz,"ImportFileCSV(u'%s/%s', linked=True)\n",cwd,csvfile);
	fprintf(vsz,"Add('page', name='page1', autoadd=False)\n");
	fprintf(vsz,"To('page1')\n");
	fprintf(vsz,"Add('graph', name='graph1', autoadd=False)\n");
	fprintf(vsz,"To('graph1')\n");
	fprintf(vsz,"Set('leftMargin', '1.7cm')\n");
	fprintf(vsz,"Set('rightMargin', '1.61cm')\n");
	fprintf(vsz,"Set('topMargin', '0.276cm')\n");
	fprintf(vsz,"Set('bottomMargin', '1.71cm')\n");
	fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
	fprintf(vsz,"To('x')\n");
	fprintf(vsz,"Set('label', u'x [cm]')\n");
	fprintf(vsz,"Set('min', -150000000.0)\n");
	fprintf(vsz,"Set('max', 150000000.0)\n");
	fprintf(vsz,"Set('autoExtend', False)\n");
	fprintf(vsz,"Set('lowerPosition', 1.0)\n");
	fprintf(vsz,"Set('upperPosition', 0.0)\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
	fprintf(vsz,"To('y')\n");
	fprintf(vsz,"Set('min', 0.0)\n");
	fprintf(vsz,"Set('max', 4.0)\n");
	fprintf(vsz,"Set('log', False)\n");
	fprintf(vsz,"Set('direction', 'vertical')\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('Label/color', u'red')\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"Set('TickLabels/color', u'red')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy1', autoadd=False)\n");
	fprintf(vsz,"To('xy1')\n");
	fprintf(vsz,"Set('yData', u'ratio2')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('PlotLine/width', u'1pt')\n");
	fprintf(vsz,"Set('PlotLine/hide', False)\n");
	fprintf(vsz,"Set('PlotLine/color', u'red')\n");
	fprintf(vsz,"Set('MarkerLine/color', u'red')\n");
	fprintf(vsz,"Set('MarkerFill/color', u'red')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('function', name='function1', autoadd=False)\n");
	fprintf(vsz,"To('function1')\n");
	fprintf(vsz,"Set('function', u'1')\n");
	fprintf(vsz,"Set('variable', u'x')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	
	
	fclose(vsz);
*/	
	free(body1);
	
	return 0;
}