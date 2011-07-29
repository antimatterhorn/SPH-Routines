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

void Fortran2(burner)(double *dt, double *temperature, double *density, 
					  double *energy, double *bund, double *temperatureout, double *energyout, double *abundout); 

void Fortran2(init_burner)();  

void Fortran2(helmeos)(double *temperature,double *den_row,double *etot_row, double *abar_row, 
							   double *zbar_row);

void Fortran2(wrapper_helmeos)(int *npart, double *den_row, 
							   double *etot_row, double *abar_row, 
							   double *zbar_row, double *temperature,
							   double *pressure);

void Fortran2(azbar)(double *xmass, double *aion, double *zion, int *ionmax, 
					 double *ymass, double *abar, double *zbar);

int gnobj, nobj;
int conf;
double x,rho,temp,mass,pressure,h,cs,u,eout,edot,dt,tempout;
double v,rbar,pri,prj,rhoi,rhoj,dw;
int i,j=0,k,na;
int nbrs=5;
int np;
double mass_in_g = 1.989E+27;        
double dist_in_cm = 6.955e7;     
double time_in_s = 1;
double energy_in_erg,dens_in_gccm,pressure_in_ergperccm,specenergy_in_ergperg;


char sdffile[80];
char csvfile[80];
char vszfile[80];

double delW(double r, double h)
{
	v = r/h;
	if (v>0 && v<=1)
	{
		return (1/h) * (2.25*v*v-3*v) * pow(h,-3.0) * pow(3.14159,-1.0);
	} else if (v>1 && v<=2)
	{
		return (1/h) * (-0.75*pow(2-v,2.0)) * pow(h,-3.0) * pow(3.14159,-1.0);
	} else {
		return 0;
	}

	return 0;
}

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
	singlPrintf("Starting Aprox13 Burner... \n"); 
	Fortran2(init_burner)();
	
	SDF *sdfp1;
	SPHbody *body1;
	
	if (argc < 2){
		printf("SDF file: ");
		gets (argv[1]);
	}
	
	if (argc > 2) nbrs = atoi(argv[2]);
	
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
	
	for(i = 0; i < nobj; i++){
		if ((fabs(body1[i].y) < (body1[i].h*2.0)) && (fabs(body1[i].z) < (body1[i].h*2.0))){
			j++;
		}
	}
	
	np = j;
	
	singlPrintf("Using %d points.\n",j);
	
	na = 13;
	
	double **outArray;
	double tempArray[na];
	
	outArray = (double **) malloc(nobj*sizeof(double *));
	for(i=0;i<j;i++){
		outArray[i]  = (double *) malloc(na*sizeof(double));
	}
	
	//0,  1,   2,   3, 4, 5,6,7,   8,   9,  10,   11,   12
	//x,rho,temp,mass,pr,cs,u,h,edot,udot,crit,cond1,cond2
	
	printf("filling array...\n");

	j=0;
	
	for(i = 0; i < nobj; i++){
		if ((fabs(body1[i].y) < (body1[i].h*2.0)) && (fabs(body1[i].z) < (body1[i].h*2.0))){
			x			= body1[i].x * dist_in_cm;
			rho			= body1[i].rho * dens_in_gccm;
			temp		= body1[i].temp;
			mass		= body1[i].mass * mass_in_g;
			pressure	= body1[i].pr * pressure_in_ergperccm;
			cs			= body1[i].vsound * dist_in_cm;
			u			= body1[i].u * specenergy_in_ergperg;
			eout		= 0;
			
			outArray[j][0] = x;
			outArray[j][1] = rho;
			outArray[j][2] = temp;
			outArray[j][3] = mass;
			outArray[j][4] = pressure;
			outArray[j][5] = cs;
			outArray[j][6] = u;
			outArray[j][7] = body1[i].h * dist_in_cm;
			for(k=8;k<na;k++) outArray[j][k] = 0;
			//outArray[j][9] = body1[i].udot * specenergy_in_ergperg;
			
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
			
			for(k=0;k<13;k++)
			{
				abund2[k] = 0;
			}
			
			if (temp > 1e6) //can change this to 1e8 later
			{
				//printf("\n\n-------before helm---------\n");
				//printf("rho\t=%e\ntemp\t=%e\npres\t=%e\nu\t=%e\nabar\t=%e\nzbar\t=%e\n",rho,temp,pressure,u,abar,zbar);
				
				Fortran2(azbar)(abund, aarray, zarray, &nion, molarabund, &abar, &zbar);
				//Fortran2(wrapper_helmeos)(&npart, &rho, &u, &abar, &zbar, &temp, &pressure); 
				//printf("-------after helm---------\n");
				//printf("rho\t=%e\ntemp\t=%e\npres\t=%e\nu\t=%e\nabar\t=%e\nzbar\t=%e\n",rho,temp,pressure,u,abar,zbar);
				
										 
										 
				Fortran2(burner)(&dt, &temp,  
								 &rho, &u, abund, &tempout, &eout, 
								 abund2);
				//printf("-------after burn---------\n");
				//printf("rho\t=%e\ntemp\t=%e\npres\t=%e\nu\t=%e\nabar\t=%e\nzbar\t=%e\neout\t=%e\n",rho,temp,pressure,u,abar,zbar,eout);
				outArray[j][8] = (eout - u)/dt;
			}
			
			j++;
		}
	}
	
	printf("sorting...\n");
	
	//sort
	for (i=0; i< (np -1); i++)    // element to be compared
    {
		for(j = (i+1); j < np; j++)   // rest of the elements
		{
			if (outArray[i][0] < outArray[j][0])          // descending order
			{
				for (k=0;k<na;k++) 
				{
					tempArray[k] = outArray[i][k];
					outArray[i][k] = outArray[j][k];
					outArray[j][k] = tempArray[k];
				}
			}
		}
	}
	
	for(j=nbrs; j<(np-nbrs); j++)
	{
		rhoi = outArray[j][1];
		pri = outArray[j][4];
		cs = outArray[j][5];
		
		for(k=-nbrs;k<(nbrs+1);k++)
		{
			rbar = fabs(outArray[j][0]-outArray[j+k][0]);
			dw = delW(rbar,outArray[j][7]);
			mass = outArray[j+k][3];
			rhoj = outArray[j+k][1];
			prj = outArray[j+k][4];
			
			//udot version
			//outArray[j][10] += -(mass/cs)*(pri/pow(rhoi,2.0))*(prj-pri)*dw;	
			outArray[j][10] += -(mass/cs)*(pri/rhoi)*(pri/pow(rhoi,2.0)+prj/pow(rhoj,2.0))*dw; //alternative formulation
			
			//gradU version
			outArray[j][9] += cs*(pri/pow(rhoi,2.0))*mass*dw;
		}

		outArray[j][11] = outArray[j][8]/outArray[j][10];
		outArray[j][12] = -outArray[j][8]/outArray[j][9];
	}
	
	printf("printing files...\n");
	
	//open the stream file
	snprintf(csvfile, sizeof(csvfile), "%s_edot.csv", argv[1]);
	FILE *stream, *fopen(), *vsz;
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");
	
	fprintf(stream,"x,rho,temp,mass,pr,cs,u,h,edot,gradU,crit,ratio1,ratio2\n");
	
	for(i = 0; i < j; i++)
	{
		fprintf(stream,"%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e,%5.4e\n", outArray[i][0],outArray[i][1],outArray[i][2],
				outArray[i][3],outArray[i][4],outArray[i][5],outArray[i][6],
				outArray[i][7],outArray[i][8],outArray[i][9],outArray[i][10],outArray[i][11],outArray[i][12]);
	}
	fclose(stream);
	
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
	
	free(body1);
	
	return 0;
}