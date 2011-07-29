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


float x,rho,mass,pressure,h,cs,vx,tpos,drtpos,u,udot;
int i,j=0,npart=1,nion=13;
double g = 4.0/3.0;
double mass_in_g = 1.989E+27;        
double dist_in_cm = 6.955e7;  
double rcon = 8.314e7;
double time_in_s = 1;
double energy_in_erg,dens_in_gccm,pressure_in_ergperccm,specenergy_in_ergperg;
double temperature, density, ex_pep, ex_xne, ex_eta, krad, kion, kopac, kcond, energy,cv;
double s2rad, scond,abar,zbar;

double zarray[13]={ 2,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28};
double aarray[13]={ 4, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56};
double abund[13],molarabund[13];

int gnobj, nobj;
int conf;

char sdffile[80];
char csvfile[80];
char vszfile[80];
char pdffile[80];



#define Fortran(x) x

#define Fortran2(x) x##_

//extern "C" {
// input:
// temp   = temperature temp (in K)
// den    = density den (in g/cm**3)
// ionmax = number of isotopes in the composition
// bund   = mass fractions of the composition
// zion   = number of protons in each isotope (charge of each isotope)
// aion   = number of protons + neutrons in each isotope (atomic weight)
// pep    = electron-positron pressure (in erg/cm**3)
// xne    = electron-positron number density (in 1/cm**3)
// eta    = electron degeneracy parameter (chemical potential / k T)
//
// output:
// orad   = radiative opacity			->	krad
// ocond  = elecrton-ion opacity		->	kion
// opac   = total opacity (in cm**2/g)	->	kopac
// sigma  = total conductivity			->	kcond
// s2rad  = effective heat conductivity

void Fortran2(sig99)(double *temperature, double *density,double *bund, double *zion, double *aion, int *ionmax,
					 double *pep, double *xne, double *eta, double *orad, double *ocond, double *opac, 
					 double *s2rad, double *scond, double *sigma); 

void Fortran2(init_helm_table)();

void Fortran2(wrapper_helmeos)(int *npart, double *den_row, 
							   double *etot_row, double *abar_row, 
							   double *zbar_row, double *temperature,
							   double *pressure, double *pep, double *xne, double *eta);

void Fortran2(azbar)(double *xmass, double *aion, double *zion, int *ionmax, 
					 double *ymass, double *abar, double *zbar);
//}



// krad		= orad
// kion		= ocond
// kopac	= opac
// kcond	= sigma

int usage()
{
	printf("\t Creates a 2d plot of values pulling in sig99\n");
	printf("\t Usage: [required] {optional,default}\n");
	printf("\t sdf_plot2d [sdf file]\n");
	return 0;
}

int main(int argc, char **argv[])
{

	
	SDF *sdfp;
	SPHbody *body;
	
	double **outArray;
	
	outArray = (double **) malloc(nobj*sizeof(double *));
	for(i=0;i<j;i++){
		outArray[i]  = (double *) malloc(11*sizeof(double));
	}
	
	
	energy_in_erg = mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
	dens_in_gccm = mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
	pressure_in_ergperccm = energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	specenergy_in_ergperg = energy_in_erg/mass_in_g;
	
	Fortran2(init_helm_table)();
	
	
	
	if (argc < 2){
		usage();
		return 0;
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
	SDFgetfloatOrDefault(sdfp, "drtpos",  &drtpos, (float)0.0);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
	
	for(i = 0; i < nobj; i++){
		if ((fabs(body[i].y) < (body[i].h*2.0)) && (fabs(body[i].z) < (body[i].h*2.0))){
			j++;
		}
	}
	
	singlPrintf("Outputting %d points.\n",j);
	

	
	j=0;
	for(i = 0; i < nobj; i++){
		if ((fabs(body[i].y) < (body[i].h*2.0)) && (fabs(body[i].z) < (body[i].h*2.0))){
			rho			= body[i].rho * dens_in_gccm;
			mass		= body[i].mass * mass_in_g;
			pressure	= body[i].pr * pressure_in_ergperccm;
			x			= body[i].x * dist_in_cm;
			vx			= fabs(body[i].vx) * dist_in_cm;
			cs			= body[i].vsound * dist_in_cm;
			u			= body[i].u * specenergy_in_ergperg;
			abar		= body[i].abar;
			zbar		= body[i].zbar;
			temperature = body[i].temp;
			cv			= 3/2*rcon/abar;
			
			abund[0] = body[i].He4;
			abund[1] = body[i].C12;
			abund[2] = body[i].O16;
			abund[3] = body[i].Ne20;
			abund[4] = body[i].Mg24;
			abund[5] = body[i].Si28;
			abund[6] = body[i].S32;
			abund[7] = body[i].Ar36;
			abund[8] = body[i].Ca40;
			abund[9] = body[i].Ti44;
			abund[10] = body[i].Cr48;
			abund[11] = body[i].Fe52;
			abund[12] = body[i].Ni56;
			
			Fortran2(azbar)(abund, aarray, zarray, &nion, molarabund, &abar, &zbar );
			
						
			Fortran2(wrapper_helmeos)(&npart, &rho, &u, &abar, &zbar, &temperature, &pressure,
									  &ex_pep, &ex_xne, &ex_eta);
			
			Fortran2(sig99)(&temperature, &rho, abund, zarray, aarray, &nion, 
							&ex_pep, &ex_xne, &ex_eta,
							&krad, &kion, &kopac, &s2rad, &scond, &kcond);
			
			
			outArray[j][0] = x;
			outArray[j][1] = rho;
			outArray[j][2] = temperature;
			outArray[j][3] = mass;
			outArray[j][4] = pressure;
			outArray[j][5] = cs;
			outArray[j][6] = vx;
			outArray[j][7] = u;
			outArray[j][8] = body[i].h * dist_in_cm;
			outArray[j][9] = kcond;
			outArray[j][10] = kcond/(rho*cv);
			j++;
		}
	}
	
	//open the stream file
	snprintf(csvfile, sizeof(csvfile), "%s_2d.csv", argv[1]);
	FILE *stream, *fopen(), *vsz;
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");

	fprintf(stream,"x,rho,temp,mass,pr,cs,vx,u,h,udot,alpha\n");
	
	for(i = 0; i < j; i++)
	{
		fprintf(stream,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", outArray[i][0],outArray[i][1],outArray[i][2],
										outArray[i][3],outArray[i][4],outArray[i][5],outArray[i][6],
										outArray[i][7],outArray[i][8],outArray[i][9],outArray[i][10]);
	}
	fclose(stream);
	
	char *cwd = getcwd(NULL, 0);
	snprintf(vszfile, sizeof(vszfile), "%s_2d.vsz", argv[1]);
	snprintf(pdffile, sizeof(pdffile), "%s_2d.png", argv[1]);
	
	vsz = fopen(vszfile,"w");
	
//	fprintf(vsz,"ImportFileCSV(u'%s/%s', linked=True)\n",cwd,csvfile);
//	fprintf(vsz,"Add('page', name='page1', autoadd=False)\n");
//	fprintf(vsz,"To('page1')\n");
//	fprintf(vsz,"Add('graph', name='graph1', autoadd=False)\n");
//	fprintf(vsz,"To('graph1')\n");
//	fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
//	fprintf(vsz,"To('x')\n");
//	fprintf(vsz,"Set('label', u'x [cm]')\n");
//	fprintf(vsz,"Set('min', -1000000000.0)\n");
//	fprintf(vsz,"Set('max', 1000000000.0)\n");
//	fprintf(vsz,"Set('autoExtend', False)\n");
//	fprintf(vsz,"Set('lowerPosition', 1.0)\n");
//	fprintf(vsz,"Set('upperPosition', 0.0)\n");
//	fprintf(vsz,"Set('Label/size', u'20pt')\n");
//	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
//	fprintf(vsz,"To('..')\n");
//	fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
//	fprintf(vsz,"To('y')\n");
//	fprintf(vsz,"Set('label', u'rho [g cm^{-3}]')\n");
//	fprintf(vsz,"Set('min', 1e3)\n");
//	fprintf(vsz,"Set('max', 1e9)\n");
//	fprintf(vsz,"Set('log', True)\n");
//	fprintf(vsz,"Set('direction', 'vertical')\n");
//	fprintf(vsz,"Set('Label/size', u'20pt')\n");
//	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
//	fprintf(vsz,"To('..')\n");
//	fprintf(vsz,"Add('xy', name='xy1', autoadd=False)\n");
//	fprintf(vsz,"To('xy1')\n");
//	fprintf(vsz,"Set('yData', u'rho')\n");
//	fprintf(vsz,"Set('markerSize', u'1pt')\n");
//	fprintf(vsz,"Set('key', u'rho')\n");
//	fprintf(vsz,"Set('PlotLine/width', u'0.5pt')\n");
//	fprintf(vsz,"Set('PlotLine/hide', True)\n");
//	fprintf(vsz,"Set('MarkerLine/color', u'red')\n");
//	fprintf(vsz,"Set('MarkerFill/color', u'red')\n");
//	fprintf(vsz,"To('..')\n");
//	fprintf(vsz,"To('..')\n");
//	fprintf(vsz,"To('..')\n");
//	fprintf(vsz,"To('..')\n");
	
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
	fprintf(vsz,"Set('min', 0.0)\n");
	fprintf(vsz,"Set('max', 1000000000.0)\n");
	fprintf(vsz,"Set('autoExtend', False)\n");
	fprintf(vsz,"Set('lowerPosition', 0.0)\n");
	fprintf(vsz,"Set('upperPosition', 1.0)\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
	fprintf(vsz,"To('y')\n");
	fprintf(vsz,"Set('label', u'rho [g cm^{-3}], temp [K]')\n");
	fprintf(vsz,"Set('min', 100000.0)\n");
	fprintf(vsz,"Set('max', 10000000000.0)\n");
	fprintf(vsz,"Set('log', True)\n");
	fprintf(vsz,"Set('direction', 'vertical')\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('Label/color', u'red')\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"Set('TickLabels/color', u'red')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='axis1', autoadd=False)\n");
	fprintf(vsz,"To('axis1')\n");
	fprintf(vsz,"Set('label', u'u [erg g^{-1}]')\n");
	fprintf(vsz,"Set('min', 1e+15)\n");
	fprintf(vsz,"Set('max', 1e+18)\n");
	fprintf(vsz,"Set('log', True)\n");
	fprintf(vsz,"Set('direction', u'vertical')\n");
	fprintf(vsz,"Set('lowerPosition', 0.0)\n");
	fprintf(vsz,"Set('upperPosition', 1.0)\n");
	fprintf(vsz,"Set('otherPosition', 1.0)\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy1', autoadd=False)\n");
	fprintf(vsz,"To('xy1')\n");
	fprintf(vsz,"Set('yData', u'rho')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('key', u'rho')\n");
	fprintf(vsz,"Set('PlotLine/width', u'0.5pt')\n");
	fprintf(vsz,"Set('PlotLine/hide', True)\n");
	fprintf(vsz,"Set('MarkerLine/color', u'red')\n");
	fprintf(vsz,"Set('MarkerFill/color', u'red')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy2', autoadd=False)\n");
	fprintf(vsz,"To('xy2')\n");
	fprintf(vsz,"Set('yData', u'u')\n");
	fprintf(vsz,"Set('marker', u'circle')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('yAxis', u'axis1')\n");
	fprintf(vsz,"Set('PlotLine/hide', True)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy3', autoadd=False)\n");
	fprintf(vsz,"To('xy3')\n");
	fprintf(vsz,"Set('yData', u'temp')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('key', u'temp')\n");
	fprintf(vsz,"Set('PlotLine/width', u'0.5pt')\n");
	fprintf(vsz,"Set('PlotLine/hide', True)\n");
	fprintf(vsz,"Set('MarkerLine/color', u'red')\n");
	fprintf(vsz,"Set('MarkerFill/color', u'red')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('label', name='label1', autoadd=False)\n");
	fprintf(vsz,"To('label1')\n");
	fprintf(vsz,"Set('label', u't=%3.2fs, %3.2eyr')\n",tpos,drtpos/3.16e8);
	fprintf(vsz,"Set('xPos', [0.060957227525199997])\n");
	fprintf(vsz,"Set('yPos', [0.90471128755160013])\n");
	fprintf(vsz,"Set('Text/size', u'20pt')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	

	fclose(vsz);
	
	return 0;
	
	
	return 0;
}