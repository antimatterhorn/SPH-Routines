/*
 *  sdf2csv.c
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "sdf_plot.h"
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
float x,rho,temp,mass,pressure,h,cs,vx,tpos;
int i,j=0,k;
double g = 4.0/3.0;
double mass_in_g = 1.989E+27;        
double dist_in_cm = 6.955e7;     
double time_in_s = 1;
double energy_in_erg,dens_in_gccm,pressure_in_ergperccm,specenergy_in_ergperg;

char sdffile[80];
char csvfile[80];
char vszfile[80];
char pdffile[80];

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
					"bdot", offsetof(SPHbody, bdot), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
	
	for(i = 0; i < nobj; i++){
		if ((fabs(body[i].y) < (body[i].h*2.0)) && (fabs(body[i].z) < (body[i].h*2.0))){
			j++;
		}
	}
	
	singlPrintf("Outputting %d points.\n",j);
	int np = j;
	int na = 10;
	
	double **outArray;
	double tempArray[na];
	
	outArray = (double **) malloc(nobj*sizeof(double *));
	for(i=0;i<j;i++){
		outArray[i]  = (double *) malloc(na*sizeof(double));
	}
	
	j=0;
	for(i = 0; i < nobj; i++){
		if ((fabs(body[i].y) < (body[i].h*2.0)) && (fabs(body[i].z) < (body[i].h*2.0))){
			//singlPrintf("h=%f y=%f z=%f\n",body[i].h,body[i].y,body[i].z);
			rho = body[i].rho * dens_in_gccm;
			mass = body[i].mass * mass_in_g;
			pressure = body[i].pr * pressure_in_ergperccm;
			x = body[i].x * dist_in_cm;
			cs = sqrt(g*pressure/rho);
			vx = fabs(body[i].vx) * dist_in_cm;
			//r = sqrt(pow(body[i].x,2.0)+pow(body[i].y,2.0)+pow(body[i].z,2.0));
			outArray[j][0] = x;
			outArray[j][1] = rho;
			outArray[j][2] = body[i].temp;
			outArray[j][3] = mass;
			outArray[j][4] = pressure;
			outArray[j][5] = cs;
			outArray[j][6] = vx;
			outArray[j][7] = body[i].h * dist_in_cm;
			outArray[j][8] = body[i].bdot;
			outArray[j][9] = body[i].u * specenergy_in_ergperg;
			j++;
		}
	}
	
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
	
	//open the stream file
	snprintf(csvfile, sizeof(csvfile), "%s_cs.csv", argv[1]);
	FILE *stream, *fopen(), *vsz;
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");

	fprintf(stream,"x,rho,temp,mass,pr,cs,vx,h,bdot,u\n");
	
	for(i = 0; i < j; i++)
	{
		fprintf(stream,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", outArray[i][0],outArray[i][1],outArray[i][2],
										outArray[i][3],outArray[i][4],outArray[i][5],outArray[i][6],outArray[i][7],outArray[i][8],outArray[i][9]);
	}
	fclose(stream);
	
	char *cwd = getcwd(NULL, 0);
	snprintf(vszfile, sizeof(vszfile), "%s_cs.vsz", argv[1]);
	snprintf(pdffile, sizeof(pdffile), "%s_cs.png", argv[1]);
	
	vsz = fopen(vszfile,"w");
	
	fprintf(vsz,"ImportFileCSV(u'%s/%s', linked=True)\n",cwd,csvfile);
	fprintf(vsz,"Add('page', name='page1', autoadd=False)\n");
	fprintf(vsz,"To('page1')\n");
	fprintf(vsz,"Add('graph', name='graph1', autoadd=False)\n");
	fprintf(vsz,"To('graph1')\n");
	fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
	fprintf(vsz,"To('x')\n");
	fprintf(vsz,"Set('label', u'x [cm]')\n");
	fprintf(vsz,"Set('min', -1000000000.0)\n");
	fprintf(vsz,"Set('max', 1000000000.0)\n");
	fprintf(vsz,"Set('autoExtend', False)\n");
	fprintf(vsz,"Set('lowerPosition', 1.0)\n");
	fprintf(vsz,"Set('upperPosition', 0.0)\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
	fprintf(vsz,"To('y')\n");
	fprintf(vsz,"Set('label', u'velocity [cm s^{-1}]')\n");
	fprintf(vsz,"Set('min', 50000000.0)\n");
	fprintf(vsz,"Set('max', 1000000000.0)\n");
	fprintf(vsz,"Set('log', True)\n");
	fprintf(vsz,"Set('direction', 'vertical')\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy1', autoadd=False)\n");
	fprintf(vsz,"To('xy1')\n");
	fprintf(vsz,"Set('yData', u'cs')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('key', u'cs')\n");
	fprintf(vsz,"Set('PlotLine/width', u'0.5pt')\n");
	fprintf(vsz,"Set('PlotLine/hide', True)\n");
	fprintf(vsz,"Set('MarkerLine/color', u'red')\n");
	fprintf(vsz,"Set('MarkerFill/color', u'red')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy2', autoadd=False)\n");
	fprintf(vsz,"To('xy2')\n");
	fprintf(vsz,"Set('yData', u'vx')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('key', u'vx')\n");
	fprintf(vsz,"Set('PlotLine/hide', True)\n");
	fprintf(vsz,"Set('MarkerLine/color', u'blue')\n");
	fprintf(vsz,"Set('MarkerFill/color', u'blue')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('key', name='key1', autoadd=False)\n");
	fprintf(vsz,"Add('label', name='label1', autoadd=False)\n");
	fprintf(vsz,"To('label1')\n");
	fprintf(vsz,"Set('label', u't=%f')\n",tpos);
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