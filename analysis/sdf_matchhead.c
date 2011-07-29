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
float x,rho,temp,mass,pressure,h,cs,vx,t1,t0,dt,u,edot;
double v,rbar,pri,prj,rhoi,rhoj,dw;
int i,j=0,k;
int nbrs=5;
int np;
double g = 4.0/3.0;
double mass_in_g = 1.989E+27;        
double dist_in_cm = 6.955e7;     
double time_in_s = 1;
double energy_in_erg,dens_in_gccm,pressure_in_ergperccm,specenergy_in_ergperg;

char sdffile[80];
char csvfile[80];
char vszfile[80];
char pdffile[80];

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
	
	SDF *sdfp1,*sdfp0;
	SPHbody *body1,*body0;
	
	if (argc < 3){
		printf("current SDF file: ");
		gets (argv[1]);
		printf("previous SDF file: ");
		gets (argv[2]);
	}
	
	if (argc > 3) nbrs = atoi(argv[3]);
	
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
	SDFgetfloatOrDefault(sdfp1, "tpos", &t1, (float)0.0);
	
	sdfp0 = SDFreadf(argv[2], (void **)&body0, &gnobj, &nobj, sizeof(SPHbody),
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
	SDFgetfloatOrDefault(sdfp0, "tpos",  &t0, (float)0.0);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
	
	for(i = 0; i < nobj; i++){
		if ((fabs(body1[i].y) < (body1[i].h*2.0)) && (fabs(body1[i].z) < (body1[i].h*2.0))){
			j++;
		}
	}
	
	np = j;
	
	dt = t1-t0;
	
	singlPrintf("Outputting %d points.\n",j);
	
	double **outArray;
	double tempArray[14];
	
	outArray = (double **) malloc(nobj*sizeof(double *));
	for(i=0;i<j;i++){
		outArray[i]  = (double *) malloc(14*sizeof(double));
	}
	
	//0,  1,   2,   3, 4, 5, 6,7,8, 9,  10,  11,12,  13
	//x,rho,temp,mass,pr,cs,vx,u,h,dW,edot,crit,id,ratio
	
	printf("filling array...\n");
	
	j=0;
	for(i = 0; i < nobj; i++){
		if ((fabs(body1[i].y) < (body1[i].h*2.0)) && (fabs(body1[i].z) < (body1[i].h*2.0))){
			rho = body1[i].rho * dens_in_gccm;
			mass = body1[i].mass * mass_in_g;
			pressure = body1[i].pr * pressure_in_ergperccm;
			x = body1[i].x * dist_in_cm;
			vx = fabs(body1[i].vx) * dist_in_cm;
			cs = body1[i].vsound * dist_in_cm;
			u = body1[i].u * specenergy_in_ergperg;
			outArray[j][0] = x;
			outArray[j][1] = rho;
			outArray[j][2] = body1[i].temp;
			outArray[j][3] = mass;
			outArray[j][4] = pressure;
			outArray[j][5] = cs;
			outArray[j][6] = vx;
			outArray[j][7] = u;
			outArray[j][8] = body1[i].h * dist_in_cm;
			outArray[j][9] = outArray[j][11] = outArray[j][13] = 0;
			outArray[j][10] = body1[i].udot * specenergy_in_ergperg;
			outArray[j][12] = body1[i].ident;
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
				for (k=0;k<14;k++) 
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
		//first get edot
		for(i = 0; i < nobj; i++)
		{
			if (outArray[j][12] = body0[i].ident) outArray[j][10] = (outArray[j][7] - body0[i].u * specenergy_in_ergperg)/dt;
			if (outArray[j][10] < 0) outArray[j][10] = 0;
		}
		
		//now calc rbar and crit condition
		//if (j>0 && j<(np-1))
//		{
//			rbar = 0.5*(fabs(outArray[j][0]-outArray[j-1][0])+fabs(outArray[j][0]-outArray[j+1][0]));
//			dw = delW(rbar,outArray[j][8]);
//			rhoi = outArray[j][1];
//			rhoj = 0.5*(outArray[j-1][1]+outArray[j+1][1]);
//			pri = outArray[j][4];
//			prj = 0.5*(outArray[j-1][4]+outArray[j+1][4]);
//			mass = 0.5*(outArray[j-1][3]+outArray[j+1][3]);
//			cs = outArray[j][5];
//		} else if (j==0) 
//		{
//			rbar = fabs(outArray[j][0]-outArray[j+1][0]);
//			dw = delW(rbar,outArray[j][8]);
//			rhoi = outArray[j][1];
//			rhoj = outArray[j+1][1];
//			pri = outArray[j][4];
//			prj = outArray[j+1][4];
//			mass = outArray[j+1][3];
//			cs = outArray[j][5];
//			
//		} else if (j==(np-1))
//		{
//			rbar = fabs(outArray[j][0]-outArray[j-1][0]);
//			dw = delW(rbar,outArray[j][8]);
//			rhoi = outArray[j][1];
//			rhoj = outArray[j-1][1];
//			pri = outArray[j][4];
//			prj = outArray[j-1][4];
//			mass = outArray[j-1][3];
//			cs = outArray[j][5];
//			
//		}
		
		rhoi = outArray[j][1];
		pri = outArray[j][4];
		cs = outArray[j][5];
		
		for(k=-nbrs;k<(nbrs+1);k++)
		{
			rbar = fabs(outArray[j][0]-outArray[j+k][0]);
			dw = delW(rbar,outArray[j][8]);
			mass = outArray[j+k][3];
			rhoj = outArray[j+k][1];
			prj = outArray[j+k][4];
			//outArray[j][11] += -(mass/cs)*(pri/pow(rhoi,2.0))*(prj-pri)*dw;
			outArray[j][11] += -(mass/cs)*(pri/rhoi)*(pri/pow(rhoi,2.0)+prj/pow(rhoj,2.0))*dw; //alternative formulation
		}
		
		//outArray[j][9] = dw;
		//outArray[j][11] = ((mass/cs)*(pri/rhoi)*(pri/pow(rhoi,2.0)+prj/pow(rhoj,2.0))*dw);
		outArray[j][13] = outArray[j][10]/outArray[j][11];
	}
	
	printf("outputting...\n");
	
	//open the stream file
	snprintf(csvfile, sizeof(csvfile), "%s_cs.csv", argv[1]);
	FILE *stream, *fopen(), *vsz;
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");
	
	fprintf(stream,"x,rho,temp,mass,pr,cs,vx,u,h,dW,edot,crit,ratio\n");
	
	for(i = 0; i < j; i++)
	{
		fprintf(stream,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", outArray[i][0],outArray[i][1],outArray[i][2],
				outArray[i][3],outArray[i][4],outArray[i][5],outArray[i][6],
				outArray[i][7],outArray[i][8],outArray[i][9],outArray[i][10],outArray[i][11],outArray[i][13]);
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
	fprintf(vsz,"Set('leftMargin', '1.7cm')\n");
	fprintf(vsz,"Set('rightMargin', '1.61cm')\n");
	fprintf(vsz,"Set('topMargin', '0.276cm')\n");
	fprintf(vsz,"Set('bottomMargin', '1.71cm')\n");
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
	fprintf(vsz,"Set('label', u'rho [g cm^{-3}]')\n");
	fprintf(vsz,"Set('min', 1000.0)\n");
	fprintf(vsz,"Set('max', 1000000000.0)\n");
	fprintf(vsz,"Set('log', True)\n");
	fprintf(vsz,"Set('direction', 'vertical')\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('Label/color', u'red')\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"Set('TickLabels/color', u'red')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='axis1', autoadd=False)\n");
	fprintf(vsz,"To('axis1')\n");
	fprintf(vsz,"Set('label', u'Temp [k]')\n");
	fprintf(vsz,"Set('min', 100000.0)\n");
	fprintf(vsz,"Set('max', 10000000000.0)\n");
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
	fprintf(vsz,"Set('yData', u'temp')\n");
	fprintf(vsz,"Set('marker', u'circle')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('yAxis', u'axis1')\n");
	fprintf(vsz,"Set('PlotLine/hide', True)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	
	
	fclose(vsz);
	
	return 0;
	
	
	return 0;
}