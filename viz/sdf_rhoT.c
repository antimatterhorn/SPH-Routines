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
float rho,temp,mass,pressure,h,cs,vx,tpos,u;
double x,y,bx,by;
int i,j,k,jst,minid,id;
int binsy,magsy,minbiny,maxbiny,  binsx,magsx,minbinx,maxbinx;
double binx,biny;
double g = 4.0/3.0;
double mass_in_g = 1.989E+27;        
double dist_in_cm = 6.955e7;     
double time_in_s = 1;
double energy_in_erg,dens_in_gccm,pressure_in_ergperccm,specenergy_in_ergperg;

char sdffile[80];
char csv1file[80];
char co1file[80];
char si1file[80];
char ni1file[80];
char csv2file[80];
char co2file[80];
char si2file[80];
char ni2file[80];
char vszfile[80];
char pdffile[80];

int main(int argc, char **argv[])
{
	FILE *stream, *streamco, *streamsi, *streamni, *fopen(), *vsz;
	
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
	
	if (argc < 3){
		printf("minid: ");
		scanf(argv[2]);
	}
		
	minid = atoi(argv[2]);
	
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
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
	
	
//----------------binning section---------------------
	
	magsx	= 4;
	binsx	= 25;	//#bins per mag
	binx	= 1/(double)binsx;
	minbinx	= 4;
	maxbinx	= minbinx+magsx;
	
	magsy	= 4;
	binsy	= 25;	//#bins per mag
	biny	= 1/(double)binsy;
	minbiny	= 6;
	maxbiny	= minbiny+magsy;
	
//----------------fill arrays for first star-----------
	
	double **gridArray;
	double **coArray;
	double **siArray;
	double **niArray;
	
	gridArray = (double **) malloc((magsy*binsy)*sizeof(double *));
	for(i=0;i<(magsy*binsy);i++){
		gridArray[i]  = (double *) malloc((magsx*binsx)*sizeof(double));
	}	
	coArray = (double **) malloc((magsy*binsy)*sizeof(double *));
	for(i=0;i<(magsy*binsy);i++){
		coArray[i]  = (double *) malloc((magsx*binsx)*sizeof(double));
	}	
	siArray = (double **) malloc((magsy*binsy)*sizeof(double *));
	for(i=0;i<(magsy*binsy);i++){
		siArray[i]  = (double *) malloc((magsx*binsx)*sizeof(double));
	}
	niArray = (double **) malloc((magsy*binsy)*sizeof(double *));
	for(i=0;i<(magsy*binsy);i++){
		niArray[i]  = (double *) malloc((magsx*binsx)*sizeof(double));
	}
	
	for (i=0;i<magsy*binsy;i++)
	{
		for (j=0;j<magsx*binsx;j++) gridArray[i][j] = coArray[i][j] = siArray[i][j] = niArray[i][j] = 0;
	}
	
	for (k=0;k<nobj;k++)
	{
		x = log10(body[k].rho*dens_in_gccm);
		y = log10(body[k].temp);
		id = body[k].ident;
		
		for (i=0;i<magsy*binsy;i++)
		{
			by = i/(double)binsy+minbiny;
			for (j=0;j<magsx*binsx;j++)
			{
				bx = j/(double)binsx+minbinx;
				
				if ((x >= bx && x < bx+binx) && (y>=by && y<by+biny))	//inside the box
				{
					if (id < minid){
						gridArray[i][j] += 1;
						coArray[i][j] = 0.5*(body[k].C12 + body[k].O16);
						siArray[i][j] = body[k].Si28;
						niArray[i][j] = body[k].Ni56;
					} 					
				}				
			}
		}
	}

//--------------------stream first star--------------------
	
	snprintf(csv1file, sizeof(csv1file), "%s_rhoT1.csv", argv[1]);
	snprintf(co1file, sizeof(csv1file), "%s_co1.csv", argv[1]);
	snprintf(si1file, sizeof(csv1file), "%s_si1.csv", argv[1]);
	snprintf(ni1file, sizeof(csv1file), "%s_ni1.csv", argv[1]);
	
	stream = fopen(csv1file,"w");
	for (i=0;i<magsy*binsy;i++)
	{
		for (j=0;j<magsx*binsx;j++) fprintf(stream,"%1.2f ", gridArray[i][j]);
		fprintf(stream,"\n");
	}	
	fclose(stream);

	stream = fopen(co1file,"w");
	for (i=0;i<magsy*binsy;i++)
	{
		for (j=0;j<magsx*binsx;j++) {
			if (coArray[i][j] > siArray[i][j] && coArray[i][j] > niArray[i][j]){
				fprintf(stream,"%3.1f ",1.0);
			} else {
				fprintf(stream,"%3.1f ",0.0);
			}			
		}
		fprintf(stream,"\n");
	}	
	fclose(stream);
	
	stream = fopen(si1file,"w");
	for (i=0;i<magsy*binsy;i++)
	{
		for (j=0;j<magsx*binsx;j++) {
			if (siArray[i][j] > coArray[i][j] && siArray[i][j] > niArray[i][j]){
				fprintf(stream,"%3.1f ",1.0);
			} else {
				fprintf(stream,"%3.1f ",0.0);
			}			
		}
		fprintf(stream,"\n");
	}	
	fclose(stream);
	
	stream = fopen(ni1file,"w");
	for (i=0;i<magsy*binsy;i++)
	{
		for (j=0;j<magsx*binsx;j++) {
			if (niArray[i][j] > siArray[i][j] && niArray[i][j] > coArray[i][j]){
				fprintf(stream,"%3.1f ",1.0);
			} else {
				fprintf(stream,"%3.1f ",0.0);
			}			
		}
		fprintf(stream,"\n");
	}	
	fclose(stream);
	

//----------------fill arrays for second star-----------
		
	for (i=0;i<magsy*binsy;i++)
	{
		for (j=0;j<magsx*binsx;j++) gridArray[i][j] = coArray[i][j] = siArray[i][j] = niArray[i][j] = 0;
	}
	
	for (k=0;k<nobj;k++)
	{
		x = log10(body[k].rho*dens_in_gccm);
		y = log10(body[k].temp);
		id = body[k].ident;
		
		for (i=0;i<magsy*binsy;i++)
		{
			by = i/(double)binsy+minbiny;
			for (j=0;j<magsx*binsx;j++)
			{
				bx = j/(double)binsx+minbinx;
				
				if ((x >= bx && x < bx+binx) && (y>=by && y<by+biny))	//inside the box
				{
					if (id >= minid){
						gridArray[i][j] += 1;
						coArray[i][j] = 0.5*(body[k].C12 + body[k].O16);
						siArray[i][j] = body[k].Si28;
						niArray[i][j] = body[k].Ni56;
					} 					
				}				
			}
		}
	}
	
	
//--------------------stream second star--------------------
	
	snprintf(csv2file, sizeof(csv1file), "%s_rhoT2.csv", argv[1]);
	snprintf(co2file, sizeof(csv1file), "%s_co2.csv", argv[1]);
	snprintf(si2file, sizeof(csv1file), "%s_si2.csv", argv[1]);
	snprintf(ni2file, sizeof(csv1file), "%s_ni2.csv", argv[1]);
	
	stream = fopen(csv2file,"w");
	for (i=0;i<magsy*binsy;i++)
	{
		for (j=0;j<magsx*binsx;j++) fprintf(stream,"%1.2f ", gridArray[i][j]);
		fprintf(stream,"\n");
	}	
	fclose(stream);
	
	stream = fopen(co2file,"w");
	for (i=0;i<magsy*binsy;i++)
	{
		for (j=0;j<magsx*binsx;j++) {
			if (coArray[i][j] > siArray[i][j] && coArray[i][j] > niArray[i][j]){
				fprintf(stream,"%3.1f ",1.0);
			} else {
				fprintf(stream,"%3.1f ",0.0);
			}			
		}
		fprintf(stream,"\n");
	}	
	fclose(stream);
	
	stream = fopen(si2file,"w");
	for (i=0;i<magsy*binsy;i++)
	{
		for (j=0;j<magsx*binsx;j++) {
			if (siArray[i][j] > coArray[i][j] && siArray[i][j] > niArray[i][j]){
				fprintf(stream,"%3.1f ",1.0);
			} else {
				fprintf(stream,"%3.1f ",0.0);
			}			
		}
		fprintf(stream,"\n");
	}	
	fclose(stream);
	
	stream = fopen(ni2file,"w");
	for (i=0;i<magsy*binsy;i++)
	{
		for (j=0;j<magsx*binsx;j++) {
			if (niArray[i][j] > siArray[i][j] && niArray[i][j] > coArray[i][j]){
				fprintf(stream,"%3.1f ",1.0);
			} else {
				fprintf(stream,"%3.1f ",0.0);
			}			
		}
		fprintf(stream,"\n");
	}	
	fclose(stream);
	
	
//-------------------stream vsz file--------------------------
	
	//char *cwd = getcwd(NULL, 0);
	snprintf(vszfile, sizeof(vszfile), "%s_rhoT.vsz", argv[1]);
	snprintf(pdffile, sizeof(pdffile), "%s_rhoT.png", argv[1]);
	
	vsz = fopen(vszfile,"w");

	//fprintf(vsz,"AddImportPath(u'%s')\n",cwd);
	fprintf(vsz,"ImportFile2D(u'%s', [u'rhoT1'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",csv1file);
	fprintf(vsz,"ImportFile2D(u'%s', [u'co1'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",co1file);
	fprintf(vsz,"ImportFile2D(u'%s', [u'ni1'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",ni1file);
	fprintf(vsz,"ImportFile2D(u'%s', [u'si1'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",si1file);
	fprintf(vsz,"ImportFile2D(u'%s', [u'rhoT2'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",csv2file);
	fprintf(vsz,"ImportFile2D(u'%s', [u'co2'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",co2file);
	fprintf(vsz,"ImportFile2D(u'%s', [u'ni2'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",ni2file);
	fprintf(vsz,"ImportFile2D(u'%s', [u'si2'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",si2file);

	fprintf(vsz,"Set('height', u'25cm')\n");
	fprintf(vsz,"Add('page', name='page1', autoadd=False)\n");
	fprintf(vsz,"To('page1')\n");
	fprintf(vsz,"Add('grid', name='grid1', autoadd=False)\n");
	fprintf(vsz,"To('grid1')\n");
	fprintf(vsz,"Set('columns', 1)\n");
	fprintf(vsz,"Add('graph', name='graph1', autoadd=False)\n");
	fprintf(vsz,"To('graph1')\n");
	fprintf(vsz,"Set('leftMargin', u'1.59cm')\n");
	fprintf(vsz,"Set('rightMargin', '0.2cm')\n");
	fprintf(vsz,"Set('topMargin', '0.2cm')\n");
	fprintf(vsz,"Set('bottomMargin', '0cm')\n");
	fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
	fprintf(vsz,"To('x')\n");
	fprintf(vsz,"Set('hide', True)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
	fprintf(vsz,"To('y')\n");
	fprintf(vsz,"Set('hide', True)\n");
	fprintf(vsz,"Set('direction', 'vertical')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='axis1', autoadd=False)\n");
	fprintf(vsz,"To('axis1')\n");
	fprintf(vsz,"Set('label', u'Density [g/cc]')\n");
	fprintf(vsz,"Set('min', 10000.0)\n");
	fprintf(vsz,"Set('max', 100000000.0)\n");
	fprintf(vsz,"Set('log', True)\n");
	fprintf(vsz,"Set('autoExtend', True)\n");
	fprintf(vsz,"Set('autoExtendZero', True)\n");
	fprintf(vsz,"Set('autoMirror', True)\n");
	fprintf(vsz,"Set('reflect', False)\n");
	fprintf(vsz,"Set('outerticks', False)\n");
	fprintf(vsz,"Set('lowerPosition', 0.0)\n");
	fprintf(vsz,"Set('upperPosition', 1.0)\n");
	fprintf(vsz,"Set('otherPosition', 0.0)\n");
	fprintf(vsz,"Set('Label/size', u'18pt')\n");
	fprintf(vsz,"Set('Label/hide', True)\n");
	fprintf(vsz,"Set('TickLabels/size', u'16pt')\n");
	fprintf(vsz,"Set('TickLabels/hide', True)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='axis2', autoadd=False)\n");
	fprintf(vsz,"To('axis2')\n");
	fprintf(vsz,"Set('label', u'Temperature [K]')\n");
	fprintf(vsz,"Set('min', 1000000.0)\n");
	fprintf(vsz,"Set('max', 10000000000.0)\n");
	fprintf(vsz,"Set('log', True)\n");
	fprintf(vsz,"Set('direction', u'vertical')\n");
	fprintf(vsz,"Set('Label/size', u'18pt')\n");
	fprintf(vsz,"Set('TickLabels/size', u'16pt')\n");
	fprintf(vsz,"To('..')\n");
    fprintf(vsz,"Add('label', name='label1', autoadd=False)\n");
    fprintf(vsz,"To('label1')\n");
    fprintf(vsz,"Set('label', u'%3.2f')\n",tpos);
    fprintf(vsz,"Set('xPos', [0.5155758712508223])\n");
    fprintf(vsz,"Set('yPos', [0.07556142227156612])\n");
    fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('image', name=u'Ni1', autoadd=False)\n");
	fprintf(vsz,"To(u'Ni1')\n");
	fprintf(vsz,"Set('data', u'rhoT1')\n");
	fprintf(vsz,"Set('colorScaling', u'log')\n");
	fprintf(vsz,"Set('transparencyData', u'ni1')\n");
	fprintf(vsz,"Set('colorMap', u'blue')\n");
	fprintf(vsz,"Set('colorInvert', True)\n");
	fprintf(vsz,"Set('transparency', 0)\n");
	fprintf(vsz,"Set('smooth', False)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('image', name=u'Si1', autoadd=False)\n");
	fprintf(vsz,"To(u'Si1')\n");
	fprintf(vsz,"Set('data', u'rhoT1')\n");
	fprintf(vsz,"Set('colorScaling', u'log')\n");
	fprintf(vsz,"Set('transparencyData', u'si1')\n");
	fprintf(vsz,"Set('colorMap', u'red')\n");
	fprintf(vsz,"Set('colorInvert', True)\n");
	fprintf(vsz,"Set('smooth', False)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('image', name=u'CO1', autoadd=False)\n");
	fprintf(vsz,"To(u'CO1')\n");
	fprintf(vsz,"Set('data', u'rhoT1')\n");
	fprintf(vsz,"Set('colorScaling', u'log')\n");
	fprintf(vsz,"Set('transparencyData', u'co1')\n");
	fprintf(vsz,"Set('colorMap', u'green')\n");
	fprintf(vsz,"Set('colorInvert', True)\n");
	fprintf(vsz,"Set('smooth', False)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	
	
	fprintf(vsz,"Add('graph', name='graph2', autoadd=False)\n");
	fprintf(vsz,"To('graph2')\n");
	fprintf(vsz,"Set('leftMargin', '1.59cm')\n");
	fprintf(vsz,"Set('rightMargin', '0.2cm')\n");
	fprintf(vsz,"Set('topMargin', '0cm')\n");
	fprintf(vsz,"Set('bottomMargin', '0.61cm')\n");
	fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
	fprintf(vsz,"To('x')\n");
	fprintf(vsz,"Set('hide', True)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
	fprintf(vsz,"To('y')\n");
	fprintf(vsz,"Set('hide', True)\n");
	fprintf(vsz,"Set('direction', 'vertical')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='axis1', autoadd=False)\n");
	fprintf(vsz,"To('axis1')\n");
	fprintf(vsz,"Set('label', u'Density [g/cc]')\n");
	fprintf(vsz,"Set('min', 10000.0)\n");
	fprintf(vsz,"Set('max', 100000000.0)\n");
	fprintf(vsz,"Set('log', True)\n");
	fprintf(vsz,"Set('autoExtend', True)\n");
	fprintf(vsz,"Set('autoExtendZero', True)\n");
	fprintf(vsz,"Set('autoMirror', True)\n");
	fprintf(vsz,"Set('reflect', False)\n");
	fprintf(vsz,"Set('outerticks', False)\n");
	fprintf(vsz,"Set('Label/size', u'18pt')\n");
	fprintf(vsz,"Set('TickLabels/size', u'16pt')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='axis2', autoadd=False)\n");
	fprintf(vsz,"To('axis2')\n");
	fprintf(vsz,"Set('label', u'Temperature [K]')\n");
	fprintf(vsz,"Set('min', 1000000.0)\n");
	fprintf(vsz,"Set('max', 10000000000.0)\n");
	fprintf(vsz,"Set('log', True)\n");
	fprintf(vsz,"Set('direction', u'vertical')\n");
	fprintf(vsz,"Set('Label/size', u'18pt')\n");
	fprintf(vsz,"Set('TickLabels/size', u'16pt')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('image', name=u'Ni2', autoadd=False)\n");
	fprintf(vsz,"To(u'Ni2')\n");
	fprintf(vsz,"Set('data', u'rhoT2')\n");
	fprintf(vsz,"Set('colorScaling', u'log')\n");
	fprintf(vsz,"Set('transparencyData', u'ni2')\n");
	fprintf(vsz,"Set('colorMap', u'blue')\n");
	fprintf(vsz,"Set('colorInvert', True)\n");
	fprintf(vsz,"Set('transparency', 0)\n");
	fprintf(vsz,"Set('smooth', False)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('image', name=u'Si2', autoadd=False)\n");
	fprintf(vsz,"To(u'Si2')\n");
	fprintf(vsz,"Set('data', u'rhoT2')\n");
	fprintf(vsz,"Set('colorScaling', u'log')\n");
	fprintf(vsz,"Set('transparencyData', u'si2')\n");
	fprintf(vsz,"Set('colorMap', u'red')\n");
	fprintf(vsz,"Set('colorInvert', True)\n");
	fprintf(vsz,"Set('smooth', False)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('image', name=u'CO2', autoadd=False)\n");
	fprintf(vsz,"To(u'CO2')\n");
	fprintf(vsz,"Set('data', u'rhoT2')\n");
	fprintf(vsz,"Set('colorScaling', u'log')\n");
	fprintf(vsz,"Set('transparencyData', u'co2')\n");
	fprintf(vsz,"Set('colorMap', u'green')\n");
	fprintf(vsz,"Set('colorInvert', True)\n");
	fprintf(vsz,"Set('smooth', False)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	
	
	
	
	

	fclose(vsz);

	free(body);
	free(gridArray);
	free(coArray);
	free(siArray);
	free(niArray);
	
	return 0;
}