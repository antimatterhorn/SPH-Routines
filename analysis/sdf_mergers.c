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
char csvfile[80];
char vszfile[80];
char pdffile[80];
char cmd[80];

int usage()
{
	printf("\t Creates a combination plot of omega, temp and a temp slice\n");
	printf("\t Usage: [required] {optional,default}\n");
	printf("\t sdf_combine [sdf file] [primary mass]\n");
	return 0;
}

int main(int argc, char **argv[])
{
	if (argc < 3)
	{
		usage();
		return 0;
	}
	
	char *cwd = getcwd(NULL, 0);	
	
	SDF *sdfp;
	SPHbody *body;
	
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
					NULL);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	
	printf("tpos = %3.5f\n",tpos);
	
	//	create the radial plot
	snprintf(cmd,sizeof(cmd),"sdf_plotrad %s",argv[1]);
	system(cmd);
	
	//	create the 2d temperature plot
	snprintf(cmd,sizeof(cmd),"sdf_2grid %s 2 1000 25 1e7",argv[1]);
	system(cmd);
	
	//	remove the individual plot files
	snprintf(cmd,sizeof(cmd),"rm %s*.vsz",argv[1]);
	system(cmd);

	//	create the combined plot file
	snprintf(vszfile, sizeof(vszfile), "%s_merger.vsz", argv[1]);
	
	FILE *fopen(), *vsz;
	vsz = fopen(vszfile,"w");
	
	fprintf(vsz,"AddImportPath('%s')\n",cwd);
	fprintf(vsz,"ImportFileCSV(u'%s_rad.csv', linked=True)\n",argv[1]);
	fprintf(vsz,"ImportFile2D(u'%s_temp.csv', [u'ds'], invertrows=False, invertcols=False, transpose=True, encoding='utf_8', linked=True)\n",argv[1]);
	fprintf(vsz,"Set('width', u'30cm')\n");
	fprintf(vsz,"Add('page', name='page1', autoadd=False)\n");
	fprintf(vsz,"To('page1')\n");
	fprintf(vsz,"Add('grid', name=u'grid1', autoadd=False)\n");
	fprintf(vsz,"To(u'grid1')\n");
	fprintf(vsz,"Set('columns', 1)\n");
	fprintf(vsz,"Set('leftMargin', '0cm')\n");
	fprintf(vsz,"Set('rightMargin', '15.3cm')\n");
	fprintf(vsz,"Set('topMargin', '0cm')\n");
	fprintf(vsz,"Set('bottomMargin', '0cm')\n");
	fprintf(vsz,"Add('graph', name='graph2', autoadd=False)\n");
	fprintf(vsz,"To('graph2')\n");
	fprintf(vsz,"Set('leftMargin', '2.14cm')\n");
	fprintf(vsz,"Set('rightMargin', '0cm')\n");
	fprintf(vsz,"Set('topMargin', u'0.15cm')\n");
	fprintf(vsz,"Set('bottomMargin', '0cm')\n");
	fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
	fprintf(vsz,"To('x')\n");
	fprintf(vsz,"Set('min', 0.0)\n");
	fprintf(vsz,"Set('log', False)\n");
	fprintf(vsz,"Set('datascale', 5.0000000000000003e-34)\n");
	fprintf(vsz,"Set('lowerPosition', 0.0)\n");
	fprintf(vsz,"Set('upperPosition', 1.0)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
	fprintf(vsz,"To('y')\n");
	fprintf(vsz,"Set('label', u'Temperature [K]')\n");
	fprintf(vsz,"Set('min', 2000000.0)\n");
	fprintf(vsz,"Set('log', False)\n");
	fprintf(vsz,"Set('direction', 'vertical')\n");
	fprintf(vsz,"Set('lowerPosition', 0.0)\n");
	fprintf(vsz,"Set('upperPosition', 1.0)\n");
	fprintf(vsz,"Set('otherPosition', 0.0)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy1', autoadd=False)\n");
	fprintf(vsz,"To('xy1')\n");
	fprintf(vsz,"Set('xData', u'mu')\n");
	fprintf(vsz,"Set('yData', u'temp')\n");
	fprintf(vsz,"Set('marker', u'circle')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('PlotLine/hide', True)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('function', name='function1', autoadd=False)\n");
	fprintf(vsz,"To('function1')\n");
	fprintf(vsz,"Set('function', u'%1.3f*2*10**33')\n",argv[2]);
	fprintf(vsz,"Set('variable', u'y')\n");
	fprintf(vsz,"Set('key', u'')\n");
	fprintf(vsz,"Set('Line/color', u'red')\n");
	fprintf(vsz,"Set('Line/width', u'2pt')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('graph', name='graph1', autoadd=False)\n");
	fprintf(vsz,"To('graph1')\n");
	fprintf(vsz,"Set('leftMargin', u'2.14cm')\n");
	fprintf(vsz,"Set('rightMargin', '0cm')\n");
	fprintf(vsz,"Set('topMargin', '0cm')\n");
	fprintf(vsz,"Set('bottomMargin', '1.48cm')\n");
	fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
	fprintf(vsz,"To('x')\n");
	fprintf(vsz,"Set('label', u'M(r) [M_{\\odot}]')\n");
	fprintf(vsz,"Set('min', 0.0)\n");
	fprintf(vsz,"Set('log', False)\n");
	fprintf(vsz,"Set('datascale', 5.0000000000000003e-34)\n");
	fprintf(vsz,"Set('lowerPosition', 0.0)\n");
	fprintf(vsz,"Set('upperPosition', 1.0)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
	fprintf(vsz,"To('y')\n");
	fprintf(vsz,"Set('label', u'\\Omega [s^{-1}]')\n");
	fprintf(vsz,"Set('min', 0.01)\n");
	fprintf(vsz,"Set('log', True)\n");
	fprintf(vsz,"Set('direction', 'vertical')\n");
	fprintf(vsz,"Set('lowerPosition', 0.0)\n");
	fprintf(vsz,"Set('upperPosition', 1.0)\n");
	fprintf(vsz,"Set('otherPosition', 0.0)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy1', autoadd=False)\n");
	fprintf(vsz,"To('xy1')\n");
	fprintf(vsz,"Set('xData', u'mu')\n");
	fprintf(vsz,"Set('yData', u'oz')\n");
	fprintf(vsz,"Set('marker', u'circle')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('key', u'SPH')\n");
	fprintf(vsz,"Set('PlotLine/width', u'2pt')\n");
	fprintf(vsz,"Set('PlotLine/hide', True)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy2', autoadd=False)\n");
	fprintf(vsz,"To('xy2')\n");
	fprintf(vsz,"Set('xData', u'mu')\n");
	fprintf(vsz,"Set('yData', u'ok')\n");
	fprintf(vsz,"Set('marker', u'none')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('key', u'Keplerian')\n");
	fprintf(vsz,"Set('PlotLine/color', u'red')\n");
	fprintf(vsz,"Set('PlotLine/width', u'3pt')\n");
	fprintf(vsz,"Set('PlotLine/hide', False)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('key', name='key1', autoadd=False)\n");
	fprintf(vsz,"To('key1')\n");
	fprintf(vsz,"Set('horzPosn', 'manual')\n");
	fprintf(vsz,"Set('vertPosn', 'manual')\n");
	fprintf(vsz,"Set('horzManual', 0.11168883058590126)\n");
	fprintf(vsz,"Set('vertManual', 0.099249445845623216)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('graph', name='graph3', autoadd=False)\n");
	fprintf(vsz,"To('graph3')\n");
	fprintf(vsz,"Set('leftMargin', '15.2cm')\n");
	fprintf(vsz,"Set('rightMargin', '0cm')\n");
	fprintf(vsz,"Set('topMargin', '0cm')\n");
	fprintf(vsz,"Set('bottomMargin', '0cm')\n");
	fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
	fprintf(vsz,"To('x')\n");
	fprintf(vsz,"Set('MajorTicks/hide', True)\n");
	fprintf(vsz,"Set('MinorTicks/hide', True)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
	fprintf(vsz,"To('y')\n");
	fprintf(vsz,"Set('direction', 'vertical')\n");
	fprintf(vsz,"Set('TickLabels/hide', True)\n");
	fprintf(vsz,"Set('MajorTicks/hide', True)\n");
	fprintf(vsz,"Set('MinorTicks/hide', True)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('label', name='label1', autoadd=False)\n");
	fprintf(vsz,"To('label1')\n");
	fprintf(vsz,"Set('label', u'%3.1f s')\n",tpos);
	fprintf(vsz,"Set('xPos', [0.050000000000000003])\n");
	fprintf(vsz,"Set('yPos', [0.90000000000000002])\n");
	fprintf(vsz,"Set('Text/size', u'18pt')\n");
	fprintf(vsz,"Set('Text/color', u'white')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('colorbar', name='colorbar1', autoadd=False)\n");
	fprintf(vsz,"To('colorbar1')\n");
	fprintf(vsz,"Set('image', u'image1')\n");
	fprintf(vsz,"Set('label', u'Temperature [K]')\n");
	fprintf(vsz,"Set('autoExtend', False)\n");
	fprintf(vsz,"Set('autoExtendZero', False)\n");
	fprintf(vsz,"Set('reflect', False)\n");
	fprintf(vsz,"Set('lowerPosition', 0.0)\n");
	fprintf(vsz,"Set('upperPosition', 1.0)\n");
	fprintf(vsz,"Set('otherPosition', 1.0)\n");
	fprintf(vsz,"Set('Label/size', u'18pt')\n");
	fprintf(vsz,"Set('Label/color', u'white')\n");
	fprintf(vsz,"Set('Label/atEdge', False)\n");
	fprintf(vsz,"Set('TickLabels/color', u'white')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('image', name='image1', autoadd=False)\n");
	fprintf(vsz,"To('image1')\n");
	fprintf(vsz,"Set('data', u'ds')\n");
	fprintf(vsz,"Set('min', 'Auto')\n");
	fprintf(vsz,"Set('max', 'Auto')\n");
	fprintf(vsz,"Set('colorScaling', u'linear')\n");
	fprintf(vsz,"Set('colorMap', u'spectrum3')\n");
	fprintf(vsz,"Set('colorInvert', True)\n");
	fprintf(vsz,"Set('smooth', True)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	
	
	
	fclose(vsz);
	
	return 0;
}