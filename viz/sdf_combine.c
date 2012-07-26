/*
 *  sdf_combine.c
 *  creates a combination plot for a given output and particle id cutoff
 *
 *
 */

#include "sdf_2grid.h"
#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
#include "Msgs.h"
#include <stdlib.h>
#include <math.h>

int gnobj, nobj;
int conf,id;
int i,j;

float tpos;
float m1,m2;

char rmvsz[80];
char cmd[80];
char vszfile[80];

int usage()
{
	printf("\t Creates a combination plot of rho-Temp, vx-cs and density\n");
	printf("\t Usage: [required] {optional,default}\n");
	printf("\t sdf_combine [sdf file] [id] [mass1] [mass2]\n");
	return 0;
}


int main(int argc, char **argv[])
{
	
	if (argc < 5)
	{
		usage();
		return 0;
	}
	
	id = atoi(argv[2]);
	m1 = atof(argv[3]);
	m2 = atof(argv[4]);
	
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
					//					"abar", offsetof(SPHbody, abar), &conf,
					//					"zbar", offsetof(SPHbody, zbar), &conf,
					//					"ax", offsetof(SPHbody, ax), &conf,
					//					"ay", offsetof(SPHbody, ay), &conf,
					//					"az", offsetof(SPHbody, az), &conf,
					//					"lax", offsetof(SPHbody, lax), &conf,
					//					"lay", offsetof(SPHbody, lay), &conf,
					//					"laz", offsetof(SPHbody, laz), &conf,
					//					"gax", offsetof(SPHbody, gax), &conf,
					//					"gay", offsetof(SPHbody, gay), &conf,
					//					"gaz", offsetof(SPHbody, gaz), &conf,
					//					"grav_mass", offsetof(SPHbody, grav_mass), &conf,
					//					"phi", offsetof(SPHbody, phi), &conf,
					//					"tacc", offsetof(SPHbody, tacc), &conf,
					//					"idt", offsetof(SPHbody, idt), &conf,
					//					"nbrs", offsetof(SPHbody, nbrs), &conf,
					//					"ident", offsetof(SPHbody, ident), &conf,
					//					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	
	printf("tpos = %3.5f\n",tpos);

	//	create the sound speed plot
	snprintf(cmd,sizeof(cmd),"sdf_plotcs %s",argv[1]);
	system(cmd);
	
	//	create the 2d density plot
	snprintf(cmd,sizeof(cmd),"sdf_gridp %s 1 1000 14.38 1e7",argv[1]);
	system(cmd);
	
	//	create the rho-T plots
	snprintf(cmd,sizeof(cmd),"sdf_rhoT %s %d",argv[1],id);
	system(cmd);
	
	//	remove the individual plot files
	snprintf(rmvsz,sizeof(rmvsz),"rm %s*.vsz",argv[1]);
	system(rmvsz);
	
	//	create the combined plot file
	snprintf(vszfile, sizeof(vszfile), "%s_comb.vsz", argv[1]);
	
	FILE *fopen(), *vsz;
	vsz = fopen(vszfile,"w");

//	fprintf(vsz,"ImportFileCSV(u'%s_rhoT.csv', linked=True)\n",argv[1]);
//	fprintf(vsz,"ImportFileCSV(u'%s_cs.csv', linked=True, dsprefix=u'cs_')\n",argv[1]);
//	fprintf(vsz,"ImportFileCSV(u'%s_rhoT2.csv', linked=True, dsprefix=u'2_')\n",argv[1]);
//	fprintf(vsz,"ImportFile2D(u'%s.csv', [u'ds'], invertrows=False, invertcols=False, transpose=True, encoding='utf_8', linked=True)\n",argv[1]);
	
	fprintf(vsz,"ImportFile2D(u'%s_rhoT1.csv', [u'rhoT1'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",argv[1]);
	fprintf(vsz,"ImportFile2D(u'%s_co1.csv', [u'co1'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",argv[1]);
	fprintf(vsz,"ImportFile2D(u'%s_ni1.csv', [u'ni1'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",argv[1]);
	fprintf(vsz,"ImportFile2D(u'%s_si1.csv', [u'si1'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",argv[1]);
	fprintf(vsz,"ImportFile2D(u'%s_rhoT2.csv', [u'rhoT2'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",argv[1]);
	fprintf(vsz,"ImportFile2D(u'%s_co2.csv', [u'co2'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",argv[1]);
	fprintf(vsz,"ImportFile2D(u'%s_ni2.csv', [u'ni2'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",argv[1]);
	fprintf(vsz,"ImportFile2D(u'%s_si2.csv', [u'si2'], invertrows=True, invertcols=False, transpose=False, encoding='utf_8', linked=True)\n",argv[1]);
	fprintf(vsz,"ImportFileCSV(u'%s_cs.csv', linked=True, dsprefix=u'cs_')\n",argv[1]);
	fprintf(vsz,"ImportFile2D(u'%s.csv', [u'ds'], invertrows=False, invertcols=False, transpose=True, encoding='utf_8', linked=True)\n",argv[1]);
	
	fprintf(vsz,"Set('width', u'25cm')\n");
	fprintf(vsz,"Set('height', u'20cm')\n");
	fprintf(vsz,"Add('page', name='page1', autoadd=False)\n");
	fprintf(vsz,"To('page1')\n");
	fprintf(vsz,"Add('grid', name='grid1', autoadd=False)\n");
	fprintf(vsz,"To('grid1')\n");
	fprintf(vsz,"Set('columns', 2)\n");
	fprintf(vsz,"Set('scaleRows', [1.0, 1.125])\n");
	fprintf(vsz,"Set('scaleCols', [1.0, 1.25])\n");
	fprintf(vsz,"Set('leftMargin', '0cm')\n");
	fprintf(vsz,"Set('rightMargin', '0cm')\n");
	fprintf(vsz,"Set('topMargin', '0cm')\n");
	fprintf(vsz,"Set('bottomMargin', '0cm')\n");
	fprintf(vsz,"Add('graph', name='graph1', autoadd=False)\n");
	fprintf(vsz,"To('graph1')\n");
	fprintf(vsz,"Set('leftMargin', '2.09cm')\n");
	fprintf(vsz,"Set('rightMargin', '0.169cm')\n");
	fprintf(vsz,"Set('topMargin', u'1.276cm')\n");
	fprintf(vsz,"Set('bottomMargin', u'0')\n");
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
	fprintf(vsz,"Set('label', u'M = %1.2f M_{\\odot}')\n",m1);
	fprintf(vsz,"Set('Text/size', u'18pt')\n");
	fprintf(vsz,"Set('xPos', [0.039439388406118166])\n");
	fprintf(vsz,"Set('yPos', [0.90394296957670195])\n");
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
	fprintf(vsz,"Set('leftMargin', u'3.53cm')\n");
	fprintf(vsz,"Set('rightMargin', '0.1cm')\n");
	fprintf(vsz,"Set('topMargin', u'1.276cm')\n");
	fprintf(vsz,"Set('bottomMargin', u'0cm')\n");
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
	fprintf(vsz,"Set('label', u'[cm/s], [g/cc], [K]')\n");
	fprintf(vsz,"Set('min', 1000000.0)\n");
	fprintf(vsz,"Set('max', 10000000000.0)\n");
	fprintf(vsz,"Set('log', True)\n");
	fprintf(vsz,"Set('autoExtend', False)\n");
	fprintf(vsz,"Set('autoExtendZero', False)\n");
	fprintf(vsz,"Set('direction', 'vertical')\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"Set('TickLabels/format', u'Auto')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('key', name='key1', autoadd=False)\n");
	fprintf(vsz,"To('key1')\n");
	fprintf(vsz,"Set('horzPosn', 'manual')\n");
	fprintf(vsz,"Set('vertPosn', 'manual')\n");
	fprintf(vsz,"Set('horzManual', 0.68175512401632088)\n");
	fprintf(vsz,"Set('vertManual', 0.60688585314552079)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy1', autoadd=False)\n");
	fprintf(vsz,"To('xy1')\n");
	fprintf(vsz,"Set('xData', u'cs_x')\n");
	fprintf(vsz,"Set('yData', u'cs_cs')\n");
	fprintf(vsz,"Set('marker', u'circle')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('key', u'c_{s}')\n");
	fprintf(vsz,"Set('PlotLine/color', u'red')\n");
	fprintf(vsz,"Set('PlotLine/width', u'2pt')\n");
	fprintf(vsz,"Set('PlotLine/hide', True)\n");
	fprintf(vsz,"Set('MarkerLine/color', u'red')\n");
	fprintf(vsz,"Set('MarkerFill/color', u'red')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy2', autoadd=False)\n");
	fprintf(vsz,"To('xy2')\n");
	fprintf(vsz,"Set('xData', u'cs_x')\n");
	fprintf(vsz,"Set('yData', u'cs_vx')\n");
	fprintf(vsz,"Set('marker', u'circle')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('key', u'v_{x}')\n");
	fprintf(vsz,"Set('PlotLine/hide', True)\n");
	fprintf(vsz,"Set('MarkerLine/color', u'blue')\n");
	fprintf(vsz,"Set('MarkerFill/color', u'blue')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy4', autoadd=False)\n");
	fprintf(vsz,"To('xy4')\n");
	fprintf(vsz,"Set('xData', u'cs_x')\n");
	fprintf(vsz,"Set('yData', u'cs_rho')\n");
	fprintf(vsz,"Set('marker', u'none')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('key', u'rho')\n");
	fprintf(vsz,"Set('PlotLine/color', u'green')\n");
	fprintf(vsz,"Set('PlotLine/width', u'2pt')\n");
	fprintf(vsz,"Set('PlotLine/hide', False)\n");
	fprintf(vsz,"Set('MarkerLine/color', u'green')\n");
	fprintf(vsz,"Set('MarkerFill/color', u'green')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy3', autoadd=False)\n");
	fprintf(vsz,"To('xy3')\n");
	fprintf(vsz,"Set('xData', u'cs_x')\n");
	fprintf(vsz,"Set('yData', u'cs_temp')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('key', u'T')\n");
	fprintf(vsz,"Set('PlotLine/hide', True)\n");
	fprintf(vsz,"Set('MarkerLine/color', u'magenta')\n");
	fprintf(vsz,"Set('MarkerFill/color', u'red')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	
	fprintf(vsz,"Add('graph', name='graph3', autoadd=False)\n");
	fprintf(vsz,"To('graph3')\n");
	fprintf(vsz,"Set('leftMargin', '2.09cm')\n");
	fprintf(vsz,"Set('rightMargin', '0.169cm')\n");
	fprintf(vsz,"Set('topMargin', u'0.276pt')\n");
	fprintf(vsz,"Set('bottomMargin', '1.71cm')\n");
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
	fprintf(vsz,"Add('label', name='label1', autoadd=False)\n");
	fprintf(vsz,"To('label1')\n");
	fprintf(vsz,"Set('label', u'M = %1.2f M_{\\odot}')\n",m2);
	fprintf(vsz,"Set('Text/size', u'18pt')\n");
	fprintf(vsz,"Set('xPos', [0.039439388406118166])\n");
	fprintf(vsz,"Set('yPos', [0.90394296957670195])\n");
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
	
	fprintf(vsz,"Add('graph', name='graph4', autoadd=False)\n");
	fprintf(vsz,"To('graph4')\n");
	fprintf(vsz,"Set('leftMargin', '3.53cm')\n");
	fprintf(vsz,"Set('rightMargin', '0.1cm')\n");
	fprintf(vsz,"Set('topMargin', u'0')\n");
	fprintf(vsz,"Set('bottomMargin', u'1.71cm')\n");
	fprintf(vsz,"Add('axis', name='axis1', autoadd=False)\n");
	fprintf(vsz,"To('axis1')\n");
	fprintf(vsz,"Set('label', u'x [cm]')\n");
	fprintf(vsz,"Set('min', -1000000000.0)\n");
	fprintf(vsz,"Set('max', 1000000000.0)\n");
	fprintf(vsz,"Set('reflect', False)\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"Set('MajorTicks/hide', False)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='axis2', autoadd=False)\n");
	fprintf(vsz,"To('axis2')\n");
	fprintf(vsz,"Set('label', u'y [cm]')\n");
	fprintf(vsz,"Set('min', -1000000000.0)\n");
	fprintf(vsz,"Set('max', 1000000000.0)\n");
	fprintf(vsz,"Set('reflect', False)\n");
	fprintf(vsz,"Set('direction', u'vertical')\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('Label/atEdge', False)\n");
	fprintf(vsz,"Set('Label/rotate', False)\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"Set('TickLabels/format', u'Auto')\n");
	fprintf(vsz,"Set('MajorTicks/hide', False)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
	fprintf(vsz,"To('x')\n");
	fprintf(vsz,"Set('TickLabels/hide', True)\n");
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
	fprintf(vsz,"Set('label', u't = %3.5f')\n",tpos);
	fprintf(vsz,"Set('Text/size', u'18pt')\n");
	fprintf(vsz,"Set('Text/color', u'black')\n");
	fprintf(vsz,"Set('xPos', [-0.3748122380805215])\n");
	fprintf(vsz,"Set('yPos', [1.9857883369206168])\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('colorbar', name='colorbar1', autoadd=False)\n");
	fprintf(vsz,"To('colorbar1')\n");
	fprintf(vsz,"Set('image', u'image1')\n");
	fprintf(vsz,"Set('label', u'rho [g/cc]')\n");
	fprintf(vsz,"Set('autoExtend', False)\n");
	fprintf(vsz,"Set('TickLabels/color', u'black')\n");
	fprintf(vsz,"Set('vertPosn', u'top')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('image', name='image1', autoadd=False)\n");
	fprintf(vsz,"To('image1')\n");
	fprintf(vsz,"Set('data', u'ds')\n");
	fprintf(vsz,"Set('min', 1.0)\n");
	fprintf(vsz,"Set('max', 30000000.0)\n");
	fprintf(vsz,"Set('colorScaling', u'log')\n");
	fprintf(vsz,"Set('colorMap', u'bartlein')\n");
	fprintf(vsz,"Set('colorInvert', True)\n");
	fprintf(vsz,"Set('transparency', 0)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	
	
	
	
	fclose(vsz);
	free(body);
	
	return 0;
}