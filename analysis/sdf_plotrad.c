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
float x,rho,temp,mass,pressure,h,cs,vx,tpos,u,prog;
double pe;
int i,j,k,l;
double g = 6.67e-8;
double mass_in_g = 1.989E+27;        
double dist_in_cm = 6.955e7;     
double time_in_s = 1;
double energy_in_erg,dens_in_gccm,pressure_in_ergperccm,specenergy_in_ergperg;

char sdffile[80];
char csvfile[80];
char vszfile[80];
char pdffile[80];

void quickSort(double **arr, double temp[], int nel, int left, int right)
{
	int m = left, n = right;
	
	double pivot = arr[(left + right) / 2][0];
	
	printf("pivot = %3.2f\n",pivot);
	
	/* partition */
	
	while (m <= n) {
		
		while (arr[m][0] < pivot) m++;
		while (arr[n][0] > pivot) n--;
		if (m <= n) {
			
			for (j=0;j<nel;j++) {
				temp[j] = arr[m][j];
				arr[m][j] = arr[n][j];
				arr[n][j] = temp[j];
			}

			m++;
			n--;
		}
	};

	
	/* recursion */
	
	if (left < n) quickSort(arr, temp, nel, left, n);
	if (m < right) quickSort(arr, temp, nel, m, right);
}

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
					NULL);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);

	int na = 15;
	int dv = 20;
	int np = nobj/dv;
	double **inArray;
	double **outArray;
	double tempArray[na];
	
	
	inArray = (double **) malloc(np*sizeof(double *));
	for(i=0;i<np;i++){
		inArray[i]  = (double *) malloc(na*sizeof(double));
	}
	
	for(i=0;i<np;i++){
		for(j=0;j<na;j++){
			inArray[i][j]  = 0.0;
		}
		
	}
	

	
	j=l=0;
	for(i = 0; i < nobj; i++){
		if (i % dv == 0 && j < np){
			rho = body[i].rho * dens_in_gccm;
			mass = body[i].mass * mass_in_g; 
			pressure = body[i].pr * pressure_in_ergperccm;
			cs = body[i].vsound * dist_in_cm;
			u = body[i].u * specenergy_in_ergperg;
			inArray[j][0] = sqrt(pow(body[i].x,2.0)+pow(body[i].y,2.0)+pow(body[i].z,2.0))*dist_in_cm;
			inArray[j][3] = mass*dv;  //multiplying by dv for missing mass
			inArray[j][14] = body[i].h - fabs(body[i].z);
			
			if(inArray[j][14] > 0){
				inArray[j][1] = rho;
				inArray[j][2] = body[i].temp;
				inArray[j][4] = pressure;
				inArray[j][5] = cs;
				inArray[j][6] = u;
				inArray[j][7] = 0.5*mass*(pow(body[i].vx,2.0)+pow(body[i].vy,2.0)+pow(body[i].vz,2.0))*pow(dist_in_cm,2.0);
				inArray[j][10] = (body[i].y*body[i].vz - body[i].z*body[i].vy) / (pow(body[i].z,2.0)+pow(body[i].y,2.0));
				inArray[j][11] = (body[i].z*body[i].vx - body[i].x*body[i].vz) / (pow(body[i].x,2.0)+pow(body[i].z,2.0));
				inArray[j][12] = (body[i].x*body[i].vy - body[i].y*body[i].vx) / (pow(body[i].x,2.0)+pow(body[i].y,2.0));
				l++;
			}

			//if(body[i].h - fabs(body[i].z) < 0) inArray[j][12] = 0;
			
			j++;
		}

	}
	
	printf("sorting...\n");
	
	//	quickSort(inArray,tempArray,na,0,nobj);
	
	for (i=0; i< (np -1); i++)    // element to be compared
    {
		prog = (double)i/(double)np*100.0;
		if (i % 100 == 0.0) printf("progress %3.2f\n",prog);
		for(j = (i+1); j < np; j++)   // rest of the elements
		{
			if (inArray[i][0] > inArray[j][0])          // descending order
			{
				for (k=0;k<na;k++) 
				{
					tempArray[k] = inArray[i][k];
					inArray[i][k] = inArray[j][k];
					inArray[j][k] = tempArray[k];
				}
			}
		}
	}
	
	mass = 0;	
	for (i=0;i<np;i++)
	{
		mass	+= inArray[i][3];
		pe		= g*mass*mass/inArray[i][0];
		inArray[i][8] = pe;
		inArray[i][9] = mass;
		inArray[i][13] = sqrt(mass*6.67e-8/pow(inArray[i][0],3.0));
		//if(inArray[i][14] < 0) inArray[i][12] = 0;
	}
	
	outArray = (double **) malloc(l*sizeof(double *));
	for(i=0;i<l;i++){
		outArray[i]  = (double *) malloc(na*sizeof(double));
	}
	
	k=0;
	for(i=0;i<np;i++){
		if(inArray[i][1]>0.0 && k<l){
			for(j=0;j<na;j++){
				outArray[k][j]  = inArray[i][j];
				
			}
			k++;
		}

		
	}
	
	
	//open the stream file
	snprintf(csvfile, sizeof(csvfile), "%s_rad.csv", argv[1]);
	FILE *stream, *fopen(), *vsz;
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");

	fprintf(stream,"r,rho,temp,mass,pr,cs,u,ke,pe,mu,ox,oy,oz,ok\n");
	
	for(i = 0; i < l; i++)
	{
		fprintf(stream,"%3.3e,%3.3e,%3.3e,%3.3e,%3.3e,%3.3e,%3.3e,%3.3e,%3.3e,%3.3e,%3.3e,%3.3e,%3.3e,%3.3e\n", 
										outArray[i][0],outArray[i][1],outArray[i][2],
										outArray[i][3],outArray[i][4],outArray[i][5],outArray[i][6],
										outArray[i][7],outArray[i][8],outArray[i][9],outArray[i][10],outArray[i][11],outArray[i][12],
										outArray[i][13]);
	}
	fclose(stream);
	
	char *cwd = getcwd(NULL, 0);
	snprintf(vszfile, sizeof(vszfile), "%s_rad.vsz", argv[1]);
	snprintf(pdffile, sizeof(pdffile), "%s_rad.png", argv[1]);
	
	vsz = fopen(vszfile,"w");
	
	fprintf(vsz,"AddImportPath(u'%s')\n",cwd);
	fprintf(vsz,"ImportFileCSV(u'%s', linked=True)\n",csvfile);
	fprintf(vsz,"Add('page', name='page1', autoadd=False)\n");
	fprintf(vsz,"To('page1')\n");
	fprintf(vsz,"Add('graph', name='graph1', autoadd=False)\n");
	fprintf(vsz,"To('graph1')\n");
	fprintf(vsz,"Set('leftMargin', '2.09cm')\n");
	fprintf(vsz,"Set('rightMargin', '0.169cm')\n");
	fprintf(vsz,"Set('topMargin', '0.276cm')\n");
	fprintf(vsz,"Set('bottomMargin', '1.71cm')\n");
	fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
	fprintf(vsz,"To('x')\n");
	fprintf(vsz,"Set('label', u'Mass [M_{\\odot}]')\n");
	fprintf(vsz,"Set('min', 0.01)\n");
	fprintf(vsz,"Set('max', 'Auto')\n");
	fprintf(vsz,"Set('log', True)\n");
	fprintf(vsz,"Set('autoExtend', False)\n");
	fprintf(vsz,"Set('datascale', 5.0000000000000003e-34)\n");
	fprintf(vsz,"Set('lowerPosition', 0.0)\n");
	fprintf(vsz,"Set('upperPosition', 1.0)\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
	fprintf(vsz,"To('y')\n");
	fprintf(vsz,"Set('label', u'\\epsilon [erg]')\n");
	fprintf(vsz,"Set('min', 'Auto')\n");
	fprintf(vsz,"Set('max', 'Auto')\n");
	fprintf(vsz,"Set('log', True)\n");
	fprintf(vsz,"Set('direction', 'vertical')\n");
	fprintf(vsz,"Set('Label/size', u'20pt')\n");
	fprintf(vsz,"Set('Label/color', u'black')\n");
	fprintf(vsz,"Set('TickLabels/size', u'18pt')\n");
	fprintf(vsz,"Set('TickLabels/color', u'black')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('key', name='key1', autoadd=False)\n");
	fprintf(vsz,"To('key1')\n");
	fprintf(vsz,"Set('horzPosn', u'left')\n");
	fprintf(vsz,"Set('vertPosn', u'top')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy1', autoadd=False)\n");
	fprintf(vsz,"To('xy1')\n");
	fprintf(vsz,"Set('xData', u'mu')\n");
	fprintf(vsz,"Set('yData', u'ke')\n");
	fprintf(vsz,"Set('markerSize', u'1pt')\n");
	fprintf(vsz,"Set('key', u'KE')\n");
	fprintf(vsz,"Set('PlotLine/hide', True)\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"Add('xy', name='xy2', autoadd=False)\n");
	fprintf(vsz,"To('xy2')\n");
	fprintf(vsz,"Set('xData', u'mu')\n");
	fprintf(vsz,"Set('yData', u'pe')\n");
	fprintf(vsz,"Set('marker', u'none')\n");
	fprintf(vsz,"Set('markerSize', u'0.75pt')\n");
	fprintf(vsz,"Set('key', u'PE')\n");
	fprintf(vsz,"Set('yAxis', u'y')\n");
	fprintf(vsz,"Set('PlotLine/hide', False)\n");
	fprintf(vsz,"Set('MarkerLine/color', u'green')\n");
	fprintf(vsz,"Set('MarkerFill/color', u'green')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	fprintf(vsz,"To('..')\n");
	

	fclose(vsz);
	
	return 0;
}