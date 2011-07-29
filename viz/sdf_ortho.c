/*
 *  sdf2csv.c
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  cxyz is the location of the camera
 *	thetaxyz is the rotation of the camera
 *	exyz is the viewers position relative to the display
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
int conf;
float xmax,ymax,rr;
float rho,temp,mass,h,tpos;
int choice,pixels;
int i,j,imi,imj,ic,jc,r;
double pos[3],d[3],b[2],c[3],e[3],theta[3],x,y,z;
char sdffile[80];
char csvfile[80];
char vszfile[80];
char pdffile[80];

double cos( double arg );
double sin( double arg );

void x2i()
{
	ic = x*(pixels/(2.0*xmax)) + pixels/2.0;	
}

void y2j()
{
	jc = -y*(pixels/(2.0*ymax)) + pixels/2.0;
}

void dmat()
{
	d[0] = cos(theta[1])*(sin(theta[2])*(pos[1]-c[1]) + cos(theta[2])*(pos[0]-c[0]))-sin(theta[1])*(pos[2]-c[2]);
	d[1] = sin(theta[0])*(cos(theta[1])*(pos[2]-c[2]) + sin(theta[1])*(sin(theta[2])*(pos[1]-c[1]) + cos(theta[2])*(pos[0]-c[0]))) + cos(theta[0])*(cos(theta[2])*(pos[1]-c[1]) - sin(theta[2])*(pos[0]-c[0]));
	d[2] = cos(theta[0])*(cos(theta[1])*(pos[2]-c[2]) + sin(theta[1])*(sin(theta[2])*(pos[1]-c[1]) + cos(theta[2])*(pos[0]-c[0]))) - sin(theta[0])*(cos(theta[2])*(pos[1]-c[1]) - sin(theta[2])*(pos[0]-c[0]));
}

void bmat()
{
	b[0] = (d[0]-e[0])*(e[2]/d[2]);
	b[1] = (d[1]-e[1])*(e[2]/d[2]);
}

int main(int argc, char **argv[])
{
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
	
	singlPrintf("%s has %d particles.\n", argv[1], gnobj);
	
	if (argc < 3){
		printf("\nOutput Parameters \n-----------------\nChoose a Variable:\n");
		printf("1) Density \n2) Temperature\n\nChoice:");
		scanf("%d", &choice);
	}
	else {
		choice = atoi(argv[2]);
	}
	
	if (argc < 4){
		printf("Pixels on a side:");
		scanf("%d", &pixels);
	}
	else {
		pixels = atoi(argv[3]);
	}
	
	if (argc < 5){
		printf("xmax:");
		scanf("%f", &xmax);
	}
	else {
		xmax = atof(argv[4]);
	}
	
	ymax = xmax;
	
	c[0]=-xmax*2.0;
	c[2]=xmax;
	//c[2] = -c[2]; //move to back corner
	e[0]=e[1]=0;
	e[2]=xmax;
	theta[1] = 3.14159*.75;
	
	
	double **dens;
	
	dens = (double **) malloc(pixels*sizeof(double *));
	for(i=0;i<pixels;i++){
		dens[i]  = (double *) malloc(pixels*sizeof(double));
	}
	
	double **img;
	
	img = (double **) malloc(pixels*sizeof(double *));
	for(i=0;i<pixels;i++){
		img[i]  = (double *) malloc(pixels*sizeof(double));
	}
	
	
	SPHbody *p;
	i=0;
	for(p = body; p < body+gnobj; p++)
	{

		
		x = p->x;
		y = p->y;
		z = p->z;
		h = p->h;
		rho = p->rho;
		temp = p->temp;
		//mass = p->mass;
		
		pos[0] = x;
		pos[1] = y;
		pos[2] = z;
		
		//h = sqrt(pow(h,2.0)-pow(z,2.0));
		dmat();
		bmat();
		x = b[0];
		y = b[1];
		x2i();
		y2j();
		//r = h*(pixels/(2.0*xmax));
		
		if (ic >= 0 && ic < pixels && jc >= 0 && jc < pixels)
		{
			dens[ic][jc] += rho;
			switch( choice )
			{
				case 1 :
					img[ic][jc] += p->rho*rho;
					break;
				case 2 :
					img[ic][jc] += p->temp*rho;
					break;
			}
		}
			
		
//		if (r==0)
//		{
//			//printf("found that r was too small, filling only 1 pixel\n");
//			if (ic >= 0 && ic < pixels && jc >= 0 && jc < pixels)
//			{
//				dens[ic][jc] += rho;
//				switch( choice )
//				{
//					case 1 :
//						img[ic][jc] += p->rho*rho;
//						break;
//					case 2 :
//						img[ic][jc] += p->temp*rho;
//						break;
//				}
//			}
//			
//		}
//		else
//		{
//			for(imi = ic-2*r; imi < ic+2*r+1; imi++)
//			{
//				for(imj = jc-2*r; imj < jc+2*r+1; imj++)
//				{
//					if (imi >= 0 && imi < pixels && imj >= 0 && imj < pixels) // inside the image
//					{
//						rr = sqrt(pow((imi-ic),2.0)+pow((imj-jc),2.0));
//						//cout << rr << "," << r << endl;
//						if (rr <= 2*r) // inside the circle
//						{							
//							dens[imi][imj] += rho;		
//							switch( choice )
//							{
//								case 1 :
//									img[imi][imj] += p->rho*rho*w(h,2*rr*(double)xmax/(double)pixels)*pi*pow(h,3.0);
//									break;
//								case 2 :
//									img[imi][imj] += p->temp*rho*w(h,2*rr*(double)xmax/(double)pixels)*pi*pow(h,3.0);
//									break;
//							}
//						}
//					}
//					
//				}
//			}	
//		}

	}
	
	//open the stream file
	snprintf(csvfile, sizeof(csvfile), "%s.csv", argv[1]);
	FILE *stream, *fopen(), *vsz;
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(csvfile,"w");
	
//	double maxv = 0.0;
//	double minv = pow(10.0,10.0);
	for(i = 0; i < pixels; i++)
	{
		for(j = 0; j < pixels; j++)
		{

			fprintf(stream,"%f ", (double)img[i][j]/(double)dens[i][j]);
		
		}
		fprintf(stream,"\n");
	}
	
	//close the stream file
	fclose(stream);
	
	char *cwd = getcwd(NULL, 0);

	
	switch (choice) {
		case 1:
			snprintf(vszfile, sizeof(vszfile), "%s_rho.vsz", argv[1]);
			snprintf(pdffile, sizeof(pdffile), "%s_rho.png", argv[1]);
			
			vsz = fopen(vszfile,"w");
			
			fprintf(vsz,"ImportFile2D(u'%s/%s', [u'ds'], invertrows=False, invertcols=False, transpose=True, linked=True)\n",cwd,csvfile);
			fprintf(vsz,"Add('page', name='page1', autoadd=False)\n");
			fprintf(vsz,"To('page1')\n");
			fprintf(vsz,"Add('graph', name='graph1', autoadd=False)\n");
			fprintf(vsz,"To('graph1')\n");
			fprintf(vsz,"Set('leftMargin', u'0cm')\n");
			fprintf(vsz,"Set('rightMargin', u'0cm')\n");
			fprintf(vsz,"Set('topMargin', u'0cm')\n");
			fprintf(vsz,"Set('bottomMargin', u'0cm')\n");
			fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
			fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
			fprintf(vsz,"To('y')\n");
			fprintf(vsz,"Set('direction', 'vertical')\n");
			fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('label', name='label1', autoadd=False)\n");
			fprintf(vsz,"To('label1')\n");
			fprintf(vsz,"Set('label', u't=%f')\n",tpos);
			fprintf(vsz,"Set('xPos', [0.050000000000000003])\n");
			fprintf(vsz,"Set('yPos', [0.90000000000000002])\n");
			fprintf(vsz,"Set('Text/size', u'18pt')\n");
			fprintf(vsz,"Set('Text/color', u'white')\n");
			fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('colorbar', name='colorbar1', autoadd=False)\n");
			fprintf(vsz,"To('colorbar1')\n");
			fprintf(vsz,"Set('image', u'image1')\n");
			fprintf(vsz,"Set('TickLabels/color', u'white')\n");
			fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('image', name='image1', autoadd=False)\n");
			fprintf(vsz,"To('image1')\n");
			fprintf(vsz,"Set('data', u'ds')\n");
			fprintf(vsz,"Set('min', 1.0)\n");
			fprintf(vsz,"Set('max', 100000.0)\n");
			fprintf(vsz,"Set('colorScaling', u'log')\n");
			fprintf(vsz,"Set('colorMap', u'heat')\n");
			break;
		case 2:
			snprintf(vszfile, sizeof(vszfile), "%s_temp.vsz", argv[1]);
			snprintf(pdffile, sizeof(pdffile), "%s_temp.png", argv[1]);
			
			vsz = fopen(vszfile,"w");
			
			fprintf(vsz,"ImportFile2D(u'%s/%s', [u'ds'], invertrows=False, invertcols=False, transpose=True, linked=True)\n",cwd,csvfile);
			fprintf(vsz,"Add('page', name='page1', autoadd=False)\n");
			fprintf(vsz,"To('page1')\n");
			fprintf(vsz,"Add('graph', name='graph1', autoadd=False)\n");
			fprintf(vsz,"To('graph1')\n");
			fprintf(vsz,"Set('leftMargin', u'0cm')\n");
			fprintf(vsz,"Set('rightMargin', u'0cm')\n");
			fprintf(vsz,"Set('topMargin', u'0cm')\n");
			fprintf(vsz,"Set('bottomMargin', u'0cm')\n");
			fprintf(vsz,"Add('axis', name='x', autoadd=False)\n");
			fprintf(vsz,"Add('axis', name='y', autoadd=False)\n");
			fprintf(vsz,"To('y')\n");
			fprintf(vsz,"Set('direction', 'vertical')\n");
			fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('label', name='label1', autoadd=False)\n");
			fprintf(vsz,"To('label1')\n");
			fprintf(vsz,"Set('label', u't=%f')\n",tpos);
			fprintf(vsz,"Set('xPos', [0.050000000000000003])\n");
			fprintf(vsz,"Set('yPos', [0.90000000000000002])\n");
			fprintf(vsz,"Set('Text/size', u'18pt')\n");
			fprintf(vsz,"Set('Text/color', u'white')\n");
			fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('colorbar', name='colorbar1', autoadd=False)\n");
			fprintf(vsz,"To('colorbar1')\n");
			fprintf(vsz,"Set('image', u'image1')\n");
			fprintf(vsz,"Set('TickLabels/color', u'white')\n");
			fprintf(vsz,"To('..')\n");
			fprintf(vsz,"Add('image', name='image1', autoadd=False)\n");
			fprintf(vsz,"To('image1')\n");
			fprintf(vsz,"Set('data', u'ds')\n");
			fprintf(vsz,"Set('min', 0.0)\n");
			fprintf(vsz,"Set('max', 4000000000.0)\n");
			fprintf(vsz,"Set('colorScaling', u'sqrt')\n");
			fprintf(vsz,"Set('colorMap', u'heat')\n");
			break;
	}

	fclose(vsz);
	
	return 0;
}