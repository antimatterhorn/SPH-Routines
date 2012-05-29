/*
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

char exec[80];
char paramfile[80];
char params[80];
char root[80];
char search[80];
char period[1];
char space[80];
char command[80];
char number[6];
int i,j,k,start,end,inc,lines,pixels,choice;

int fi,fc;
float fx,fr,ft;

double **par;

char vszfile[80];

int main()
{
	

    
    
    
    printf("exec: ");	
	gets (exec);
	printf("on root: ");
	gets (root);
	printf("param file of the form (iter xmax rhomax tempmax): ");
	gets (paramfile);
    printf("lines in param file: ");
    scanf("%d",&lines);
    printf("pixels: ");
    scanf("%d",&pixels);
	printf("start end increment: ");
	scanf("%d %d %d",&start,&end,&inc);
    
    par = (double **) malloc(lines*sizeof(double *));
	for(i=0;i<lines;i++){
		par[i]  = (double *) malloc(4*sizeof(double));
	}
	
    FILE *file1;
    file1 = fopen(paramfile,"r");
	i=0;
	while (fscanf(file1,"%d %d %f %f %f",&fc,&fi,&fx,&fr,&ft) !=EOF)
	{
		printf("reading %d:%d:\t%3.2e %3.2e %3.2e\n",
			   fc,
               fi,
			   fx,
			   fr,
			   ft);
		
		par[i][0] = fi;
		par[i][1] = fx;
		par[i][2] = fr;
		par[i][3] = ft;
		i++;
	}	
	fclose(file1);

    
	for(i=start;i<(end+1);i+=inc)
	{
		for(j=0;j<lines;j++)
        {
            if(par[j][0] <= i)
            {
                k   = j+1;
                fx  = (i-par[j][0])*(par[k][1]-par[j][1])/(par[k][0]-par[j][0])+par[j][1];
                fr  = (i-par[j][0])*(par[k][2]-par[j][2])/(par[k][0]-par[j][0])+par[j][2];
                ft  = (i-par[j][0])*(par[k][3]-par[j][3])/(par[k][0]-par[j][0])+par[j][3];
                
                snprintf(params,sizeof(params),"%d %d %3.2f %3.2e %3.2e",fc,pixels,fx,fr,ft);
            }
        }
        
        
        snprintf(command, sizeof(command),"%s %s.%04d %s",exec,root,i,params);
		printf(":> %s\n",command);
        
		//if(i<1000){
//			snprintf(command, sizeof(command),"%s %s.0%d %s",exec,root,i,params);
//		}else {
//			snprintf(command, sizeof(command),"%s %s.%d %s",exec,root,i,params);
//		}
		system(command);
	}
    
    system("rm -r png/");
	system("mkdir png");
	
	char *cwd = getcwd(NULL, 0);
	
	snprintf(vszfile, sizeof(vszfile), "exportImages.vsz");
	
	FILE *vsz;
	vsz = fopen(vszfile,"w");
	
	fprintf(vsz,"import os.path\n");
	fprintf(vsz,"os.getcwd()\n");
	fprintf(vsz,"import glob\n");
	fprintf(vsz,"examples = glob.glob(os.path.join('%s', '%s*.vsz'))\n",cwd,root);
	fprintf(vsz,"for filename in examples:\n");
	fprintf(vsz,"   # this will delete all the widgets ready for the next file\n");
	fprintf(vsz,"   To('/')\n");
	fprintf(vsz,"   for child in GetChildren():\n");
	fprintf(vsz,"       Remove(child)\n");
	fprintf(vsz,"\n");
	fprintf(vsz,"   execfile(filename)\n");
	fprintf(vsz,"\n");
	fprintf(vsz,"   Export('%s/png/%s.png' %s os.path.basename(filename))",cwd,"%s","%");
	
	fclose(vsz);
	

	
	return 0;
}