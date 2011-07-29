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
char params[80];
char root[80];
char search[80];
char period[1];
char space[80];
char command[80];
char number[6];
int i,start,end,inc;

char vszfile[80];

int main()
{
	strcpy(space," ");
	strcat(exec,"sdf_opacity ");
	printf("on file: ");
	gets (root);
	//strcat(space,params);
	//strcpy(params,space);
	printf("start end increment: ");
	scanf("%d %d %d",&start,&end,&inc);
	
	for(i=start;i<(end+1);i+=inc)
	{
		strcpy(command,exec);
		strcat(command,root);
		strcat(command,space);
		strcat(command,"200 20 5e6 ");
		sprintf(number,"%d",i);
		//strcat(command,space);
		strcat(command,number);
		strcat(command," 0 0.8");
		printf(command);
		system(command);
	}
	
	char *cwd = getcwd(NULL, 0);
	
	snprintf(vszfile, sizeof(vszfile), "exportImages.vsz");
	
	FILE *vsz;
	vsz = fopen(vszfile,"w");
	
	fprintf(vsz,"import os.path\n");
	fprintf(vsz,"os.getcwd()\n");
	fprintf(vsz,"import glob\n");
	fprintf(vsz,"examples = glob.glob(os.path.join('%s', 'fltmass_sph*.vsz'))\n",cwd);
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