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
	printf("exec: ");	
	gets (exec);
	printf("on root: ");
	gets (root);
	printf("params: ");
	gets (params);
	printf("start end increment: ");
	scanf("%d %d %d",&start,&end,&inc);
	
	for(i=start;i<(end+1);i+=inc)
	{
		snprintf(command, sizeof(command),"%s %s.%04d %s",exec,root,i,params);
		
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