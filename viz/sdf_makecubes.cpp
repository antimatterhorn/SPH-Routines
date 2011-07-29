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


char root[80];
char search[80];
char period[1];
char command[80];
char number[6];
int i,start,end,inc;

int main()
{
		
	printf("root: ");
	gets (root);
	strcpy(period,".");
	strcat(root,period);
	printf("start end increment: ");
	scanf("%d %d %d",&start,&end,&inc);
	
	for(i=start;i<(end+1);i+=inc)
	{
		strcpy(command,"sdf_cube ");
		strcat(command,root);
		sprintf(number,"%d",i);
		strcat(command,number);
		strcat(command," 100 0.04");
		system(command);
	}
	
	return 0;
}