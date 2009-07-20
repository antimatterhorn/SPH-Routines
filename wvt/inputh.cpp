/*
 *  inputh.cpp
 *  
 *
 *  Created by Cody Raskin on 7/15/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "inputh.h"

using namespace std;

const double pi = 3.1415926;
const double onethird = 1.0 / 3.0;
const double twothirds = 2.0 / 3.0;
const double fourthirds = 4.0 / 3.0;

int npart;
char filename[80];
char output[80];

int main()
{
	FILE *file2;
	float d1,d2,d3;
	file2 = fopen("inputmodel.dat","r");
	while (fscanf(file2,"%g %g %*g %*g %*g", &d1,&d2)!=EOF) {
		d3 = pow(d2,-onethird);
		printf("%e %e\n",d1,d3);
	}
		
	fclose(file2);
	return 0;
}


