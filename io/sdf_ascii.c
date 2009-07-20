#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
//#include "singlio.h"
#include "Msgs.h"
#include <stdlib.h>

typedef struct {
	double x, y, z;             /* position of body */
	float mass;           /* mass of body */
	float vx, vy, vz;     /* velocity of body */
	float u;              /* specific energy of body*/
	float h;              /* smoothing length of body */
	float rho;            /* density of body */
	float pr;            /* pressure of body */
	float temp;           /* temperature of body */	
} SPHbody;

int main()
{
	
	int i,j;

	char name [80];
	char asciifile[80];
	printf("filename: ");
	gets (name);
	singlPrintf("Reading \"%s\"\n", name);

	SDF *sdfp;
	SPHbody *inArray;
	int gnobj, nobj;
	int conf;
	sdfp = SDFreadf(name, (void **)&inArray, &gnobj, &nobj, sizeof(SPHbody),
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
					"temp", offsetof(SPHbody, temp), &conf,
					NULL);

	
	
	//Msgf(("Data read, SPHnobj=%d, SPHgnobj=%d\n", nobj, gnobj));
	//SDFout *outArray = malloc(sizeof(SDFout) * nobj);

	
	double **outArray;
	
	outArray = (double **) malloc(nobj*sizeof(double *));
	for(i=0;i<nobj;i++){
		outArray[i]  = (double *) malloc(4*sizeof(double));
	}
	
	
	for(i = 0; i < nobj; i++){
		outArray[i][0] = inArray[i].x;
		outArray[i][1] = inArray[i].y;
		outArray[i][2] = inArray[i].z;
		outArray[i][3] = inArray[i].rho;
	}
	
	snprintf(asciifile, sizeof(asciifile), "%s.ascii", name);
	FILE *stream, *fopen();
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(asciifile,"w");
	
	for(i = 0; i < nobj; i++)
	{
		fprintf(stream,"%f %f %f %f\n", outArray[i][0],outArray[i][1],outArray[i][2],outArray[i][3]);
	}
	fclose(stream);
	

	return 0;
}