#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
#include "Msgs.h"
#include <stdlib.h>


/* This code will setup the initial conditions SDF file		*/
/* for a two-body system, given positions and velocities	*/
/* of each.													*/

#define SPHOUTBODYDESC \
"struct {\n\
	double x, y, z;		/* position of body */\n\
	float mass;			/* mass of body */\n\
	float vx, vy, vz;		/* velocity of body */\n\
	float u;			/* internal energy */\n\
	float h;			/* smoothing length */\n\
	float rho;			/* density */\n\
	float pr;			/* pressure */\n\
	float drho_dt;              /* time derivative of rho */\n\
	float udot;			/* time derivative of u */\n\
	float temp;           /* temperature of body */\n\
	float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */\n\
	float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */\n\
	float vsound;		/* sound speed */\n\
	float abar;           /* avg number of nucleons per particle of body */\n\
	float zbar;           /* avg number of protons per particle of body */\n\
	float ax, ay, az;		/* acceleration */\n\
	float lax, lay, laz;	/* acceleration at tpos-dt */\n\
	float phi;			/* potential */\n\
	float idt;			/* timestep */\n\
	unsigned int nbrs;          /* number of neighbors */\n\
	unsigned int ident;		/* unique identifier */\n\
	unsigned int windid;        /* wind id */\n\
}"


typedef struct {
	double x, y, z;             /* position of body */
	float mass;           /* mass of body */
	float vx, vy, vz;     /* velocity of body */
	float u;              /* specific energy of body*/
	float h;              /* smoothing length of body */
	float rho;            /* density of body */
	float pr;            /* pressure of body */
	float drho_dt;        /* drho/dt of body */
	float udot;           /* du/dt of body */
	float temp;           /* temperature of body */
	float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */
	float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */
	float vsound;			/* sound speed */
	float abar;           /* avg number of nucleons per particle of body */
	float zbar;           /* avg number of protons per particle of body */
	float ax, ay, az;     /* acceleration of body */
	float lax, lay, laz;  /* last acceleration of body */
	//float gax, gay, gaz;  /* gravity acceleration of body */
	//float grav_mass;      /* gravitational mass of body */
	float phi;            /* potential at body location */
	//float tacc;           /* time of last acceleration update of body */
	float idt;
	unsigned int nbrs;     /* number of neighbors */
	unsigned int ident;    /* unique identifier */
	unsigned int windid;   /* wind id */
	//unsigned int useless;  /* to fill the double block (4int=1double)*/		
} SPHbody;

int iter = 0;
float dt = 0;
float eps = 0;
float Gnewt = 0;
float tolerance = 0;
float frac_tolerance = 0;
int ndim = 3;
float tpos = 0;
float tvel = 0;
float gamma = 0;
double ke = 0;
double pe = 0;
double te = 0;
double tmass;

float r,rmax,rowtot;
float rf,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12;
int lines;

int gnobj1, nobj1;
int conf,i,j,k;

FILE *file1;
SDF *sdfp;
SPHbody *body1;
double **compos;

char outfile[80];

int usage()
{
	printf("\t Interpolates on a composition file of the form r/rmax i0 i1 i2 i3...i12\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t sdf_composition [sdf file] [composition file] [lines]\n");
	return 0;
}

int main(int argc, char **argv[])
{
	


	
	if (argc < 4){
		usage();
		return 0;
	}
	else {
		lines = atoi(argv[3]);
	}
	
	
	
	
	
	compos = (double **) malloc(lines*sizeof(double *));
	for(i=0;i<lines;i++){
		compos[i]  = (double *) malloc(14*sizeof(double));
	}

	sdfp = SDFreadf(argv[1], (void **)&body1, &gnobj1, &nobj1, sizeof(SPHbody),
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
					"vsound", offsetof(SPHbody, vsound), &conf,
					"abar", offsetof(SPHbody, abar), &conf,
					"zbar", offsetof(SPHbody, zbar), &conf,
					"ax", offsetof(SPHbody, ax), &conf,
					"ay", offsetof(SPHbody, ay), &conf,
					"az", offsetof(SPHbody, az), &conf,
					"lax", offsetof(SPHbody, lax), &conf,
					"lay", offsetof(SPHbody, lay), &conf,
					"laz", offsetof(SPHbody, laz), &conf,
					//"gax", offsetof(SPHbody, gax), &conf,
					//"gay", offsetof(SPHbody, gay), &conf,
					//"gaz", offsetof(SPHbody, gaz), &conf,
					//"grav_mass", offsetof(SPHbody, grav_mass), &conf,
					"phi", offsetof(SPHbody, phi), &conf,
					//"tacc", offsetof(SPHbody, tacc), &conf,
					"idt", offsetof(SPHbody, idt), &conf,
					"nbrs", offsetof(SPHbody, nbrs), &conf,
					"ident", offsetof(SPHbody, ident), &conf,
					"windid", offsetof(SPHbody, windid), &conf,
					//"useless", offsetof(SPHbody, useless), &conf,
					NULL);
	//now to fill the misc. floats and ints//
	SDFgetfloatOrDefault(sdfp, "dt",  &dt, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "eps",  &eps, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "Gnewt",  &Gnewt, (float)0.000393935);
	SDFgetfloatOrDefault(sdfp, "tolerance",  &tolerance, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "frac_tolerance",  &frac_tolerance, (float)0.0);
	//SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "tvel",  &tvel, (float)0.0);
	SDFgetfloatOrDefault(sdfp, "gamma",  &gamma, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "ke",  &ke, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "pe",  &pe, (float)0.0);
	SDFgetdoubleOrDefault(sdfp, "te",  &te, (float)0.0);
	//SDFgetintOrDefault  (sdfp, "iter",  &iter, 0);
	
	printf("Opening %s\n",argv[2]);
	
	file1 = fopen(argv[2],"r");
	i=0;
	while (fscanf(file1,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f",&rf,&i0,&i1,&i2,&i3,&i4,&i5,&i6,&i7,&i8,&i9,&i10,&i11,&i12) !=EOF)
	{
		printf("reading %d:\t%3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f\n",
			   i,
			   rf,
			   i0,
			   i1,
			   i2,
			   i3,
			   i4,
			   i5,
			   i6,
			   i7,
			   i8,
			   i9,
			   i10,
			   i11,
			   i12);
		
		compos[i][0] = rf;
		compos[i][1] = i0;
		compos[i][2] = i1;
		compos[i][3] = i2;
		compos[i][4] = i3;
		compos[i][5] = i4;
		compos[i][6] = i5;
		compos[i][7] = i6;
		compos[i][8] = i7;
		compos[i][9] = i8;
		compos[i][10] = i9;
		compos[i][11] = i10;
		compos[i][12] = i11;
		compos[i][13] = i12;
		i++;
	}
	
	fclose(file1);
	
	printf("composition:\n");
	printf("\t\tR/Rm He4  C12  O16  Ne20 Mg24 Si28 S32  Ar36 Ca40 Ti44 Cr48 Fe52 Ni56\n");
	
	for(i=0;i<lines;i++)
	{
		rowtot = 0;
		for(j=0;j<13;j++)
		{
			rowtot += compos[i][j+1];
		}
		for(j=0;j<13;j++)
		{
			compos[i][j+1] = compos[i][j+1]/rowtot;
		}
		printf("\t\t%3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f\n",
			   compos[i][0],
			   compos[i][1],
			   compos[i][2],
			   compos[i][3],
			   compos[i][4],
			   compos[i][5],
			   compos[i][6],
			   compos[i][7],
			   compos[i][8],
			   compos[i][9],
			   compos[i][10],
			   compos[i][11],
			   compos[i][12],
			   compos[i][13]);		
	}
	
	rmax = 0;
	for(j=0;j<nobj1;j++)
	{
		r = sqrt(pow(body1[j].x,2.0)+pow(body1[j].y,2.0)+pow(body1[j].z,2.0));
		if (r>rmax) rmax = r;
	}
	
	for(j=0;j<nobj1;j++)
	{
		r = sqrt(pow(body1[j].x,2.0)+pow(body1[j].y,2.0)+pow(body1[j].z,2.0));
		rf = r/rmax;
		for(i=0;i<lines;i++)
		{
			if(compos[i][0] < r/rmax)
			{
				k = i+1;
				body1[j].He4 = (rf-compos[i][0])*(compos[k][1]-compos[i][1])/(compos[k][0]-compos[i][0])+compos[i][1];
				body1[j].C12 = (rf-compos[i][0])*(compos[k][2]-compos[i][2])/(compos[k][0]-compos[i][0])+compos[i][2];
				body1[j].O16 = (rf-compos[i][0])*(compos[k][3]-compos[i][3])/(compos[k][0]-compos[i][0])+compos[i][3];
				body1[j].Ne20 = (rf-compos[i][0])*(compos[k][4]-compos[i][4])/(compos[k][0]-compos[i][0])+compos[i][4];
				body1[j].Mg24 = (rf-compos[i][0])*(compos[k][5]-compos[i][5])/(compos[k][0]-compos[i][0])+compos[i][5];
				body1[j].Si28 = (rf-compos[i][0])*(compos[k][6]-compos[i][6])/(compos[k][0]-compos[i][0])+compos[i][6];
				body1[j].S32 = (rf-compos[i][0])*(compos[k][7]-compos[i][7])/(compos[k][0]-compos[i][0])+compos[i][7];
				body1[j].Ar36 = (rf-compos[i][0])*(compos[k][8]-compos[i][8])/(compos[k][0]-compos[i][0])+compos[i][8];
				body1[j].Ca40 = (rf-compos[i][0])*(compos[k][9]-compos[i][9])/(compos[k][0]-compos[i][0])+compos[i][9];
				body1[j].Ti44 = (rf-compos[i][0])*(compos[k][10]-compos[i][10])/(compos[k][0]-compos[i][0])+compos[i][10];
				body1[j].Cr48 = (rf-compos[i][0])*(compos[k][11]-compos[i][11])/(compos[k][0]-compos[i][0])+compos[i][11];
				body1[j].Fe52 = (rf-compos[i][0])*(compos[k][12]-compos[i][12])/(compos[k][0]-compos[i][0])+compos[i][12];
				body1[j].Ni56 = (rf-compos[i][0])*(compos[k][13]-compos[i][13])/(compos[k][0]-compos[i][0])+compos[i][13];
			}
		}
	}
	
		
	//write the output file
	snprintf(outfile, sizeof(outfile), "compos_%s", argv[1]);
	//*fp = fopen(outfile, "w");
	
	SDFwrite(outfile, gnobj1, 
			 nobj1, body1, sizeof(SPHbody),
			 SPHOUTBODYDESC,
			 "npart", SDF_INT, gnobj1,
			 "iter", SDF_INT, iter,
			 "dt", SDF_FLOAT, dt,
			 "eps", SDF_FLOAT, eps,
			 "Gnewt", SDF_FLOAT, Gnewt,
			 "tolerance", SDF_FLOAT, tolerance,
			 "frac_tolerance", SDF_FLOAT, frac_tolerance,
			 "ndim", SDF_INT, ndim,
			 "tpos", SDF_FLOAT, tpos,
			 "tvel", SDF_FLOAT, tvel,
			 //"t_wind", SDF_FLOAT, twind_out,
			 //"R0", SDF_FLOAT, output_R0,
			 //"Omega0", SDF_FLOAT, cosmo.Omega0,
			 //"H0", SDF_FLOAT, cosmo.H0,
			 //"Lambda_prime", SDF_FLOAT, cosmo.Lambda,
			 //"hubble", SDF_FLOAT, output_h,
			 //"redshift", SDF_FLOAT, output_z,
			 "gamma", SDF_FLOAT, gamma,
			 //"centmass", SDF_FLOAT, centmass, 
			 "ke", SDF_DOUBLE, ke,
			 "pe", SDF_DOUBLE, pe,
			 "te", SDF_DOUBLE, te,
			 NULL);
	free(body1);
    free(compos);
	
	return 0;
}


