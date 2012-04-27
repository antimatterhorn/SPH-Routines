/*
 *  sdf2csv.h
 *  loads a SDF file and outputs a CSV file suitable for plotting
 *	in an application like Veusz
 *
 *  Created by Cody Raskin on 7/14/09.
 *
 */

typedef struct {
	double x, y, z;		/* position of body */
	float mass;	   /* mass of body */
	float vx, vy, vz;     /* velocity of body */
	float u;      	   /* specific energy of body*/
	float h;      	   /* smoothing length of body */
	float rho;            /* density of body */
	float pr;            /* pressure of body */
	float drho_dt;        /* drho/dt of body */
	float udot;           /* du/dt of body */
	float temp;           /* temperature of body */
	float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */
	float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */
	float vsound;		/* sound speed */
	float abar;           /* avg number of nucleons per particle of body */
	float zbar;           /* avg number of protons per particle of body */
	float ax, ay, az;     /* acceleration of body */
	float lax, lay, laz;  /* last acceleration of body */
	float phi;            /* potential at body location */
	float dt;           /* timestep */
	float bdot;			/* energy generation rate */
	float kappa;			/* opacity */
	float sigma;			/* conductivity */
	unsigned int openup;	/* exposed upward */
	unsigned int opendown;	/* exposed downward */
	unsigned int openout;	/* exposed outward */
	float durad;			/* du radiated away */
	unsigned int nbrs;     /* number of neighbors */
	unsigned int ident;	   /* unique identifier */
	unsigned int windid;   /* wind id */
	unsigned int useless;	/* to fill double block */	
} SPHbody;


const double pi = 3.1415926;
const double onethird = 1.0 / 3.0;
const double twothirds = 2.0 / 3.0;
const double fourthirds = 4.0 / 3.0;
