/*
 *  sdf2csv.h
 *  loads a SDF file and outputs a CSV file suitable for plotting
 *	in an application like Veusz
 *
 *  Created by Cody Raskin on 7/14/09.
 *
 */

typedef struct {
	double x, y, z;             /* position of body */
	double mass;           /* mass of body */
	double vx, vy, vz;     /* velocity of body */
	double u;              /* specific energy of body*/
	double h;              /* smoothing length of body */
	double rho;            /* density of body */
	double pr;            /* pressure of body */
	double drho_dt;        /* drho/dt of body */
	double udot;           /* du/dt of body */
	double temp;           /* temperature of body */
	double He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */
	double Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */
	double abar;           /* avg number of nucleons per particle of body */
	double zbar;           /* avg number of protons per particle of body */
	double ax, ay, az;     /* acceleration of body */
	double lax, lay, laz;  /* last acceleration of body */
	double gax, gay, gaz;  /* gravity acceleration of body */
	double grav_mass;      /* gravitational mass of body */
	double phi;            /* potential at body location */
	double tacc;           /* time of last acceleration update of body */
	double idt;
	unsigned int nbrs;     /* number of neighbors */
	unsigned int ident;    /* unique identifier */
	unsigned int windid;   /* wind id */
	//unsigned int useless;  /* to fill the double block (4int=1double)*/		
} SPHbody;


const double pi = 3.1415926;
const double onethird = 1.0 / 3.0;
const double twothirds = 2.0 / 3.0;
const double fourthirds = 4.0 / 3.0;
