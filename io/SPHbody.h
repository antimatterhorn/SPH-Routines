#define SPHOUTBODYDESC \
"struct {\n\
	double x, y, z;		/* position of body */\n\
	float mass;	   /* mass of body */\n\
	float vx, vy, vz;     /* velocity of body */\n\
	float u;      	   /* specific energy of body*/\n\
	float h;      	   /* smoothing length of body */\n\
	float rho;            /* density of body */\n\
	float pr;            /* pressure of body */\n\
	float drho_dt;        /* drho/dt of body */\n\
	float udot;           /* du/dt of body */\n\
	float temp;           /* temperature of body */\n\
	float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */\n\
	float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */\n\
	float vsound;		/* sound speed */\n\
	float abar;           /* avg number of nucleons per particle of body */\n\
	float zbar;           /* avg number of protons per particle of body */\n\
	float ax, ay, az;     /* acceleration of body */\n\
	float lax, lay, laz;  /* last acceleration of body */\n\
	float phi;            /* potential at body location */\n\
	float dt;           /* timestep */\n\
	float bdot;			/* energy generation rate */\n\
	float kappa;			/* opacity */\n\
	float sigma;			/* conductivity */\n\
	unsigned int nbrs;     /* number of neighbors */\n\
	unsigned int ident;	   /* unique identifier */\n\
	unsigned int windid;   /* wind id */\n\
	unsigned int useless;	/* to fill double block */\n\
}"

typedef struct {
	double x, y, z;             /* position of body */							//3
	float mass;           /* mass of body */
	float vx, vy, vz;     /* velocity of body */								//4
	float u;              /* specific energy of body*/
	float h;              /* smoothing length of body */						//6
	float rho;            /* density of body */
	float pr;            /* pressure of body */
	float drho_dt;        /* drho/dt of body */
	float udot;           /* du/dt of body */
	float temp;           /* temperature of body */								//11
	float He4, C12, O16, Ne20, Mg24, Si28, S32; /* abundances of body */		//18
	float Ar36, Ca40, Ti44, Cr48, Fe52, Ni56; /* abundances of body */			//24
	float vsound;
	float abar;           /* avg number of nucleons per particle of body */
	float zbar;           /* avg number of protons per particle of body */
	float ax, ay, az;     /* acceleration of body */							//30
	float lax, lay, laz;  /* last acceleration of body */						//33
	//float gax, gay, gaz;  /* gravity acceleration of body */
	//float grav_mass;      /* gravitational mass of body */
	float phi;            /* potential at body location */					
	//float tacc;           /* time of last acceleration update of body */
	float idt;
	float bdot;
	float kappa;			/* opacity */
	float sigma;			/* conductivity */									//38
	unsigned int nbrs;     /* number of neighbors */
	unsigned int ident;    /* unique identifier */
	unsigned int windid;   /* wind id */
	unsigned int useless;  /* to fill the double block (4int=1double)*/		
} SPHbody;


const double pi = 3.1415926;
const double onethird = 1.0 / 3.0;
const double twothirds = 2.0 / 3.0;
const double fourthirds = 4.0 / 3.0;