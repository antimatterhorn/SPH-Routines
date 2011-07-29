#include <stdio.h>

/* This code calculates the necessary positions for a two-body system	*/
/* for roche overflow of body 1.										*/

double dist_in_cm;
double mass_in_g;
double mass_in_msun;
double dens_in_gccm;
double energy_in_erg;
double time_in_s;
double specenergy_in_ergperg;
double pressure_in_ergperccm;
double Gnewt;

float m1,m2,x1,x2,r1,r2;
float a,b,c;
float tol;
int i,j;

float lhs(float x)
{
	a = (m1+m2)*Gnewt / pow(x*(1+m1/m2),3);
	b = x-r1;
	c = m1*Gnewt / pow(r1,2);
	return a*b+c;
}

float rhs(float x)
{
	return m2*Gnewt / pow(x*(1+m1/m2)-r1,2);
}

int main(int argc, char **argv[])
{
	
	time_in_s				= 1;
	dist_in_cm				= 6.955e7;
	mass_in_g				= 1.989e27;
	mass_in_msun			= mass_in_g / 1.99e33;
	dens_in_gccm			= mass_in_g/dist_in_cm/dist_in_cm/dist_in_cm;
    energy_in_erg			= mass_in_g*dist_in_cm*dist_in_cm/time_in_s/time_in_s;
    specenergy_in_ergperg	= energy_in_erg/mass_in_g;
    pressure_in_ergperccm	= energy_in_erg/dist_in_cm/dist_in_cm/dist_in_cm;
	Gnewt					= 0.000393935;
	
	if (argc < 2){
		printf("m1(sol): ");
		scanf("%f", &m1);
	}else {
		m1 = atof(argv[1]);
	}
	
	if (argc < 3){
		printf("m2(sol): ");
		scanf("%f", &m2);
	}else {
		m2 = atof(argv[2]);
	}

	m1 = m1 / mass_in_msun;
	m2 = m2 / mass_in_msun;
	
	printf("m1(code) = %f\nm2(code) = %f\n",m1,m2);
	
	if (argc < 4){
		printf("r1(code): ");
		scanf("%f", &r1);
	}else {
		r1 = atof(argv[3]);
	}
	
	if (argc < 5){
		printf("tol(<0.1): ");
		scanf("%f", &tol);
	}else {
		tol = atof(argv[4]);
	}
	
	j = 1;
	
	for (i = 0; i < 100000;i++)
	{
		x1 = (float)i / (float)100;
		if (lhs(x1)/rhs(x1) < (1+tol) && lhs(x1)/rhs(x1) > (1-tol))
		{
			x2 = x1*m1/m2;
			printf("solution %d\n",j);
			printf("\tx1 = %3.3f\n\tx2 = %3.3f\n",x1,x2);
			j++;
		}
	}
	
	//printf("x1 = %f\nx2 = %f\n",x2*(m2/m1),x2);
	
	return 0;
}


