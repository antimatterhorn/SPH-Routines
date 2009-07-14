#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>

/* This code will solve for the initial separation required */
/* to bring two bodies of masses m1,m2 to within a desired  */
/* separation on the first approach.						*/

using namespace std;

const double pi = 3.1415926;
const double G = 3.93935e-7; //in SPH code units

const double onethird = 1.0 / 3.0;
const double twothirds = 2.0 / 3.0;
const double fourthirds = 4.0 / 3.0;

double b = 25.0; //Rsun
double db = 0.1;
double m1,m2;
double r12;
double r1[2],r2[2];
double v1[2],v2[2];
double v = 50.0;
double h = 0.5;
double k1,k2,k3,k4,KE,PE,TE0,TE1;
double ca = 0.005;
double rmin = 50.0;
double c1 = (0.001 - 0.5) / 9.0;
double c2 = 0.5 - 10.0*c1;


double a10(double dr)
{
	return -G*m2*(r1[0]+dr-r2[0]) / pow(pow(r1[0]+dr-r2[0],2.0)+pow(r1[1]-r2[1],2.0),1.5);
}

double a11(double dr)
{
	return -G*m2*(r1[1]+dr-r2[1]) / pow(pow(r1[0]-r2[0],2.0)+pow(r1[1]+dr-r2[1],2.0),1.5);
}

double a20(double dr)
{
	return G*m1*(r1[0]-r2[0]-dr) / pow(pow(r1[0]-r2[0]-dr,2.0)+pow(r1[1]-r2[1],2.0),1.5);
}

double a21(double dr)
{
	return G*m1*(r1[1]-r2[1]-dr) / pow(pow(r1[0]-r2[0],2.0)+pow(r1[1]-r2[1]-dr,2.0),1.5);
}

int main()
{
	//do preliminary setup
	
	cout << "Desired Closest Approach (Rsun): ";
	cin >> ca;
	
	cout << "Guess: ";
	cin >> b;
	
	cout << "Velocity (km/s): ";
	cin >> v;

	cout << "Accretor Mass (Msun) = ";
	cin >> m1;
	cout << "Impactor Mass (Msun) = ";
	cin >> m2;

	do
	{
		r1[0] = 0;
		r1[1] = 0;
		r2[0] = -10;
		r2[1] = b;
		
		r12 = sqrt(pow(r1[0]-r2[0],2.0)+pow(r1[1]-r2[1],2.0));
	
		v1[0] = 0;
		v1[1] = 0;
		v2[0] = 0.000001*v;	//  km/s -> Rsun/s
		v2[1] = 0;
		
		rmin = 50.0;
		h = 0.5;
	
		do //RK
		{
			//TE0 = KE+PE;
			
			r1[0] += h*v1[0];
			r2[0] += h*v2[0];
			r1[1] += h*v1[1];
			r2[1] += h*v2[1];
			
			r12 = sqrt(pow(r1[0]-r2[0],2.0)+pow(r1[1]-r2[1],2.0));
			
			k1 = h*a10(0);
			k2 = h*a10(0.5*k1);
			k3 = h*a10(0.5*k2);
			k4 = h*a10(k3);
			
			v1[0] += onethird*(0.5*k1 + k2 + k3 + 0.5*k4);
			
			k1 = h*a20(0);
			k2 = h*a20(0.5*k1);
			k3 = h*a20(0.5*k2);
			k4 = h*a20(k3);
			
			v2[0] += onethird*(0.5*k1 + k2 + k3 + 0.5*k4);
			
			k1 = h*a11(0);
			k2 = h*a11(0.5*k1);
			k3 = h*a11(0.5*k2);
			k4 = h*a11(k3);
			
			v1[1] += onethird*(0.5*k1 + k2 + k3 + 0.5*k4);
			
			k1 = h*a21(0);
			k2 = h*a21(0.5*k1);
			k3 = h*a21(0.5*k2);
			k4 = h*a21(k3);
			
			v2[1] += onethird*(0.5*k1 + k2 + k3 + 0.5*k4);
			

			
			//KE = 0.5*m1*(pow(v1[0],2.0)+pow(v1[1],2.0)) + 0.5*m2*(pow(v2[0],2.0)+pow(v2[1],2.0));
			//PE = -G*m1*m2/r12;
			//TE1 = KE+PE;
			
			h = c1*fabs(r2[0]-r1[0]) +c2;
			
			if (r12 < rmin)
			{
				rmin = r12;
			}

		} while (r12 <= rmin);

	b-= db;
	cout << fixed << setprecision(1) << b << " --> " << fixed << setprecision(5) << rmin << endl;

	} while (rmin > ca);
	
	cout << "Impact Parameter = " << b << "-" << b+db << endl;
	
	return 0;
}


