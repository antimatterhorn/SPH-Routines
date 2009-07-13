#include <iostream>
#include <math.h>
#include <fstream>

//using std::cout;
//using std::cin;
//using std::endl;

using namespace std;

const double pi = 3.1415926;
const double G = 3.93935e-7; //in SPH code units

const double onethird = 1.0 / 3.0;
const double twothirds = 2.0 / 3.0;
const double fourthirds = 4.0 / 3.0;
const double onesixth = 1.0 / 6.0;

double b = 1.0; //Rsun
double m1,m2;
double r12;
double r1[2],r2[2];
double v1[2],v2[2];
double v = 50.0;
double h = 0.5;
double k1,k2,k3,k4,KE,PE,TE0,TE1;
double cm[2];
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
	
	cout << "Impact Parameter (Rsun): ";
	cin >> b;
	
	cout << "Velocity (km/s): ";
	cin >> v;
	
	r1[0] = 0;
	r1[1] = 0;
	r2[0] = -10;
	r2[1] = b;
	
	r12 = sqrt(pow(r1[0]-r2[0],2.0)+pow(r1[1]-r2[1],2.0));
	
	v1[0] = 0;
	v1[1] = 0;
	v2[0] = 0.000001*v;	//  km/s -> Rsun/s
	v2[1] = 0;
	
	cout << "Accretor Mass (Msun) = ";
	cin >> m1;
	cout << "Impactor Mass (Msun) = ";
	cin >> m2;
	
	KE = 0.5*m1*(pow(v1[0],2.0)+pow(v1[1],2.0)) + 0.5*m2*(pow(v2[0],2.0)+pow(v2[1],2.0));
	PE = -G*m1*m2/r12;
	
	cout << "E = " << KE+PE << endl;
	
	
	//open the stream file
	ofstream pout("positions.csv");
	pout << "x1,y1,x2,y2" << endl;
	
	ofstream eout("energy.csv");
	eout << "t,E" << endl;
	
	double t = 0;
	int i = 0;
	do 
	{
        TE0 = KE+PE;
		
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
		

		
		KE = 0.5*m1*(pow(v1[0],2.0)+pow(v1[1],2.0)) + 0.5*m2*(pow(v2[0],2.0)+pow(v2[1],2.0));
		PE = -G*m1*m2/r12;
		TE1 = KE+PE;
		
		t+=h;
		i+=1;
		
		if (fmod(i, 20.0) ==0)
		{
			pout << r1[0] << "," << r1[1] << "," 
			<< r2[0] << "," << r2[1] << endl;
			eout << i << "," << TE1 << endl;
		}
		
		h = c1*fabs(r2[0]-r1[0]) +c2;
				
	//} while (t<pow(10,6.0));
	} while (r12 > 1.0);

	pout.close();
	eout.close();
	
	//recenter on CM
	cm[0] = (m1*r1[0]+m2*r2[0])/(m1+m2);
	cm[1] = (m1*r1[1]+m2*r2[1])/(m1+m2);
	r2[0] -= cm[0];
	r2[1] -= cm[1];
	r1[0] -= cm[0];
	r1[1] -= cm[1];
	
	cout << "x1=" << r1[0] << "\ty1=" << r1[1] << 
	"\tv1x=" << v1[0] << "\tv1y=" << v1[1] << endl;
	cout << "x2=" << r2[0] << "\ty2=" << r2[1] << 
	"\tv2x=" << v2[0] << "\tv2y=" << v2[1] << endl;
	
	cout << "E = " << TE1 << endl;
	
	return 0;
}


