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

double b = 1.0; //Rsun
double m1,m2;
double r12;
double r1[2],r2[2];
double v1[2],v2[2];
double v = 50;
double h = 1.0;

double a1(int j)
{
	return -G*m2*(r1[j]-r2[j]) / pow(r12,3.0);
}

double a2(int j)
{
	return G*m1*(r1[j]-r2[j]) / pow(r12,3.0);
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
	
	v1[0] = 0;
	v1[1] = 0;
	v2[0] = 0.000001*v;	//  km/s -> Rsun/s
	v2[1] = 0;
	
	cout << "Accretor Mass (Msun) = ";
	cin >> m1;
	cout << "Impactor Mass (Msun) = ";
	cin >> m2;

	//open the stream file
	ofstream pout("positions.csv");
	pout << "x1,y1,x2,y2" << endl;
	
	double i = 0;
	do //Euler's Method
	{
        
		r1[0] += h*v1[0];
		r2[0] += h*v2[0];
		r1[1] += h*v1[1];
		r2[1] += h*v2[1];
		
		r12 = sqrt(pow(r1[0]-r2[0],2.0)+pow(r1[1]-r2[1],2.0));
		
		v1[0] += h*a1(0);
		v2[0] += h*a2(0);
		
		v1[1] += h*a1(1);
		v2[1] += h*a2(1);
		
		i += h;

		pout << r1[0] << "," << r1[1] << "," 
			<< r2[0] << "," << r2[1] << endl;
		
	} while (r12 > 0.1);

	pout.close();
	
	//recenter on CM
	double cm[2];
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
	
	return 0;
}


