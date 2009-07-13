#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>

//using std::cout;
//using std::cin;
//using std::endl;

using namespace std;

const double pi = 3.1415926;
const double G = 3.93935e-7; //in SPH code units

const double onethird = 1.0 / 3.0;
const double twothirds = 2.0 / 3.0;
const double fourthirds = 4.0 / 3.0;

double b = 20.0; //Rsun
double db = 0.1;
double m1,m2;
double r12;
double r1[2],r2[2];
double v1[2],v2[2];
double v = 50;
double h = 1.0;
double ca = 0.005;
double rmin = 50.0;

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
	

	cout << "Accretor Mass (Msun) = ";
	cin >> m1;
	cout << "Impactor Mass (Msun) = ";
	cin >> m2;
	
	cout << "x1(Rsun) =";
	cin >> r1[0];
	cout << "y1 =";
	cin >> r1[1];
	cout << "vx1(Rsun/s) =";
	cin >> v1[0];
	cout << "vy1 =";
	cin >> v1[1];
	
	cout << "x2 =";
	cin >> r2[0];
	cout << "y2 =";
	cin >> r2[1];
	cout << "vx2 =";
	cin >> v2[0];
	cout << "vy2 =";
	cin >> v2[1];


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
		
			if (r12 < rmin)
			{
				rmin = r12;
			}
			//cout << r12 << "\n";
		
		} while ((r12 <= rmin) or (r12 > 0.1));

	cout << "Closest Approach =" << r12 << "\n";
		
	return 0;
}


