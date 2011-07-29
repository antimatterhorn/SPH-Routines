#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

/* This code will read in the inputh.dat and inputmodel.dat files */
/* to produce the correct cutoff that corresponds to the desired mass */

using namespace std;

const double pi = 3.1415926;

const double onethird = 1.0 / 3.0;
const double twothirds = 2.0 / 3.0;
const double fourthirds = 4.0 / 3.0;

double mass,total,sub,r,dr;
int n,i;

int main()
{
	FILE *file2;
	
	ifstream file("inputh.dat") ;
	string line ;
	vector<string> lines ;
	while( getline( file, line ) ) lines.push_back( line ) ;
	n = lines.size();	
	cout << "inputh.dat contains " << n << " line(s)" << endl;
	
	cout << "Desired mass: ";
	cin >> mass;
		
	double input[n][3];
	/* inputh.dat will contain r and h, while inputmodel will contain r, rho, u ... */
	/* here we want just r, h, rho */
	
	float d1,d2,d3;
	file2 = fopen("inputh.dat","r");
	i=0;
	while (fscanf(file2,"%g %g", &d1,&d2)!=EOF) {
		//printf("%e %e\n",d1,d2);
		input[i][0] = d1;
		input[i][1] = d2;
		i++;
	}
	fclose(file2);
	
	file2 = fopen("inputmodel.dat", "r");
	i=0;
	while (fscanf(file2,"%g %g %*g %*g %*g", &d1,&d2)!=EOF) {
		if (d1 != input[i][0]) {
			/* there's a problem here */
			printf("r in inputmodel.dat does not match r in inputh.dat for row %d\n",i);
		}
		
		input[i][2] = d2;	
		i++;
	}
	fclose(file2);
	
//	for (i=0;i<n;i++)
//	{
//		fscanf(file2,"%f %f %f %f %f\n",&d1,&d2,&d3,&d3,&d3);
//		if (d1 != input[i][0]) {
//			/* there's a problem here */
//			//printf("r in inputmodel.dat does not match r in inputh.dat for row %d\n",i);
//		}
//		input[i][2] = d2;
//	}
//	fclose(file2);
	for (i=0;i<n-1;i++) {
		r = (input[i+1][0] + input[i][0])*0.5;
		dr = (input[i+1][0] - input[i][0]);
		sub = 4*pi*(input[i+1][2]+input[i][2])*dr*0.5*pow(r,2.0);
		//printf("%e %e %e\n",input[i][0],input[i][1],input[i][2]);
		total += sub;
	}
	printf("total mass:%e\n",total);
	return 0;
}


