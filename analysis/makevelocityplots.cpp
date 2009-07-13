#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>

//using std::cout;
//using std::cin;
//using std::endl;

using namespace std;

const double pi = 3.1415926;
const double onethird = 1.0 / 3.0;
const double twothirds = 2.0 / 3.0;
const double fourthirds = 4.0 / 3.0;

int npart;
double x,vx,z,h,rho,temp,pr,cs;
char filename[80];
char output[80];

int main()
{
	FILE *file2;
	
	//do preliminary setup
	cout << "Input file name: ";
	cin >> filename;	
	ifstream file(filename) ;
	string line ;
	vector<string> lines ;
	while( getline( file, line ) ) lines.push_back( line ) ;
	//cout << "#lines: " << lines.size() << '\n' ;
	npart = lines.size();	
	cout << filename << " contains " << npart << " particle(s)" << endl;
	cout << "Output file name: ";
	cin >> output;
	
	file2 = fopen(filename,"r");
	//here is where the loading of the data should happen (x,vx,z,h,rho,temp,pr)
	double data[npart][7];
	float d1,d2,d3,d4,d5,d6,d7;
	for(int i=0;i<npart;i++)
	{
		fscanf(file2,"%f%f%f%f%f%f%f\n",&d1,&d2,&d3,&d4,&d5,&d6,&d7);
		//printf("%6.4e %6.4e %6.4e %6.4e %6.4e %6.4e %6.4e \n",
		//	d1,d2,d3,d4,d5,d6,d7);
		data[i][0]=d1;
		data[i][1]=d2;
		data[i][2]=d3;
		data[i][3]=d4;
		data[i][4]=d5;
		data[i][5]=d6;
		data[i][6]=d7;
	}
	fclose(file2);
	
	//false data for now
	cout << "Data Loaded\n";
	
	double out[npart][3];
			
	for(int p = 0; p < npart; p++)
	{
		cout << p << "/" << npart-1 << endl;

			
			x = data[p][0];
			vx = data[p][1];
			z = data[p][2];
			h = data[p][3];
			rho = data[p][4];
			temp = data[p][5];
			pr = data[p][6];
			
			cs = sqrt(fourthirds*pr/rho)*696*pow(10.0,3.0);
		
		if (fabs(z) < h)
		{
			out[p][0] = x;
			out[p][1] = fabs(vx)*696*pow(10.0,3.0);
			out[p][2] = cs;	
			
			
		}
		

		

	}
	
	//open the stream file
	ofstream fout(output);
	
	fout << "x,vx,cs\n";

	for(int i = 0; i < npart; i++)
	{
		fout << out[i][0] << "," << out[i][1] << "," << out[i][2] << endl;
	}
	
	//close the stream file
	fout.close();

	return 0;
}


