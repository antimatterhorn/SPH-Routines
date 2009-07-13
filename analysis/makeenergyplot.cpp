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

const double econv = (5.99e44)/(5.87e-9);

int npart;
double t,p;
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
	file.close();
	cout << filename << " contains " << npart << " line(s)" << endl;
	cout << "Output file name: ";
	cin >> output;
	
	file2 = fopen(filename,"r");
	//here is where the loading of the data should happen (t,ke,pe,U,te)
	double data[npart][5];
	float d1,d2,d3,d4,d5;
	for(int i=0;i<npart;i++)
	{
		fscanf(file2,"%f %f %f %f %f\n",&d1,&d2,&d3,&d4,&d5);
		//cout << d1 << " " << d2 << " " << d3 << " " << d4 << endl;
		data[i][0]=d1;
		data[i][1]=d2*econv;
		data[i][2]=d3*econv;
		data[i][3]=d4;
		data[i][4]=(d2+d3)*econv;
	}
	fclose(file2);
	
	//false data for now
	cout << "Data Loaded\n";
	
	//open the stream file
	ofstream fout(output);
	
	fout << "t,ke,pe,u,te\n";

	for(int i = 0; i < npart; i++)
	{
		fout << data[i][0] << "," << data[i][1] << "," << data[i][2] << "," << data[i][3] << "," << data[i][4] << endl;
	}
	
	//close the stream file
	fout.close();

	return 0;
}


