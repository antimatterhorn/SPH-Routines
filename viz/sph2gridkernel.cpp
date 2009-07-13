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

int choice,pixels,npart;
double xmax,ymax;
double x,y,z,h,rho,temp,mass;
char filename[80];
char output[80];

int x2i(double xx)
{
	int i;
	xx = xx*(pixels/(2.0*xmax));
	i = xx + pixels/2.0;	
	return i;
}

int y2j(double yy)
{
	int j;
	yy = -yy*(pixels/(2.0*ymax));
	j = yy + pixels/2.0;
	return j;
}

double w(double hi, double ri)
{
	double v = ri/hi;
	double coef = pow(pi,-1.0)*pow(hi,-3.0);
	if (v <= 1 and v >= 0)
	{
		return coef*(1.0-1.5*v*v+0.75*v*v*v);
	}
	else if (v > 1 and v <= 2)
	{
		return coef*0.25*pow(2.0-v,3.0);
	}
	else
	{
		return 0;
	}
}


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
	
	cout << "\nOutput Parameters \n-----------------\nChoose a Variable:\n";
	cout << "1) Density \n2) Temperature\n\nChoice:";
	cin >> choice;
	choice += 3;
	
	cout << "Pixels on a side:";
	cin >> pixels;
	
	cout << "Xmax:";
	cin >> xmax;
	ymax = xmax;

	//file.close;
	file2 = fopen(filename,"r");
	//here is where the loading of the data should happen (x,y,z,h,rho,temp,mass)
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
	
	double **dens;
	
	dens = (double **) malloc((size_t)pixels*sizeof(double *));
	for(int i=0;i<pixels;i++){
		dens[i]  = (double *) malloc((size_t)pixels*sizeof(double));
	}
	
	double **img;
	
	img = (double **) malloc((size_t)pixels*sizeof(double *));
	for(int i=0;i<pixels;i++){
		img[i]  = (double *) malloc((size_t)pixels*sizeof(double));
	}
	
	//double img[pixels][pixels];
	//double dens[pixels][pixels];
	
	//fill arrays with 0s?
	for(int i = 0; i < pixels; i++)
	{
		for(int j = 0; j < pixels; j++)
		{
			img[i][j]=dens[i][j]=0;
		}
	}
		
	for(int p = 0; p < npart; p++)
	{
		cout << p << "/" << npart-1 << endl;
		
		if (fabs(data[p][2]) < data[p][3]) // |z| < h
		{
			
			x = data[p][0];
			y = data[p][1];
			z = data[p][2];
			h = data[p][3];
			rho = data[p][4];
			temp = data[p][5];
			mass = data[p][6];
			
			h = sqrt(pow(h,2.0)-pow(z,2.0));
			
			int ic = x2i(x);
			int jc = y2j(y);
			int r = h*(pixels/(2.0*xmax));
			
			if (r==0)
			{
				cout << "found that r was too small, filling only 1 pixel" << endl;
				if (ic >= 0 and ic < pixels and jc >= 0 and jc < pixels)
				{
					img[ic][jc] += data[p][choice]*rho;
					dens[ic][jc] += rho;	
				}
				
			}
			else
			{
				for(int imi = ic-2*r; imi < ic+2*r+1; imi++)
				{
					for(int imj = jc-2*r; imj < jc+2*r+1; imj++)
					{
						if (imi >= 0 and imi < pixels and imj >= 0 and imj < pixels) // inside the image
						{
							double rr = sqrt(pow((imi-ic),2.0)+pow((imj-jc),2.0));
							//cout << rr << "," << r << endl;
							if (rr <= 2*r) // inside the circle
							{							
								img[imi][imj] += data[p][choice]*rho*w(h,2*rr*(double)xmax/(double)pixels)*pi*pow(h,3.0);
								dens[imi][imj] += rho;						
							}
						}
			
					}
				}	
			}
		}
	}
	
	//open the stream file
	ofstream fout(output);
	double maxv = 0.0;
	double minv = pow(10.0,10.0);
	for(int i = 0; i < pixels; i++)
	{
		for(int j = 0; j < pixels; j++)
		{
			if (dens[i][j] != 0)
			{
				fout << img[i][j]/dens[i][j] << " ";
				if (img[i][j]/dens[i][j] < minv) {minv = img[i][j]/dens[i][j];}
				if (img[i][j]/dens[i][j] > maxv) {maxv = img[i][j]/dens[i][j];}
			}
			else
			{fout << 0.0 << " ";}
						
		}
		fout << endl;
	}
	
	//close the stream file
	fout.close();
	std::cout << "min,max = " << minv << " " << maxv << endl;
	return 0;
}


