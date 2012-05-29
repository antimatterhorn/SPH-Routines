/*
 *  sdf_masscutoff.c
 *  
 *
 *  Created by Cody Raskin on 7/14/09.
 *  reads an sdf file and outputs a subset with a chosen mass
 *
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>


double abundin[13],abundout[13];

double za[13] = {2.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28.};
double aa[13] = {4.,12.,16.,20.,24.,28.,32.,36.,40.,44.,48.,52.,56.};
//double ma[13] = {4.00260325415,12.,15.99491461956,20.,24.,27.9769265325,32.,36.,39.96259098,44.,48.,52.,55.942132};
double ea[13] = {0.005974,0.017909,0.023871,0.029837,0.035796,0.041753,0.047716,0.053679,0.059641,0.065606,0.071568,0.077528,0.083489};

double igex[14] = {0.,3.642e6,4.098e6,4.963e6,6.1e6,6.763e6,7.499e6,8.313e6,9.354e6,1.150e7,
                    1.294e7,1.477e7,1.898e7,1.e10};
double igey[14] = {0.,0.003,0.014,0.033,0.052,0.071,0.108,0.184,0.373,0.691,0.831,0.917,1.0,1.0};

double imex[30] = {0.,4.9e4,6.4e4,7.89e4,1.49e5,2.314e5,2.642e5,2.973e5,3.601e5,5.943e5,1.279e6,
                    1.549e6,1.795e6,2.020e6,2.306e6,2.595e6,2.877e6,3.285e6,3.696e6,4.890e6,6.191e6,
                    6.763e6,7.499e6,8.192e6,9.082e6,1.133e7,1.294e7,1.499e7,1.955e7,1.e10};
double imey[30] = {0.,0.,0.010,0.040,0.260,0.388,0.407,0.415,0.422,0.430,0.441,0.452,0.475,0.532,
                    0.660,0.834,0.944,0.993,0.997,0.970,0.955,0.933,0.895,0.815,0.706,0.331,0.169,
                    0.089,0.,0.};

double carbx[12] = {0.,3.95e4,5.2995e4,6.51e4,8.13e4,1.21e5,1.444e5,1.698e5,2.247e5,2.681e5,3.933e5,1.e10};
double carby[12] = {0.5,0.498,0.479,0.449,0.388,0.180,0.101,0.063,0.021,0.010,0.0,0.};

double oxygx[20] = {0.,5.222e4,6.419e4,1.210e5,1.423e5,1.698e5,2.214e5,2.603e5,3.107e5,5.205e5,
                    1.260e6,1.397e6,1.549e6,1.743e6,2.050e6,2.306e6,2.963e6,3.285e6,3.696e6,1.e10};
double oxygy[20] = {0.5,0.517,0.539,0.634,0.641,0.634,0.596,0.581,0.577,0.573,0.566,0.562,0.547,0.528,
                    0.468,0.350,0.059,0.010,0.,0.};

double igeabund[13] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0};
double imeabund[13] = {0.0,0.0,0.00003,0.0,0.00002,0.57207,0.28019,0.04144,0.02738,0.00020,0.00329,0.07202,0.00337};
double carbabund[13] = {0.0,0.56717,0.43024,0.00118,0.00037,0.00034,0.00016,0.00003,0.00002,0.0,0.0,0.00045,0.00003};
double oxygabund[13] = {0.0,0.00106,0.70143,0.00013,0.06546,0.15346,0.07014,0.00701,0.00019,0.00009,0.00007,0.00031,0.00064};

int usage()
{
	printf("\t Returns a burn table for 13 isotopes.\n");
	printf("\t Usage: [required] {optional}\n");
	printf("\t abundtest [1=go]\n");
	return 0;
}

int burn(double rhoin)
{
    int i;
    double atot=0,mix[4];
    
    for(i=0;i<4;i++) mix[i] = 0;
    
    for(i=0;i<12;i++) if(carbx[i]>rhoin) break;
    mix[0] = (carby[i]-carby[i-1])/(carbx[i]-carbx[i-1])*(rhoin-carbx[i-1])+carby[i-1];
    
    for(i=0;i<20;i++) if(oxygx[i]>rhoin) break;
    mix[1] = (oxygy[i]-oxygy[i-1])/(oxygx[i]-oxygx[i-1])*(rhoin-oxygx[i-1])+oxygy[i-1];
    
    for(i=0;i<30;i++) if(imex[i]>rhoin) break;
    mix[2] = (imey[i]-imey[i-1])/(imex[i]-imex[i-1])*(rhoin-imex[i-1])+imey[i-1];
    
    for(i=0;i<14;i++) if(igex[i]>rhoin) break;
    mix[3] = (igey[i]-igey[i-1])/(igex[i]-igex[i-1])*(rhoin-igex[i-1])+igey[i-1];
    
    for(i=0;i<4;i++) atot += mix[i];
    for(i=0;i<4;i++) mix[i] = mix[i]/atot;
    
    for(i=0;i<13;i++) 
		abundout[i] = (2.0*mix[0])*carbabund[i]+
		(mix[1]-mix[0])*oxygabund[i]+
		(mix[2])*imeabund[i]+
		(mix[3])*igeabund[i];  

    atot = 0;
    for(i=0;i<13;i++)
    {
        if(abundout[i]<0) abundout[i] = 0;
        atot += abundout[i];
    }
        
    if(atot>0) for(i=0;i<13;i++) abundout[i] = abundout[i]/atot;
    
    return 0;
}

void quicksort(double** data,int m,int n)
{
    int key,i,j,k;
    double *temp;
    if( m < n)
    {
        key = m;
        i = m;
        j = n;
        while(i < j)
        {
            while((i < n) && (data[i][0] <= data[key][0]))
                i++;
            while(data[j][0] > data[key][0])
                j--;
            if( i < j)
            {
                temp = data[i];
                data[i] = data[j];
                data[j] = temp;
            }
                
        }
        temp = data[key];
        data[key] = data[j];
        data[j] = temp;

        quicksort(data,m,j-1);
        quicksort(data,j+1,n);
    }
}

int main(int argc, char **argv[])
{
	if(argc < 2)
	{
		usage();
		return 0;
	}
	
	int i,j,k,id;    
    double drho = (1.e8)/10000.;
    double rho;
    
    for(i=0;i<10000;i++)
    {
        rho += drho;
        burn(rho);
        printf("%3.2e ",rho);
        for(j=0;j<13;j++) printf("%1.4f ",abundout[j]);
        printf("\n");
    }
    
    return 0;
}