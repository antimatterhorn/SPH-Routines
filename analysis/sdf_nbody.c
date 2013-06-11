#include "../SPHbody.h"
#include <stdio.h>
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include <stddef.h>
#include "Msgs.h"
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)



char sphfile[80];
char outfile[80];
char outfile2[80];
char cmd[80];

int usage()
{
	printf("\t Runs an n-body calculation on sph particles\n");
	printf("\t Usage: [required] {optional,default}\n");
	printf("\t sdf_nbody [sdf file] [central mass]\n");
	return 0;
}


int main(int argc, char **argv[])
{
	int gnobj, nobj;
	int conf;
	float tpos,tend,dt,tout,dtn;
	float ax,ay,az,r,rx,ry,rz,cmass,rmin,amax,vmin,vdotr;
	int i,j,itr=0,ditr;
	double mass_in_g = 1.989E+27;        
	double dist_in_cm = 6.955e7;  
	double Gnewt;
	double time_in_s = 1;
	double energy_in_erg,dens_in_gccm,pressure_in_ergperccm,specenergy_in_ergperg;
	
	if (argc < 3)
	{
		usage();
		return 0;
	}
	
    
	/* try sample data first */
    /*
    cmass = atof(argv[2]);
    nobj = 2;
    SPHbody *body = malloc(nobj*sizeof(SPHbody));
    
    Gnewt = 39.53;
    
    //earth
	body[0].x = 1.0;
	body[0].y = body[0].z = 0.0;
	body[0].vy = 2*3.14159;
	body[0].vx = body[0].vz = 0.0;
	
    //mars
	body[1].x = 1.523;
	body[1].y = body[1].z = 0.0;
	body[1].vy = 2*3.14159*1.523/1.88;
	body[1].vx = body[1].vz = 0.0;
    
	body[0].mass = 3e-6;
    body[1].mass = 3e-6;

    */
    
    cmass = atof(argv[2])*1.0005e6;
     
    SDF *sdfp;
	SPHbody *body;
	
    
     
	sdfp = SDFreadf(argv[1], (void **)&body, &gnobj, &nobj, sizeof(SPHbody),
					"x", offsetof(SPHbody, x), &conf,
					"y", offsetof(SPHbody, y), &conf,
					"z", offsetof(SPHbody, z), &conf,
					"mass", offsetof(SPHbody, mass), &conf,
					"vx", offsetof(SPHbody, vx), &conf,
					"vy", offsetof(SPHbody, vy), &conf,
					"vz", offsetof(SPHbody, vz), &conf,
					NULL);
	SDFgetfloatOrDefault(sdfp, "tpos",  &tpos, (float)0.0);
    
	SDFclose(sdfp);
     
    Gnewt = 0.000393935;

    dt = 1.0e-5;
	/* test params */
	//tpos = 0;
	tout = 10.0;
	
		
	printf("tpos = %3.5f\n",tpos);
    printf("tend = ");
    scanf("%f",&tend);
	
	
	snprintf(outfile, sizeof(outfile), "nbody.dat");
	snprintf(outfile2, sizeof(outfile2), "nbody2.dat");
	FILE *stream, *fopen(), *stream2;
	/* declare a stream and prototype fopen */ 
	
	stream = fopen(outfile,"w");	
	stream2 = fopen(outfile2,"w");
	fprintf(stream,"t,dt,x,y,z,vx,vy,vz,ax,ay,az\n");
	fprintf(stream2,"t,dt,x,y,z,vx,vy,vz,ax,ay,az\n");
    
    ditr = 100;
    
    fprintf(stream,"%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e\n",
            tpos,dt,
            body[0].x,body[0].y,body[0].z,
            body[0].vx,body[0].vy,body[0].vz,
            body[0].ax,body[0].ay,body[0].az);
    fprintf(stream2,"%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e\n",
            tpos,dt,
            body[1].x,body[1].y,body[1].z,
            body[1].vx,body[1].vy,body[1].vz,
            body[1].ax,body[1].ay,body[1].az);
    
    while (tpos<tend) {
        itr++;
		tpos += dt;
		
		dtn = 1e10;
        rmin= 1e10;
		amax= 0;
		vmin= 1e10;
#pragma omp parallel for \
        shared(body,dt,nobj) \
        default(none)
        for(i=0;i<nobj;i++)
        {			
            body[i].x += body[i].vx*dt;
            body[i].y += body[i].vy*dt;
            body[i].z += body[i].vz*dt;
        }
#pragma omp parallel for \
    private(j,rx,ry,rz,r,ax,ay,az,vdotr) \
    shared(body,dt,nobj,cmass,Gnewt,dtn,rmin,amax,vmin) \
    default(none)
        for(i=0;i<nobj;i++)
        {       
            ax = ay = az = 0;
            
            for(j=0;j<nobj;j++)
            { 
                if(i!=j)
                {
                    rx = body[j].x - body[i].x;
                    ry = body[j].y - body[i].y;
                    rz = body[j].z - body[i].z;
                    r = sqrt(pow(rx,2.0)+
                             pow(ry,2.0)+
                             pow(rz,2.0));
                    
                    rmin = MIN(rmin,r);
                    
                    ax += body[j].mass/pow(r,3.0)*rx;
                    ay += body[j].mass/pow(r,3.0)*ry;
                    az += body[j].mass/pow(r,3.0)*rz;
                }
            }
            
            /* central mass attraction */
            ax -= cmass*body[i].x/pow(sqrt(pow(body[i].x,2.0)+
                                           pow(body[i].y,2.0)+
                                           pow(body[i].z,2.0)),3.0);
            ay -= cmass*body[i].y/pow(sqrt(pow(body[i].x,2.0)+
                                           pow(body[i].y,2.0)+
                                           pow(body[i].z,2.0)),3.0);
            az -= cmass*body[i].z/pow(sqrt(pow(body[i].x,2.0)+
                                           pow(body[i].y,2.0)+
                                           pow(body[i].z,2.0)),3.0);
            
            ax = ax*Gnewt;
            ay = ay*Gnewt;
            az = az*Gnewt;
			
			body[i].ax = ax;
			body[i].ay = ay;
			body[i].az = az;
			
			amax = MAX(amax,sqrt(ax*ax+ay*ay+az*az)*dt/
						   sqrt(pow(body[i].vx,2.0)+
								pow(body[i].vy,2.0)+
								pow(body[i].vz,2.0)));
            
            body[i].vx += ax*dt;
            body[i].vy += ay*dt;
            body[i].vz += az*dt;

            //vdotr = body[i].vx*rx + body[i].vy*ry + body[i].vz*rz;
            vdotr = sqrt(pow(body[i].vx,2.0)+
                         pow(body[i].vy,2.0)+
                         pow(body[i].vz,2.0));
            if(vdotr > 0) 
            {
                dtn = MIN(dtn,0.1*rmin/vdotr);
            }
        }

        if(itr % ditr == 0) 
        {
			printf("iter %04d: tpos = %3.2e \t dt=%3.2e-->%3.2e \t %3.2e\n",itr,tpos,dt,dtn,amax);
            fprintf(stream,"%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e\n",
                    tpos,dtn,
                    body[0].x,body[0].y,body[0].z,
                    body[0].vx,body[0].vy,body[0].vz,
                    body[0].ax,body[0].ay,body[0].az);
            fprintf(stream2,"%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e\n",
                    tpos,dtn,
                    body[1].x,body[1].y,body[1].z,
                    body[1].vx,body[1].vy,body[1].vz,
                    body[1].ax,body[1].ay,body[1].az);
		
		}
		dt = dtn;

    }
    fclose(stream);
	fclose(stream2);
    
    snprintf(outfile, sizeof(outfile), "nbody.%04d", itr);
	SDFwrite(outfile, nobj, 
			 nobj, body, sizeof(SPHbody),
			 SPHOUTBODYDESC,
			 "npart", SDF_INT, nobj,
			 "iter", SDF_INT, itr,
			 "dt", SDF_FLOAT, dt,
			 "Gnewt", SDF_FLOAT, Gnewt,
			 "ndim", SDF_INT, 3,
			 "tpos", SDF_FLOAT, tpos,
			 NULL);	
    
    
    free(body);
    
	return 0;
}