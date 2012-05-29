########################################################
# Input & output file names, settings
########################################################
int do_restart		=	0;
int independent_dt 	= 	0;		#newly added
char outfile[] 		= 	"fltmass";	#output name root
int nsteps 		= 	10000;
int save_first 		= 	1;
int output_freq 	= 	50;
int timer_freq 		= 	1000;
#char Msg_turn_on[] 	= 	"shrink.c,pqsort.c";
#int Msg_memfile 	= 	65536;



########################################################
# Unit conversion factors
########################################################
double DIST_IN_CM	=	6.955e7;	#0.001Rsun in cm
double TIME_IN_S 	=	1.;
double MASS_IN_G 	=	1.989E+27;	#10^-6 Msun in gm



########################################################
# BURNER SETTINGS
########################################################
int do_burner		=	1;		#not in mainline
double mintemp_burner	=	1e8; 		#not in mainline
int do_compos		=	0;
float ufloor		=	1.0e-6;
float dumax		=	0.3;		#maximum allowed changed of
      						#internal energy from burning
float dtumin		=	1.0;
# APROX13	##### He,C,O,Ne,Mg,Si,S,Ar,Ca,Ti,Cr,Fe,Ni
struct {float composition;}[13] = {0.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
# ISO7 	      	##### He4, C12, O16, Ne20, Mg24, Si28, Ni56
# struct {float composition;}[7] = {0.,1.,1.,0.,0.,0.,0.};



########################################################
# Input set #1, ignored if do_restart is set
########################################################
int do_sph		=	1;
int do_nbody 		= 	0;
char SPHdatafile[] 	= 	"ic.sdf";	#the initial conditions file
char Nbodyfile[] 	= 	"headon";
double x_offset		=	0.0;
double y_offset		=	0.0;
double z_offset		=	0.0;
double vx_offset	=	0.0;
double vy_offset	=	0.0;
double vz_offset	=	0.0;



########################################################
# Boundary SPH Input set, ignored if do_restart is set
########################################################
int do_sphboundary	=	0;
char SPHboundaryfile[] 	= 	"Hecore1kb.sdf";
double x_offset3	=	0.;
double y_offset3	=	0.;
double z_offset3	=	0.;
double vx_offset3	=	0.;
double vy_offset3	=	0.;
double vz_offset3	=	0.;
int do_posfixed		=	1;
int do_sphfixed		=	1;



########################################################
# Add a point mass (N-body only)
########################################################
int do_pointmass	=	0;
double x_offpoint	=	80.;
double y_offpoint	=	0.;
double z_offpoint	=	0.;
double vx_offpoint	=	-1e-5;
double vy_offpoint	=	0.;
double vz_offpoint	=	0.;
double mass_point	=	0.39197292;
double eps_point	=	-0.1;
int fixed_point		=	0;



########################################################
# Add an SPH+Nbody particle
########################################################
int do_sphpointmass	=	0;
double x_offsphpoint	=	83.;
double y_offsphpoint	=	0.;
double z_offsphpoint	=	0.;
double vx_offsphpoint	=	0.;
double vy_offsphpoint	=	6.68673e-5;
double vz_offsphpoint	=	0.;
double mass_sphpoint	=	1e-30;
double gravmass_sphpoint=	0.053;
double h_sphpoint	=	0.078297300;
double u_sphpoint	=	1e-30;
double rho_sphpoint	=	1e-30;
int fixed_sphpoint	=	0;
int sphfixed_sphpoint	=	1;
int dual_sphpoint	=	1;
int special0_sphpoint	=	0;
int special1_sphpoint	=	1;



########################################################
# Major options 
########################################################
int exact_rho		=	1;
int do_plummer 		= 	0;
int do_grav 		= 	1;
int do_gravSPH 		= 	1;
int has_grav_data 	= 	0;
int timeout 		= 	1200;
int set_id		=	1;
int do_kernelsymmetrymaxh=	0;	# Avoid particle diffusion, 
    					# 	hij=max(hi,hj)


########################################################
# Pick the Equation of State
########################################################
int do_eosidealgas	=	0;
int do_eospolytrope	=	0;
int do_eosnadyozhin	=	0;
int do_eoshelmholtz	=	1;	#uses helm_table.dat
double helm_small_temp	=	1e4; 	#low temp to not solve for states
       					#    not on the grid


########################################################
# Energy Transport
########################################################
int do_cooling		=	0;
int do_diffusion	=	0;
int do_thermdiff	=	1;	#CR's thermal diffusion routine
int do_radcool		=	1;	#CR's radiative cooling routine
float oang		=	1.04720;#opening angle for do_radcool

int do_thermalize	= 	1;	#this turns off most other physics!!
float therm_temp	= 	2.0e6;	#temperature about which to artificially thermalize
inter thermiter		= 	1000;


########################################################
# Timestep settings
########################################################
int adaptive_dt		=	1;
int do_timedrift	=	0;	#CR's physics speed-up routine
double dt 		= 	0.004;
double dark_dt 		= 	1e30;
double dt_max 		= 	32.0;
int tlow_cut 		= 	1;
int independent_dt	=	1;
int dark_independent_dt	=	0;
int maxnumdtlevels	=	1;	# Max nr. of individual dt levels, 
    					#     remember 2^(n-1) subs
double maxdtfrac	=	0.125; 	# At least maxdtfrac particles 
       					#    advanced at system step 
double min_dt		=	0.001;



########################################################
# Smoothing length settings
########################################################
int nbrcut_max		=	40;
int nbrcut_min 		= 	25;
double nbrcut_fac 	= 	0.02;
double min_h 		= 	1e-8;
double max_h 		=	10.0;	#Particle will detach from neighbors
       					#  if it grows larger than this


########################################################
# Gravity MAC, error tolerances & constants
########################################################
double Gnewt		=	3.93935e-4;
int do_Arel 		= 	1;
int do_DL 		= 	0;
int do_SPHDL	 	=	1;	# New criterion: correct gravity 
    					#     (body-body) within 2h enforced
double epsilon		=	1e-30;	# Softening
#double errtol 		= 	1.0e-3;
					#Note: AREL mac will choose the LEAST 
					#      stringent constraint of the two.
					#Watchout: errtol is rescaled by 
					#	total/sysradius^2! (not appropriate 
					#	for explosions)
double errtol 	     	=     	10.;
#double frac_tol	=	2.0e-8;
double frac_tol		=	1e-2;
double default_nterms 	=	10;



########################################################
# Misc. parameters
########################################################
double core_radius	=	1095.42;
double gamma		=	1.5;
int do_boundary		=	0;
double r_limit		=	1000.;
int do_mergeparticles	=	0;



########################################################
# Artificial viscosity parameters
########################################################
double visc_alpha	=	1.;
double visc_beta	=	2.;
double visc_epsilon	=	0.01;



########################################################
# Add an additional external Potential
########################################################
int AddExternalPhi	=	0;



##############################################################
# This block sets parameters to produce initial conditions
##############################################################
int do_static		=	0;
double static_viscos 	= 	0.15;	# Viscosity to damp sound waves 
       					# 	    on boundaries
int do_static_u		=	0;	# Keep temperature uniform
double static_u 	=	0.100128486;
int do_reset		=	0;	# Resets all velocities and time to 0
int do_recenterofmass	=	0;
int do_initial	  	=   	0;	# Add centrifugal force for 
    					# corotating objects? omega is 2*pi/T
int do_rotation		=	0;  
double omega		=	0.2113;
int do_removeangmom	=	0;
double domegadt		=	-7.1059004e-05;         # Driving for 1% per orbit,
double tstop_removeangmom=	59.471702;   # for 2 orbits. Don't forget the "-"!

# New Spline Kernel Coefficients from Chris
# if omitted, the default is the Monaghan kernel 
#int kernel_ncoef1 = 5;
#struct {double kernel_coef1;}[5] = {1.2798, 0.0, -3.25703558174757, 3.0425, -0.83492};
#int kernel_ncoef2 = 5;
#struct {double kernel_coef2;}[5] = {2.1164, -3.3596555, 1.797, -0.32368, 0.000272};

#From Gabe
#int kernel_ncoef1 = 5;
#struct {double kernel_coef1;}[5] = {1.27951495845619, 0.0, -3.25703558174757, 3.03498707759319, -0.83541795014743};
#int kernel_ncoef2 = 4;
#struct {double kernel_coef2;}[4] = {2.11493290860362, -3.34167180058971, 1.75547211913700, -0.30668472299652};


# SDF-EOH
