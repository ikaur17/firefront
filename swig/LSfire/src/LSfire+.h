#include <omp.h>
#include <math.h>

#include "lsm_grid.h"
#include "lsm_data_arrays.h"

#include "fireLib.h"

#include "parameters.h"

/* general input data */

unsigned int flag_sanitycheck ; // for comparison purposes
unsigned int flag_backward_compatibility ; 

unsigned int Nthreads_max, Nthreads_tot, Nthreads ;

double t_max ; // max simulation time
LSMLIB_REAL deltat ; // time step
unsigned int flag_deltat_const, flag_t_end ;

double h0 ; // smoothing length
double D0 ; // diffusion coefficient ( h0^2 = 4*D = 2*sigma^2/t )
int ntotal ; //number of burnt points from forefire;
double t_ff ; //present time step from forefire;

unsigned int flag_hotfront ;
double tc_h_up, tc_f_up, tc_up ; // time scale fuel (upwind)
double tc_h_dw, tc_f_dw, tc_dw ; // time scale fuel (downwind)
int Nneighbour ;

unsigned int flag_init_circle, flag_init_rectangle ;
int init_circle_N, init_rectangle_N ;
double init_circle_xc, init_circle_yc, init_circle_r ;
double init_rectangle_x0, init_rectangle_y0, init_rectangle_Lx, init_rectangle_Ly ;
int init_circle_inside, init_rectangle_inside ;

/* grid variables */

Grid *gridLS ;

int Nx, Ny, N ;       // number of internal cells (x-dir, y-dir, total) 
int Nx_g, Ny_g, N_g ; // number of cells including ghostboxes (x-dir, y-dir, total)

double x_0, x_1, y_0, y_1 ;  // coordinates of the domain
double x_0_g, x_1_g, y_0_g, y_1_g ; // coordinates of the domain including ghostboxes

double Lx, Ly, A ; // dimensions of the internal domain 

int i0_g, i1_g, j0_g, j1_g ;
int i0_f, i1_f, j0_f, j1_f ;
int i0_i, i1_i, j0_i, j1_i ;

double dx, dy ;  // cell sizes (x-direction, y-direction) 
unsigned int flag_d_center ; // offset per tenere conto che i vaalori di phi si riferiscono ai centroidi (non ai nodi)
double d_center_x, d_center_y ;

unsigned int flag_kernel ;
unsigned int flag_kernel_LSM, flag_kernel_Gauss, flag_kernel_Gauss_num ;
unsigned int  flag_kernel_Gauss_FF, flag_kernel_Lognorm_FF ;
unsigned int flag_kernel_Weibull, flag_kernel_Lognorm, flag_kernel_Mainardi ;
double dr, deltat_kernel, dCsi ; // kernel discretization
double weibull_lambda, weibull_h ;
double lognorm_mu, lognorm_s, lognorm_Fr, lognorm_mean_l ;
double Io, If ;
double rho_amb, T_amb, c_pg, g_acc ;

int Nquad ; // number of nodes of quadrature rule
double *x_q, *w_q ; // nodes and weights of quadrature rule
double *xf, *yf ;
double *fhistory ; // coordinates of the firefront output

unsigned int flag_adiff ; // calculations for anomalous diffusion 
double Mbeta_beta, Mbeta_u, Mbeta_D ; // parameter for anomalous diffusion
unsigned int Mbeta_N ; // number of node for Mainardi function
double *x_Mbeta, *y_Mbeta ; // data of the Mainardi function

/* data structures and arrays: LSM */

LSM_DataArrays *dataLS ;

LSMLIB_REAL *ros ;   // ptr to rate-of-spread map
LSMLIB_REAL *nx ;
LSMLIB_REAL *ny ;

LSMLIB_REAL *phi, *phi_ff ;
LSMLIB_REAL *phieff ;
//LSMLIB_REAL *psi ;
double *psi ;
LSMLIB_REAL *tburn_phi ;
LSMLIB_REAL *tburn_phieff ;
LSMLIB_REAL *tburn_psi ;
LSMLIB_REAL *kernel, *kernel_old ; //AM, 9/5/2013
LSMLIB_REAL *r, *ri, *Dr, *DR, *ka, *kb, *kc, *th_ka, *th_kb, *th_kc, *th_m, *th_y ; 

/* data structures and arrays: FireLib */

FuelCatalogPtr catalog ; // fuel catalog handle

size_t Model ;  // NFFL

double LowHeatComb; // fuel low heat of combustion
double omega0 ;     // oven-dry mass of fuel consumed per unit area
double ros_byram ;  // ros given by Byram formula (given as input data)

double WindSpd ;    // wind speed
double WindDir ;    // degrees clockwise from north
double Slope ;      // fraction rise / reach
double Aspect ;     // degrees clockwise from north
double M1 ;         // 1-hr dead fuel moisture
double M10 ;        // 10-hr dead fuel moisture
double M100 ;       // 100-hr dead fuel moisture
double Mherb ;      // Live herbaceous fuel moisture
double Mwood ;      // Live woody fuel moisture

unsigned int flag_wind_const ; // wind as a function of time 

size_t *fuelMap ;       // ptr to fuel model map
LSMLIB_REAL *ignMap ;   // ptr to ignition time map (minutes)
//LSMLIB_REAL *flameMap ; // ptr to flame length map (feet)
LSMLIB_REAL *slpMap ;   // ptr to slope map (rise/reach)
LSMLIB_REAL *aspMap ;   // ptr to aspect map (degrees from north)
LSMLIB_REAL *wspdMap ;  // ptr to wind speed map (ft/min)
LSMLIB_REAL *wdirMap ;  // ptr to wind direction map (deg from north)
LSMLIB_REAL *m1Map ;    // ptr to 1-hr dead fuel moisture map
LSMLIB_REAL *m10Map ;   // ptr to 10-hr dead fuel moisture map
LSMLIB_REAL *m100Map ;  // ptr to 100-hr dead fuel moisture map
LSMLIB_REAL *mherbMap ; // ptr to live herbaceous fuel moisture map
LSMLIB_REAL *mwoodMap ; // ptr to live stem fuel moisture map

unsigned int flag_fireMap_const ;
unsigned int flag_byram, flag_ros_const , flag_update_wind ; 

unsigned int flag_fire_obstacle, flag_fire_obstacle_1, flag_fire_obstacle_2 ;
LSMLIB_REAL x0_fire_obstacle_1, x1_fire_obstacle_1 ;
LSMLIB_REAL y0_fire_obstacle_1, y1_fire_obstacle_1 ;
LSMLIB_REAL x0_fire_obstacle_2, x1_fire_obstacle_2 ;
LSMLIB_REAL y0_fire_obstacle_2, y1_fire_obstacle_2 ;

/* data saving */

unsigned int flag_saveformat_binary ;
unsigned int flag_saveformat_text ;
unsigned int flag_saveformat_ghostboxes ;

int save_each_iteration ;
char *filename_ws ;

unsigned int flag_save_data ; 
unsigned int flag_save_phi ;
unsigned int flag_save_phieff ;
unsigned int flag_save_psi ;
unsigned int flag_save_tburn ;
unsigned int flag_save_tburn_phi ;
unsigned int flag_save_tburn_phieff ;
unsigned int flag_save_tburn_psi ;
unsigned int flag_save_nx ;
unsigned int flag_save_ny ;
unsigned int flag_save_kernel ;
unsigned int flag_save_ros ;
unsigned int flag_save_fuelMap ;
unsigned int flag_save_windMap ;

char* filename_data ;
char* filename_phi ;
char* filename_phieff ;
char* filename_psi ;
char* filename_tburn_phi ;
char* filename_tburn_phieff ;
char* filename_tburn_psi ;
char* filename_nx ;
char* filename_ny ;
char* filename_kernel ;
char* filename_ros ;
char* filename_fuelMap ;
char* filename_windMap ;

char* filename_extension_txt ;
char* filename_extension_bin ;


//int save_to_file(char*, LSMLIB_REAL*, char*, char*) ;

void createCircle(LSMLIB_REAL*, LSMLIB_REAL, LSMLIB_REAL, LSMLIB_REAL,
	int, Grid*) ;
	
void createRectangle(
  LSMLIB_REAL *phi,
  LSMLIB_REAL corner_x, LSMLIB_REAL corner_y,
  LSMLIB_REAL side_length_x, LSMLIB_REAL side_length_y,
  int inside_flag,
  Grid *grid);
  
void signedLinearExtrapolationBC(LSMLIB_REAL*, Grid*, int) ;

int get_parameters(void) ;  
int get_input_fire(int, char**) ;
int get_input_adiff(int, char**) ;


int sanity_check(unsigned int, LSMLIB_REAL*, LSMLIB_REAL*, LSMLIB_REAL*, 
	LSMLIB_REAL*, LSMLIB_REAL*) ;

double fun_WindSpd(LSMLIB_REAL, LSMLIB_REAL, LSMLIB_REAL, LSMLIB_REAL) ;
double fun_WindDir(LSMLIB_REAL, LSMLIB_REAL, LSMLIB_REAL) ;
double fun_WindDir_der(LSMLIB_REAL, LSMLIB_REAL, LSMLIB_REAL) ;


void LSM2_signedLinearExtrapolationBC(double*, int*, int*, int*, int*,
    int*, int*, int*, int*) ;
                
void LSM2_LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS(double*, int*, int*, int*, int*, unsigned int*) ;

void LSM2_LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS(double*, int*, int*,
    int*, int*, double*, double*, int*, int*, int*, int*, double*,
    double*, int*, int*, int*, int*, double*, int*, int*, int*, int*,
    int*, int*, int*, int*, unsigned int*) ;
    
void LSM2_LSM2D_HJ_ENO1(double*, double*, int*, int*, int*, int*,
    double*, double*, int*, int*, int*, int*, double*, int*, int*, 
    int*, int*, double*, int*, int*, int*, int*, int*, int*, int*, int*,
    double*, double*, unsigned int*) ;
    
void LSM2_LSM2D_TVD_RK2_STAGE1(double*, int*, int*, int*, int*, double*, int*,
    int*, int*, int*, double*, int*, int*, int*, int*, int*, int*, int*,
    int*, double*, unsigned int*) ;
    
void LSM2_LSM2D_TVD_RK2_STAGE2(double*, int*, int*, int*, int*, double*,
    int*, int*, int*, int*, double*, int*, int*, int*, int*, double*,
    int*, int*, int*, int*, int*,  int*, int*, int*, double*, unsigned int*) ;
    
void LSM2_LSM2D_AVERAGE_GRAD_PHI(double*, double*, int*, int*, int*,
    int*, double*, double*, int*, int*, int*, int*, double*, double*,
    int*, int*, int*, int*, int*, int*, int*, int*, unsigned int*) ;
    
void LSM2_LSM2D_COMPUTE_UNIT_NORMAL(double*, double*, int*, int*, int*,
    int*, double*, double*, int*, int*, int*, int*, int*, int*, int*,
    int*, unsigned int*) ;
