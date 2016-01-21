//      LSFire+.c
//
//      Copyright 2011 CRS4  (LSFire.c)
//      Copyright 2013 CIRAM (LSFire+.c)
//
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <ctype.h>
#include <unistd.h>
     
#include "lsm_data_arrays.h"
#include "lsm_tvd_runge_kutta2d.h"
#include "lsm_level_set_evolution2d.h"
#include "lsm_spatial_derivatives2d.h"
#include "lsm_geometry2d.h"
// which, in turn, #include "lsm_grid.h", 
//                 #include "lsm_file.h"
// which, in turn, #include "LSMLIB_config.h"

#include "fireLib.h"
	
#include "LSfire+.h"

int i0_save, i1_save, j0_save, j1_save ;
int ic ;
int *nm, *nm_kernel ;
int Nr, NCsi ;
double dr, t_kernel ;
double sint, cost, Ut ;
double *Csi, *G1 ;
char filename_full[60] ;

double ros_max, ros_increment_max ;
char filename_xarray[60], filename_yarray[60], filename_fhistory[60] ;
double *xf, *yf;
double *xff, *yff;
double *fhistory ;
int ntotal  ;
int ib, j ;
double sum_quad;
int nvert;
#define DEBUG (0)

int pnpoly(int, double*, double*, double, double);
#define min(x,y) (x<y?x:y)
#define max(x,y) (x>y?x:y) 

#define deg2sin(angle) ( sin(angle*DEG2RAD) ) // degrees clockwise from north
#define deg2cos(angle) ( cos(angle*DEG2RAD) )


/***********************************************************************
 * Set the flag of the type of Kernel                                  *
 **********************************************************************/
int get_flag_kernel(unsigned int flag_kernel)
{
	flag_kernel_LSM       = 0 ;
	flag_kernel_Gauss     = 0 ;
	flag_kernel_Gauss_num = 0 ;
	flag_kernel_Weibull   = 0 ;
	flag_kernel_Lognorm   = 0 ;
	flag_kernel_Mainardi  = 0 ;
	flag_kernel_Gauss_FF  = 0 ; 
	flag_kernel_Lognorm_FF  = 0 ; 
	
	if (flag_kernel == 0)
		flag_kernel_LSM = 1 ;
		
	else if (flag_kernel == 10) {
		flag_kernel_Gauss = 1 ;
		}
	
	else if (flag_kernel == 11)
		flag_kernel_Gauss_num = 1 ;
	
	else if (flag_kernel == 20)
		flag_kernel_Weibull = 1 ;
	
	else if (flag_kernel == 30 || flag_kernel == 31 || flag_kernel == 32)
		flag_kernel_Lognorm = 1 ;
	
	else if (flag_kernel == 50)
        	flag_kernel_Gauss_FF = 1 ;

        else if (flag_kernel == 60 || flag_kernel == 61 || flag_kernel == 62)
                flag_kernel_Lognorm_FF = 1;

	else if (flag_kernel == 40)
		flag_kernel_Mainardi = 1 ;
        
    return 0 ;


	
} /*** end of: get_flag_kernel() **************************************/


/**********************************************************************
 * calculate lognorm parameters for firespotting                      *
 **********************************************************************/
int get_lognorm_par(double If, double Ut)
{
	// Ref.: N.Sardoy et al., Combustion and Flame 154 (2008), 478-488
	
	// If, fire intensity [kW/m]
	// Ut, wind speed [m/s]
	
	double Iff = If / 1e3 ; // [kW/m] -> [MW/m]
    double Lc ;
    
	// Calculate Fr
	
	if (flag_kernel == 30)
	{
		Lc = pow( If/(rho_amb*c_pg*T_amb*pow(g_acc, 0.5)), 2.0/3.0) ; // here we need If, not Iff
        lognorm_Fr = Ut / pow(g_acc*Lc, 0.5) ; // [non-dimensional]
    }
    else if (flag_kernel == 31)
    {
		lognorm_Fr = 0.0 ;
	}
	else if (flag_kernel == 32)
	{
		lognorm_Fr = INFINITY ;
	}
		
	/* Calculate mu, s, and <l> */
	
	if (lognorm_Fr < 1)// buoyancy driven regime
	{
        // note that Iff is the intensity in [MW/m]
		lognorm_mu = 1.47 * pow(Iff, 0.54) * pow(Ut, -0.55) + 1.14 ;
		lognorm_s  = 0.86 * pow(Iff, -0.21) * pow(Ut, 0.44) + 0.19 ;
	}
	else // lognorm_Fr >= 1 , wind-driven regime
	{
		lognorm_mu = 1.32 * pow(Iff, 0.26) * pow(Ut, 0.11) - 0.02 ;
		lognorm_s  = 4.95 * pow(Iff, -0.01) * pow(Ut, -0.02) - 3.48 ;
	}
    
    /* Note: lognorm_mu, lognorm_s, lognorm_mean_l as calculated above
     * require that l is given in [m]. Since we need to give the
     * argument of the pdf in [ft], these parameters have to be
     * suitably adjusted.
     * 
     * A translation of lognorm_mu is necessary:
     *     mu = <ln(l)>
     *        = <ln(l*M2FT/M2FT)> 
     *        = <ln(l*M2FT)> - ln(M2FT)
     * and so:
     *     lognorm_mu -> lognorm_u + ln(M2FT)
     * 
     * Being lognorm_s connected to the standard deviation (it is a 
     * shape factor), no translation/scaling is required.
     * 
     * For what concerns <l>, it can be computed by means of the 
     * formula (see wikipedia):
     *     <l> = exp(mu+0.5*sigma^2)
     * If we compute <l> before translating mu, the result will be
     * in [m], and we will have to multiply for the conversion factor:
     *     lognorm_mean_l ->  lognorm_mean_l * M2FT
     * The alternative is to compute <l> after the translation of mu.
     * In this case the result will be automatically given in [ft].
     * --- Undy, 19.10.2014 ---
     */
    lognorm_mu += LN_M2FT ;
    
    lognorm_mean_l = exp(lognorm_mu + 0.5*lognorm_s*lognorm_s) ;   
    
    
    if (DEBUG) 
	    printf("*** Ut = %f, Lc = %f > Fr = %f > mu = %f [->%f], s = %f, <l> = %f\n", 
		Ut, Lc, lognorm_Fr, lognorm_mu, lognorm_mu+LN_M2FT, lognorm_s, lognorm_mean_l) ;
	
	return 0 ;
	
} /*** end of: get_lognorm_par() ****************************************/



/**********************************************************************
 * preprocessing of the input data                                    *
 **********************************************************************/
int preprocessing()
{
    
	get_flag_kernel(flag_kernel) ;
	
	// define ther presence of some obstacles	
	if ( flag_fire_obstacle < 3 && flag_fire_obstacle_1 == 0 && flag_fire_obstacle_2 == 0 )
		flag_fire_obstacle = 0 ;
	else if ( flag_fire_obstacle < 3 && flag_fire_obstacle_1 == 1 && flag_fire_obstacle_2 == 0 )
		flag_fire_obstacle = 1 ;
	else if ( flag_fire_obstacle < 3 && flag_fire_obstacle_1 == 1 && flag_fire_obstacle_2 == 1 )
		flag_fire_obstacle = 2 ;
	else
		flag_fire_obstacle = 3 ;
		
        
    // check if the wind direction changes
    flag_update_wind = (!flag_fireMap_const) || (!flag_wind_const) ;
    
		
	// check if necessary to save time
	flag_save_tburn = flag_save_tburn_phi || flag_save_tburn_phieff || flag_save_tburn_psi ;
	
	
	// backward compatibility with LSfire
	if ( flag_backward_compatibility )
		ic = 0 ;
	else
		ic = 1 ;
        
    if ( flag_backward_compatibility && flag_d_center!=0 )	
	{
		printf("*** WARNING: 'flag_d_center' automatically modified [%d->0] ***\n", flag_d_center) ;
		flag_d_center = 0 ;
	}
    
    if (flag_d_center) {
        d_center_x = 0.5*dx ;
        d_center_y = 0.5*dy ;
    }
    else {
        d_center_x = 0.0 ;
        d_center_y = 0.0 ;
	}
    	
	if ( flag_backward_compatibility && flag_byram!=0 )	
	{
		printf("*** WARNING: 'flag_byram' automatically modified [%d->-1] ***\n", flag_byram) ;
		flag_byram = 0 ;
	}
	
	if (flag_sanitycheck || flag_backward_compatibility)
	{
		if (flag_t_end != -1)
		{
			printf("*** WARNING: 'flag_t_end' automatically modified [%d->-1] ***\n", flag_t_end) ;
			flag_t_end = -1 ;
		}
	}	
	
	if ( flag_backward_compatibility && Nneighbour!=3 )
	{
		printf("*** WARNING: 'Nneighbour' automatically modified [%d->-1] ***\n", Nneighbour) ;
		Nneighbour = 3 ;
	}
    
    /* preprocessing */
	
//	if (flag_kernel_Lognorm )
	if (flag_kernel_Lognorm || flag_kernel_Lognorm_FF)
	{
        double Ut = WindSpd * FTMIN2MS ;
		get_lognorm_par(If, Ut) ;
        
        printf("*** Lognorm parameters: If = %f > Ut = %f > Fr = %f > mu = %f, s = %f, <l> = %f\n", 
		    If, Ut, lognorm_Fr, lognorm_mu, lognorm_s, lognorm_mean_l) ;
	}
	else
	{
		lognorm_mu     = 0.0 ;
		lognorm_s      = 0.0 ;
		lognorm_Fr     = 0.0 ;
        lognorm_mean_l = 0.0 ;
	}
    
    return 0 ;
	
} /*** end of: preprocessing() ****************************************/


/**********************************************************************
 * Initialization of the grid                                         *
 **********************************************************************/
int init_grid()
{    
	LSMLIB_REAL x_lo[2] ;
	LSMLIB_REAL x_hi[2] ;
	
	x_1 = x_0 + Nx * dx ;
	y_1 = y_0 + Ny * dy ;
		
	x_lo[0] = x_0 ; x_lo[1] = y_0 ;
	x_hi[0] = x_1 ; x_hi[1] = y_1 ;
	    
    /* setup grid and data structure */   

    gridLS = createGridSetDx(2, dx, x_lo, x_hi, 1) ;
    
    
    if (Nx != gridLS->grid_dims[0])
    {
		printf("*** WARNING: Nx ricalcolato [%d --> %d] ***\n", Nx, gridLS->grid_dims[0]) ;
		Nx = gridLS->grid_dims[0] ;
	}
	
    if (Ny != gridLS->grid_dims[1]) 
    {
		printf("*** WARNING: Ny ricalcolato [%d --> %d] ***\n", Ny, gridLS->grid_dims[1]) ;
		Ny = gridLS->grid_dims[1] ;
	}
	
	if (x_1 != gridLS->x_hi[0])
    {
		printf("*** WARNING: x_1 ricalcolato [%f --> %f] ***\n", x_1, gridLS->x_hi[0]) ;
		x_1 = gridLS->x_hi[0] ;
	}
	
	if (y_1 != gridLS->x_hi[1])
    {
		printf("*** WARNING: y_1 ricalcolato [%f --> %f] ***\n", y_1, gridLS->x_hi[1]) ;
		y_1 = gridLS->x_hi[1] ;
	}
	
	N = Nx * Ny ;
	
	Lx = x_1 - x_0 ;
	Ly = y_1 - y_0 ;
	A = Lx*Ly ;
		    
    Nx_g = gridLS->grid_dims_ghostbox[0] ;
    Ny_g = gridLS->grid_dims_ghostbox[1] ;
    N_g  = gridLS->num_gridpts ;
    
    x_0_g = gridLS->x_lo_ghostbox[0] ; x_1_g = gridLS->x_hi_ghostbox[0] ;
    y_0_g = gridLS->x_lo_ghostbox[1] ; y_1_g = gridLS->x_hi_ghostbox[1] ;
    
    // gridLS->grid_dims[0]            Nx
    // gridLS->grid_dims[1]            Ny
    // gridLS->grid_dims_ghostbox[0]   Nx_g
    // gridLS->grid_dims_ghostbox[1]   Ny_g
    // gridLS->num_gridpts             Nx_g*Ny_g
    
    // indixes of the extended domain (including ghostboxes)
    i0_g = gridLS->ilo_gb ; i1_g = gridLS->ihi_gb ;
    j0_g = gridLS->jlo_gb ; j1_g = gridLS->jhi_gb ; 
	
	// indexes of the domain for the calculation of derivatives (fillboxes)
	i0_f = gridLS->ilo_fb ; i1_f = gridLS->ihi_fb ; 
	j0_f = gridLS->jlo_fb ; j1_f = gridLS->jhi_fb ; 
	
	// indexes of the internal domain
	i0_i = i0_g + (Nx_g-Nx)/2 ; i1_i = i1_g - (Nx_g-Nx)/2 ;
    	j0_i = j0_g + (Ny_g-Ny)/2 ; j1_i = j1_g - (Ny_g-Ny)/2 ;
    
	
	/* preprocessing */
	
	// calculate appropriate indeces of saved data
    
    if (flag_saveformat_ghostboxes)
    {
		i0_save = i0_g ; i1_save = i1_g ; 
		j0_save = j0_g ; j1_save = j1_g ;
	}
	else
	{
		i0_save = i0_i ; i1_save = i1_i ; 
		j0_save = j0_i ; j1_save = j1_i ;
	}
	
	// kernel calculation
    if (flag_kernel_Gauss_num)
		Nr = (int) sqrt((x_1-x_0)*(x_1-x_0) + (y_1-y_0)*(y_1-y_0)) / dr ;
	else
		Nr = 0 ;
			
	nm = (int*) calloc(2, sizeof(int)) ;
	
	nm[0] = j1_save - j0_save + 1 ;
	nm[1] = i1_save - i0_save + 1 ;	
	
	nm_kernel = (int*) calloc(2, sizeof(int)) ;
	
	nm_kernel[0] = 2 ;
	nm_kernel[1] = Nr ;
    
    return 0 ;
	
} /*** end of: init_grid() ********************************************/



/**********************************************************************
 * Initialization of the LSM data structures                          *
 **********************************************************************/
int init_lsm()
{    
    dataLS = allocateLSMDataArrays() ;
    
    allocateMemoryForLSMDataArrays(dataLS, gridLS) ;
        
    if (flag_init_circle)
        createCircle(dataLS->phi, init_circle_xc, init_circle_yc, 
			init_circle_r, init_circle_inside, gridLS) ;
            
	if (flag_init_rectangle)
		createRectangle(dataLS->phi, init_rectangle_x0, init_rectangle_y0,
			init_rectangle_Lx, init_rectangle_Ly, 
            init_rectangle_inside, gridLS) ;
    
    for (int idx = 0; idx < N_g; idx++)
    {
		//if (dataLS->phi[idx] < 0.) printf("%f, ",dataLS->phi[idx]);
		//if (dataLS->phi[idx] > 0.5) printf("%f, ",dataLS->phi[idx]);
		
		if (dataLS->phi[idx] < 0.) {
			dataLS->phi[idx] = 0. ;
        }
		else {
			dataLS->phi[idx] = 1.;
        }
			
		//dataLS->phi[n] = 1.0 - dataLS->phi[n] ;
	}
	    
    	nx     = (LSMLIB_REAL*) calloc(N_g, sizeof(LSMLIB_REAL)) ;
	ny     = (LSMLIB_REAL*) calloc(N_g, sizeof(LSMLIB_REAL)) ;
	phieff = (LSMLIB_REAL*) calloc(N_g, sizeof(LSMLIB_REAL)) ;
	phi_ff = (LSMLIB_REAL*) calloc(N_g, sizeof(LSMLIB_REAL)) ;
	psi    = (LSMLIB_REAL*) calloc(N_g, sizeof(LSMLIB_REAL)) ;
	//psi    = (double*) calloc(N_g, sizeof(double)) ;


	if (flag_save_tburn_phi)
		tburn_phi    = (LSMLIB_REAL*) calloc(N_g, sizeof(LSMLIB_REAL)) ;
	if (flag_save_tburn_phieff)
		tburn_phieff = (LSMLIB_REAL*) calloc(N_g, sizeof(LSMLIB_REAL)) ;
	if (flag_save_tburn_psi)
		tburn_psi    = (LSMLIB_REAL*) calloc(N_g, sizeof(LSMLIB_REAL)) ;
	
	
	for (int idx = 0; idx < N_g; idx++)
	{
		phieff[idx] = dataLS->phi[idx] ;
		phi_ff[idx] = 0.0;
		psi[idx]    = 0.0 ;
		if (flag_save_tburn_phi)    tburn_phi[idx]    = T_NO_BURN ;
		if (flag_save_tburn_phieff) tburn_phieff[idx] = T_NO_BURN ;
		if (flag_save_tburn_psi)    tburn_psi[idx]    = T_NO_BURN ;
	}
	

	
	if (flag_kernel_Gauss_num)
	{
		kernel     = (LSMLIB_REAL*) calloc(Nr, sizeof(LSMLIB_REAL)) ;
		kernel_old = (LSMLIB_REAL*) calloc(Nr, sizeof(LSMLIB_REAL)) ;
			
		r  = (LSMLIB_REAL*) calloc(Nr, sizeof(LSMLIB_REAL)) ; 
		ri = (LSMLIB_REAL*) calloc(Nr-1, sizeof(LSMLIB_REAL)) ;
		Dr = (LSMLIB_REAL*) calloc(Nr-1, sizeof(LSMLIB_REAL)) ;
		DR = (LSMLIB_REAL*) calloc(Nr, sizeof(LSMLIB_REAL)) ;
		ka = (LSMLIB_REAL*) calloc(Nr, sizeof(LSMLIB_REAL)) ;
		kb = (LSMLIB_REAL*) calloc(Nr-1, sizeof(LSMLIB_REAL)) ;
		kc = (LSMLIB_REAL*) calloc(Nr-1, sizeof(LSMLIB_REAL)) ;
		th_m = (LSMLIB_REAL*) calloc(Nr, sizeof(LSMLIB_REAL)) ;
		th_y = (LSMLIB_REAL*) calloc(Nr, sizeof(LSMLIB_REAL)) ; 
		th_ka = (LSMLIB_REAL*) calloc(Nr, sizeof(LSMLIB_REAL)) ;
		th_kb = (LSMLIB_REAL*) calloc(Nr, sizeof(LSMLIB_REAL)) ; 
		th_kc = (LSMLIB_REAL*) calloc(Nr, sizeof(LSMLIB_REAL)) ; 
	}
	else
	{
		flag_save_kernel = 0 ;
	}
    
    return 0 ;
    
} /*** end of: init_lsm() *********************************************/



/**********************************************************************
 * Initialization of the obstacle(s)                                  *
 **********************************************************************/
int init_obstacle()
{          
    if ( flag_fire_obstacle == 1 || flag_fire_obstacle == 2 )
    {
		//#pragma omp parallel for private(x0_fire_obstacle, \
        //    x1_fire_obstacle, y0_fire_obstacle, y1_fire_obstacle)
		for (int j = j0_g; j < j1_g + ic; j++)
        {
			for (int i = i0_g; i < i1_g + ic; i++)
            {
				double x = x_0_g + i * dx + d_center_x  ;
				double y = y_0_g + j * dy + d_center_y  ;
				
				if (x1_fire_obstacle_1 < x0_fire_obstacle_1) 
				{
					x0_fire_obstacle_1 = -1.0*INFINITY;
					x1_fire_obstacle_1 = INFINITY ;
				}
				if (y1_fire_obstacle_1 < y0_fire_obstacle_1) 
				{
					y0_fire_obstacle_1 = -1.0*INFINITY ;
					y1_fire_obstacle_1 = INFINITY ;
				}
				
				if ( (y>y0_fire_obstacle_1) && (y<y1_fire_obstacle_1) &&
					 (x>x0_fire_obstacle_1) && (x<x1_fire_obstacle_1) )
                {
                    int idx = i + j * Nx_g ;
//		    printf("index for obstacle %f %f %d\n", x*FT2M, y*FT2M, idx);
                    fuelMap[idx]  = 0 ;
                }
			}
		}
	}
	
	
	if ( flag_fire_obstacle == 2 )
    {
		for (int j = j0_g; j < j1_g + ic; j++) 
        {
			for (int i = i0_g; i < i1_g + ic; i++)
            {
				double x = x_0_g + i * dx + d_center_x  ;
				double y = y_0_g + j * dy + d_center_y  ;
				
				if (x1_fire_obstacle_2 < x0_fire_obstacle_2) 
				{
					x0_fire_obstacle_2 = -1.0*INFINITY;
					x1_fire_obstacle_2 = INFINITY ;
				}
				if (y1_fire_obstacle_2 < y0_fire_obstacle_2) 
				{
					y0_fire_obstacle_2 = -1.0*INFINITY ;
					y1_fire_obstacle_2 = INFINITY ;
				}
				
				if ( (y>y0_fire_obstacle_2) && (y<y1_fire_obstacle_2) &&
					 (x>x0_fire_obstacle_2) && (x<x1_fire_obstacle_2) )
                {
                    int idx = i + j * Nx_g ;
                    fuelMap[idx]  = 0 ;
                }
			}
		}
	}
	
	if ( flag_fire_obstacle == 3)
	{
		double x1, x2, y1, y2, y3, x1l, x1r, x2l, x2r, y3b, y3t, tga ;
		double dx1, dx2, dy1 ;
		
		x1 = 1000.0 ; dx1 = 400.0 ; 
		x2 = 3500.0 ; dx2 = 600.0 ; dy1 = 200.0 ; 
		y1 = 3000.0 ; y2 = 4000.0 ; y3 = 6000.0 ;
		 
		x1l = x1 - dx1/2.0 ; x1r = x1 + dx1/2.0 ;
		x2l = x2 - dx2/2.0 ; x2r = x2 + dx2/2.0 ;
		y3b = y3 - dy1/2.0 ; y3t = y3 + dy1/2.0 ;
		tga = tan(M_PI/3.0) ;
		
		for (int j = j0_g; j < j1_g + ic; j++)
		{
			for (int i = i0_g; i < i1_g + ic; i++) 
			{
				double x = x_0_g + i * dx + d_center_x  ;
				double y = y_0_g + j * dy + d_center_y  ;
				
				if ( (x>=x1l && x<=x1r && y>=y2 && y<=y3t) ||
                     (x>=x1l && x<=x2r && y>=y3b && y<= y3t) ||
                     (x>=x2l && x<=x2r && y>=y1 && y<=y3t) )
                {
                    int idx = i + j * Nx_g ;
                    fuelMap[idx] = 0 ;
                }
			}
		}		
	}
    
    return 0 ;
    
} /*** end of: init_obstacle() ****************************************/



/**********************************************************************
 * Initialization of the terrain and wind data arrays                 *
 **********************************************************************/
int init_fire() {

    int N_map ;
       
    catalog = Fire_FuelCatalogCreateStandard("Standard", NMAX_CATALOG) ;
    
    
    // evaluate memory requirement for terrain and wind data arrays,
    // and rate of spread...
    
    if (flag_fireMap_const) 
		N_map = 1 ;
	else 
		N_map = N_g ;
    

    // ...allocate memory...
    
	ignMap   = (LSMLIB_REAL*) calloc(N_map, sizeof(LSMLIB_REAL)) ;
	//flameMap = (LSMLIB_REAL*) calloc(N_map, sizeof(LSMLIB_REAL)) ;
	slpMap   = (LSMLIB_REAL*) calloc(N_map, sizeof(LSMLIB_REAL)) ;
	aspMap   = (LSMLIB_REAL*) calloc(N_map, sizeof(LSMLIB_REAL)) ;
	wspdMap  = (LSMLIB_REAL*) calloc(N_map, sizeof(LSMLIB_REAL)) ;
	wdirMap  = (LSMLIB_REAL*) calloc(N_map, sizeof(LSMLIB_REAL)) ;
	m1Map    = (LSMLIB_REAL*) calloc(N_map, sizeof(LSMLIB_REAL)) ;
	m10Map   = (LSMLIB_REAL*) calloc(N_map, sizeof(LSMLIB_REAL)) ;
	m100Map  = (LSMLIB_REAL*) calloc(N_map, sizeof(LSMLIB_REAL)) ;
	mherbMap = (LSMLIB_REAL*) calloc(N_map, sizeof(LSMLIB_REAL)) ;
	mwoodMap = (LSMLIB_REAL*) calloc(N_map, sizeof(LSMLIB_REAL)) ;
	
	ros = (LSMLIB_REAL*) calloc(N_g, sizeof(LSMLIB_REAL) ) ;
		
        
	// ...and initialize.

	
	if (flag_fireMap_const) 
	{	
		ignMap[0]   = INFINITY ;
		//flameMap[0] = 0. ;
		slpMap[0]   = Slope ;
		aspMap[0]   = Aspect  ; // [degrees]
        if (flag_wind_const) 
        {
            wspdMap[0]  = WindSpd ; // [ft/min]
            wdirMap[0]  = WindDir ; // [degrees]
        }
        else 
        {
            wspdMap[0]  = fun_WindSpd(0.0, 0.0, 0.0, WindSpd) ; // [ft/min]
            wdirMap[0]  = fun_WindDir(0.0, 0.0, 0.0) ; // [degrees]
        }
		m1Map[0]    = M1 ;
		m10Map[0]   = M10 ;
		m100Map[0]  = M100 ;
		mherbMap[0] = Mherb ;
		mwoodMap[0] = Mwood ; 
		
		
		//sint = sin(wdirMap[0]*DEG2RAD) ; //DegreesToRadians(WindDir)
		//cost = cos(wdirMap[0]*DEG2RAD) ; //DegreesToRadians(WindDir)
        //Ut = WindSpd * FTMIN2MS ;
        
        sint = deg2sin(wdirMap[0]) ;
        cost = deg2cos(wdirMap[0]) ;
		Ut = wspdMap[0] * FTMIN2MS ;
		
	} 
	else 
	{
		//#pragma omp parallel for
		for (int j = j0_g; j < j1_g + ic; j++) 
		{
			for (int i = i0_g; i < i1_g + ic; i++) 
			{
				int idx = i + j * Nx_g ;
                
                double x = x_0_g + i * dx + d_center_x  ;
				double y = y_0_g + j * dy + d_center_y  ;
            
				ignMap[idx]   = INFINITY ;
				//flameMap[idx] = 0. ;
				slpMap[idx]   = Slope ;
				aspMap[idx]   = Aspect  ;
                if (flag_wind_const)
                {
                    wspdMap[idx]  = WindSpd ;
                    wdirMap[idx]  = WindDir ;
                }
                else 
                {
                    wspdMap[idx]  = fun_WindSpd(0.0, x, y, WindSpd) ;
                    wdirMap[idx]  = fun_WindDir(0.0, x, y) ;
                }
				m1Map[idx]    = M1 ;
				m10Map[idx]   = M10 ;
				m100Map[idx]  = M100 ;
				mherbMap[idx] = Mherb ;
				mwoodMap[idx] = Mwood ;
				
			}
		}
	}

	
	// evaluate memory requirement for fuelMap, allocate, and initialize
	
/*	if (flag_fireMap_const && !flag_fire_obstacle && !flag_save_fuelMap)
	{
		printf("assigning size to fuelMap\n");
		fuelMap  = (size_t*) calloc(1, sizeof(size_t)) ;
		
		fuelMap[0]  = Model ;
	}		
	else
	{*/
		fuelMap  = (size_t*) calloc(N_g, sizeof(size_t)) ;
		//#pragma omp parallel for private(i, idx)
		for (int j = j0_g; j < j1_g + ic; j++) 
		{
			for (int i = i0_g; i < i1_g + ic; i++) 
			{
				int idx = i + j * Nx_g ;
				fuelMap[idx]  = Model ;
			}
		}
//	}
		
	flag_ros_const = 0 ;
			
	// setup of the obstacle	
	if (flag_fire_obstacle) init_obstacle() ;
	   
    return 0 ;
    
} /*** end of: init_fire() ********************************************/










/**********************************************************************
 * Initialization of nodes/weights of Gauss-Laguerre generalized rule *
 **********************************************************************/
int init_quad() {
	
	char filename_x[MAX_CHAR_NFILE], filename_w[MAX_CHAR_NFILE] ;
	FILE *fp_x, *fp_w ;
	
	
	if (flag_kernel_Lognorm || flag_kernel_Weibull || flag_kernel_Lognorm_FF) 
	{		
		x_q = (double*) calloc(Nquad, sizeof(double)) ;
		w_q = (double*) calloc(Nquad, sizeof(double)) ;
	}
	
	
	if (flag_kernel_Weibull) // generalized Gauss-Laguerre quadrature rule
	{
		sprintf(filename_x, "quadlib/quad_laguerre_%d_x.txt", Nquad) ;
		sprintf(filename_w, "quadlib/quad_laguerre_%d_w.txt", Nquad) ;
	}
	else if (flag_kernel_Lognorm || flag_kernel_Lognorm_FF) // generalized Gauss-Hermite quadrature rule
	{
		sprintf(filename_x, "quadlib/quad_hermite_%d_x.txt", Nquad) ;
		sprintf(filename_w, "quadlib/quad_hermite_%d_w.txt", Nquad) ;
	}
	
	
	fp_x = fopen(filename_x, "r") ;
	fp_w = fopen(filename_w, "r") ;
	
	
	if (flag_kernel_Weibull)
		printf("*** reading nodes/weights for Guass-Laguerre quadrature with %d nodes... ", Nquad) ;
	else if (flag_kernel_Lognorm || flag_kernel_Lognorm_FF)
		printf("*** reading nodes/weights for Guass-Hermite quadrature with %d nodes...", Nquad) ;

	
	for (int iq = 0; iq < Nquad; iq++)
	{
		fscanf(fp_x, "%lf", &x_q[iq]) ;
		fscanf(fp_w, "%lf", &w_q[iq]) ;
		if (DEBUG) printf("%02d: (%e, %e)\n", iq+1, x_q[iq], w_q[iq]) ;
	}
	
	printf("successfully read gauss-hermite quad coeff ***\n") ;
	
	fclose(fp_w) ;
	fclose(fp_x) ;
	
    return 0 ;

} /*** end of: init_quad() ********************************************/



/**********************************************************************
 * Initialization of data for ADIFF simulations                       *
 **********************************************************************/
int init_adiff() 
{
	if (flag_fire_obstacle_1 || flag_fire_obstacle_2 || flag_fire_obstacle)
	{
		flag_fire_obstacle_1 = flag_fire_obstacle_2 = flag_fire_obstacle = 0 ;
		printf("*** WARNING: adiff > all obstacles are getting ignored ***\n") ;
	}
	
	if (flag_save_psi)
	{
		flag_save_psi = 0 ;
		printf("*** WARNING: adiff > flag_save_psi set to 0 ***\n") ;
	}
	
	if (flag_save_tburn_psi)
	{
		flag_save_tburn_psi = 0 ;
		printf("*** WARNING: adiff > flag_save_tburn_psi set to 0 ***\n") ;
	}
	
	if (Mbeta_beta > 0 && !flag_kernel_Mainardi)
	{
		flag_kernel = 40 ;
		get_flag_kernel(flag_kernel) ;
		
		
		printf("*** WARNING: adiff > flag_kernel set to 40 (Mainardi function) ***\n") ;
	}
	else if (Mbeta_beta == 0)
	{
		flag_kernel = 0 ;
		get_flag_kernel(flag_kernel) ;
		printf("*** WARNING: adiff > flag_kernel set to 0 (pure LSM) ***\n") ;
	}
		
	flag_ros_const = 1 ;
	for (int j = j0_g; j < j1_g + ic; j++) 
    {
		for (int i = i0_g; i < i1_g + ic; i++)
        {
			int idx = i + j * Nx_g ;
            ros[idx]  = Mbeta_u ;
        }
    }
	
	return 0 ;
}



/**********************************************************************
 * Initialization of data for anomalous diffusion (Mainardi function) *
 **********************************************************************/
int init_anomalous() {
	
	char filename[MAX_CHAR_NFILE] ;
	double beta ;
	FILE *fp ;
	
	if (Mbeta_beta==0.)
		sprintf(filename, "../adlib/Mbeta_r_0.txt") ;
	else if (Mbeta_beta==1./4.)
		sprintf(filename, "../adlib/M_data_radial_beta_1_4.txt") ;
	else if (Mbeta_beta==1./2.)
		sprintf(filename, "../adlib/M_data_radial_beta_1_2.txt") ;
	else if (Mbeta_beta==3./4.)
		sprintf(filename, "../adlib/M_data_radial_beta_3_4.txt") ;
	else
	{
		printf("*** ERROR: no data available for Mbeta when beta = %f ***\n", Mbeta_beta) ;
		return -1 ;
	}
				
	fp = fopen(filename, "r") ;
	
	printf("*** reading data for Mainardi function with ") ;
	
	fscanf(fp, "%lf", &beta) ;
	fscanf(fp, "%d",  &Mbeta_N) ;
		
	printf("beta = %f (nodes: %d)... ", beta, Mbeta_N) ;
	
	x_Mbeta = (double*) calloc(Mbeta_N, sizeof(double)) ;
	y_Mbeta = (double*) calloc(Mbeta_N, sizeof(double)) ;
		
	for (int i = 0; i < Mbeta_N; i++)
	{
		fscanf(fp, "%lf %lf", &x_Mbeta[i], &y_Mbeta[i]) ;
		if (DEBUG) printf("%02d: (%22.15f, %18.15f)\n", i+1, x_Mbeta[i], y_Mbeta[i]) ;
	}
	
	fclose(fp) ;
	
	//printf(" ok ***\n") ;
	
	if (beta != Mbeta_beta) 
		printf("*** WARNING: read beta (%f) does not match input value (%f) ***\n", beta, Mbeta_beta) ;
	
	//M_interp_test(1) ;
	
	return 0 ;

} /*** end of: init_anomalous() ***************************************/



/***********************************************************************
 * Interpolation of r                                                  *
 **********************************************************************/
double M_interp(double rho, int verbosity)
{
	// returns: y_Mbeta[0]                    if rho<x_Mbeta[0]
	//          y_Mbeta[Mbeta_N-1]            if rho>x_Mbeta[Mbeta_N-1]
	//          interp(x_Mbeta, y_Mbeta, rho) otherwise
	//
	// other potential value to return (index):
	//			-1   if rho<x_Mbeta[0]
	//          -2   if rho>x_Mbeta[Mbeta_N-1]
	//           i   if M[i] <= rho < M[i+1]
	//
	// [implemented by Undy on October 22, 2013]
	
	unsigned int idx, idx0, idx1, i ;
	double y ;
	
	//if (verbosity) printf("\n") ;
	
	if (rho <= x_Mbeta[0])
	{
		idx = -1 ;
		y = y_Mbeta[0] ;		
		if (verbosity>0)
			printf("M_interp: *** below lower bound > idx=%d > r: %19.16f < %19.16f * y: %19.16f=%19.16f ***\n", 
			idx, rho, x_Mbeta[0], y, y_Mbeta[0]) ;
	}
	else if (rho >= x_Mbeta[Mbeta_N-1])
	{
		idx = -2 ;
		y = y_Mbeta[Mbeta_N-1] ;
		if (verbosity>0)
			printf("M_interp: *** above upper bound > idx=%d > r: %19.16f > %19.16f * y: %19.16f=%19.16f***\n", 
			idx, rho, x_Mbeta[Mbeta_N-1], y, y_Mbeta[Mbeta_N-1]) ;
	}
	else
	{
		idx0 = 0 ; idx1 = Mbeta_N - 1 ; i = 0 ;
		do
		{
			idx = (idx0 + idx1) / 2 ; // indice ultimo elemento I blocco
			
			if (verbosity>1)
				printf("search: [#%d] %d (%d) %d : [%f (%f) %f] searching for %f\n", 
				i, idx0, idx, idx1, x_Mbeta[idx0], x_Mbeta[idx], x_Mbeta[idx1], rho );
			
			if ( rho < x_Mbeta[idx] )
				idx1 = idx ;
			else
				idx0 = idx ;
				
			i++ ;
		
		} while ( idx1 > idx0+1 ) ;
		
		idx = idx0 ;
		y = y_Mbeta[idx] + 
            (y_Mbeta[idx+1]-y_Mbeta[idx]) * 
            (rho-x_Mbeta[idx])/(x_Mbeta[idx+1]-x_Mbeta[idx]) ;
		
		if (verbosity>2)
			printf("M_interp: idx=%d | r:%19.16f <=%19.16f <%19.16f | y:(%19.16f)%19.16f (%19.16f) [%d iter]\n", 
			idx0, x_Mbeta[idx], rho, x_Mbeta[idx+1], y_Mbeta[idx], y, y_Mbeta[idx+1], i) ;
	}
    
	return y ;
	
} /*** end of: M_interp() *********************************************/



/***********************************************************************
 * Small utility to test the interpolation function                    *
 **********************************************************************/
int M_interp_test(int verbosity)
{
	double ytest ;
    
    // just for testing purposes //
	
	ytest = M_interp(-0.01, verbosity) ;
	ytest = M_interp(0.0, verbosity) ;
	ytest = M_interp(0.000001, verbosity) ;
	ytest = M_interp(0.4100000000000003, verbosity) ;
	ytest = M_interp(1., verbosity) ;
	ytest = M_interp(2.5, verbosity) ;
	ytest = M_interp(M_PI, verbosity) ;
	ytest = M_interp(11.028000000000352, verbosity) ;
	ytest = M_interp(1e7, verbosity) ;
	
	printf("M_interp: *** y: %19.16f\n", ytest) ;
	
	return 0 ;
	
} /*** end of: M_interp_test() ****************************************/



/**********************************************************************
 * Evaluation of the rate of spread                                   *
 **********************************************************************/
int get_RoS(LSMLIB_REAL t)
{
    ros_max = 0.0 ;
    ros_increment_max = 0.0 ;
    
    //#pragma omp parallel for
    // THIS IS NOT THREAD-SAFE BECAUSE OF CALLS TO FIRELIB
	for (int j = j0_g; j < j1_g + ic; j++) 
	{
		for (int i = i0_g; i < i1_g + ic; i++) 
		{
			int idx = i + j * Nx_g ;
            
            double x = x_0_g + i * dx + d_center_x ;
            double y = y_0_g + j * dy + d_center_y ;
           
			int idxMap, idxMap_fuel ;
            
            if (flag_fireMap_const)
				idxMap = 0 ;
            else
                idxMap = idx ;
                
            if (flag_fireMap_const && !flag_fire_obstacle)
				idxMap_fuel = 0 ;
			else
				idxMap_fuel = idx ;
           
			size_t modelNumber = fuelMap[idxMap_fuel] ;
            
            double moisture[6] ; // fuel moisture content at current cell
		
            moisture[0] = m1Map[idxMap] ;
			moisture[1] = m10Map[idxMap] ;
			moisture[2] = m100Map[idxMap] ;
			moisture[3] = m100Map[idxMap] ;
			moisture[4] = mherbMap[idxMap] ;
			moisture[5] = mwoodMap[idxMap] ;
            
            
            double wspd, wdir, wdir_der ;
            
            if (flag_wind_const)
            {
                wspd = WindSpd ;
                wdir = WindDir ;
                wdir_der = 0.0 ;
            }
            else 
            {
                wspd     = fun_WindSpd(t, x, y, WindSpd) ;
                wdir     = fun_WindDir(t, x, y) ;
                wdir_der = fun_WindDir_der(t, x, y) ;
                
                double Ut = wspd * FTMIN2MS ;
                
                get_lognorm_par(If, Ut) ;
            }
            
                
			//slopeN = slpMap[idxMap] * 
            //    (deg2sin(aspMap[idxMap])*nx[idx]+deg2cos(aspMap[idxMap])*ny[idx]) ;
            //if (slopeN < 0.0) slopeN = 0.0 ;
                    
            double r0, rw, fw, rs, fs ; 
            double f_max = 1.0 ;

            
            if (flag_byram == 1)
			{
                double windN = wspd ;
                
                double theta = aspMap[idxMap]*DEG2RAD ;
				double slopeN = slpMap[idxMap] * (sin(theta)*nx[idx]+cos(theta)*ny[idx]) ;
				if (slopeN < 0.) slopeN = 0.0 ;
                
				Fire_SpreadNoWindNoSlope(catalog, modelNumber, moisture) ;
				r0 = Fuel_Spread0(catalog, modelNumber) ;
				
				Fire_SpreadWindSlopeMax(catalog, modelNumber, windN, 0., 0., 0.) ;
				rw = Fuel_SpreadMax(catalog, modelNumber) ;
				fw = (rw - r0) / r0 ;
				
				Fire_SpreadWindSlopeMax(catalog, modelNumber, 0., 0., slopeN, 0.) ;
				rs = Fuel_SpreadMax(catalog, modelNumber) ;
				fs = (rs - r0) / r0 ;
				
				f_max = 1.0 + fw + fs ;
			}
			            
            // windN is the wind component normal to the front
            
            
			double nS_dot_n = deg2sin(aspMap[idxMap])*nx[idx] + deg2cos(aspMap[idxMap])*ny[idx] ;
            double slopeN = slpMap[idxMap] * nS_dot_n ;
			if (slopeN < 0.) slopeN = 0.0 ;
            
            double nU_dot_n = deg2sin(wdir)*nx[idx] + deg2cos(wdir)*ny[idx] ;
            double windN = wspd * nU_dot_n ;
			if (windN < 0.) windN = 0.0 ;
		
			Fire_SpreadNoWindNoSlope(catalog, modelNumber, moisture) ;
			r0 = Fuel_Spread0(catalog, modelNumber) ;
			
			Fire_SpreadWindSlopeMax(catalog, modelNumber, windN, 0., 0., 0.) ;
			rw = Fuel_SpreadMax(catalog, modelNumber) ;
			fw = (rw - r0) / r0 ;
			
			Fire_SpreadWindSlopeMax(catalog, modelNumber, 0., 0., slopeN, 0.) ;
			rs = Fuel_SpreadMax(catalog, modelNumber) ;
			fs = (rs - r0) / r0 ;
			
			if (flag_byram == 0)
				ros[idx] = r0 * (1.0 + fw + fs) ;
			else
				ros[idx] = ros_byram * (1.0 + fw + fs) / f_max ;
                
            ros_max = max(ros_max, ros[idx]) ;
                
            /* additional terms (implemented on October 18, 2014): 
             * 
             * V_l_1 : d<l>/dt n_U * dn/dt = 0
             *         
             *         -> for the moment this is zero because <l>
             *            does not vary in time
             * 
             * V_l_2 : <l> * dn_U/dt.n =
             *         <l> * [d(sin(z*DEG2RAD))/dt, d(cos(z*DEG2RAD))/dt] =
             *         DEG2RAD * dz/dt * [cos(z*DEG2RAD), -sin(z*DEG2RAD)]
             * 
             *         -> z is the wind angle in degrees (clockwise 
             *            from north), provided by fun_WindDir()
             *         -> dz/dt is provided by fun_WindDir_der()
             */
             
            double V_l_1 = 0.0 ;
            
            double dn_U_dt_dot_n = DEG2RAD * wdir_der * 
                (deg2cos(wdir)*nx[idx] - deg2sin(wdir)*ny[idx] ) ;
                
            double V_l_2 = lognorm_mean_l * dn_U_dt_dot_n ;
            
            V_l_2 = max(0.0, V_l_2) ;
            
            double ros_old = ros[idx] ;
                                   
            ros[idx] += ( V_l_1 +  V_l_2 ) ;
            
            ros_increment_max = max(ros_increment_max, (ros[idx]/ros_old-1.0)) ;
// 3% rate of spread ------------------------------            
//	    ros[idx] = 0.03*wspd;
//	    if (fuelMap[idxMap_fuel] > 0)ros[idx] = 0.05*M2FT*60.0;
//-------------------------------------------------
		}
	}
    
    //printf("ros_increment_max at t=%f: %f\n", t, ros_increment_max) ;
    
    //printf("V_l_2: %f, %f (%f, %f)\n", V_l_2_min, V_l_2_max, ros_min, ros_max);
    
    return 0 ;

} /*** end of: get_RoS() **********************************************/






/**********************************************************************
 * Time iteration of the LSM method (RK2)                             *
 **********************************************************************/
int lsm2_RK2(LSMLIB_REAL dt, LSMLIB_REAL t) 
{
	//int iflagtest = 0 ;		
	
	for (int istage = 1; istage < 3; istage++)
	{			
		// set the rhs to zero
        LSM2_LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS(
            dataLS->lse_rhs, &i0_g, &i1_g, &j0_g, &j1_g, &Nthreads) ;

		// calculate the upwind approximation of Grad(phi)
        
		if (istage==1)
			LSM2_LSM2D_HJ_ENO1(
                dataLS->phi_x_plus, dataLS->phi_y_plus, &i0_g, 
				&i1_g, &j0_g, &j1_g, dataLS->phi_x_minus, 
				dataLS->phi_y_minus, &i0_g, &i1_g, &j0_g, &j1_g,
				dataLS->phi, &i0_g, &i1_g, &j0_g, &j1_g, dataLS->D1, 
				&i0_g, &i1_g, &j0_g, &j1_g, &i0_f, &i1_f, &j0_f, &j1_f,	
				&dx, &dy, &Nthreads) ;
		else
			LSM2_LSM2D_HJ_ENO1(
                dataLS->phi_x_plus, dataLS->phi_y_plus, &i0_g,
				&i1_g, &j0_g, &j1_g, dataLS->phi_x_minus, 
				dataLS->phi_y_minus, &i0_g, &i1_g, &j0_g, &j1_g,
				dataLS->phi_stage1,	&i0_g, &i1_g, &j0_g, &j1_g, 
				dataLS->D1,	&i0_g, &i1_g, &j0_g, &j1_g, &i0_f, &i1_f,
				&j0_f, &j1_f, &dx, &dy, &Nthreads) ;
			
		LSM2_LSM2D_AVERAGE_GRAD_PHI(
            dataLS->phi_x, dataLS->phi_y, &i0_g,
			&i1_g, &j0_g, &j1_g, dataLS->phi_x_plus, dataLS->phi_y_plus,
			&i0_g, &i1_g, &j0_g, &j1_g, dataLS->phi_x_minus, 
			dataLS->phi_y_minus, &i0_g, &i1_g, &j0_g, &j1_g, &i0_f, 
			&i1_f, &j0_f, &j1_f, &Nthreads) ;
		
		// compute the normal to the surface phi=0
		LSM2_LSM2D_COMPUTE_UNIT_NORMAL(
            nx, ny, &i0_g, &i1_g, &j0_g, &j1_g,
			dataLS->phi_x, dataLS->phi_y, &i0_g, &i1_g, &j0_g, &j1_g,
			&i0_f, &i1_f, &j0_f, &j1_f, &Nthreads) ;
		
		// compute the spread rate
        
		LSM2_signedLinearExtrapolationBC(
           nx, &i0_g, &i1_g, &j0_g, &j1_g, &i0_f, &i1_f, &j0_f, &j1_f) ;
		LSM2_signedLinearExtrapolationBC(
           ny, &i0_g, &i1_g, &j0_g, &j1_g, &i0_f, &i1_f, &j0_f, &j1_f) ;
		//signedLinearExtrapolationBC(ny, gridLS, 9) ;
		//signedLinearExtrapolationBC(ny, gridLS, 9) ;
       
		if (!flag_ros_const) get_RoS(t) ;
        
		// compute the rhs of the level set equation for a normal but 
		// non uniform propagation
		LSM2_LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS(
            dataLS->lse_rhs, &i0_g, 
			&i1_g, &j0_g, &j1_g, dataLS->phi_x_plus, dataLS->phi_y_plus,
			&i0_g, &i1_g, &j0_g, &j1_g, dataLS->phi_x_minus, 
			dataLS->phi_y_minus, &i0_g, &i1_g, &j0_g, &j1_g, ros, &i0_g,
			&i1_g, &j0_g, &j1_g, &i0_f, &i1_f, &j0_f, &j1_f, &Nthreads) ;
		
		//~ for (int idx = 0; idx < N_g; idx++)
			//~ dataLS->lse_rhs[idx] += 0.0 * (dataLS->phi_x_plus[idx] - 
				//~ dataLS->phi_x_minus[idx] + dataLS->phi_y_plus[idx] - 
				//~ dataLS->phi_y_minus[idx]) ;                               // why?
		        
		// advance phi and fill boundary cells
		if (istage==1)  
		{
			LSM2_LSM2D_TVD_RK2_STAGE1(dataLS->phi_stage1, &i0_g, &i1_g,
				&j0_g, &j1_g, dataLS->phi, &i0_g, &i1_g, &j0_g, &j1_g,
				dataLS->lse_rhs, &i0_g, &i1_g, &j0_g, &j1_g, &i0_f, 
				&i1_f, &j0_f, &j1_f, &dt, &Nthreads) ;
						
			LSM2_signedLinearExtrapolationBC(
                dataLS->phi_stage1, 
                &i0_g, &i1_g, &j0_g, &j1_g, &i0_f, &i1_f, &j0_f, &j1_f) ;
			//signedLinearExtrapolationBC(dataLS->phi_stage1, gridLS, 9) ;
		}
		else 
		{
			LSM2_LSM2D_TVD_RK2_STAGE2(
                dataLS->phi_next, &i0_g, &i1_g, &j0_g, 
				&j1_g, dataLS->phi_stage1, &i0_g, &i1_g, &j0_g, &j1_g,
				dataLS->phi, &i0_g, &i1_g, &j0_g, &j1_g, 
				dataLS->lse_rhs, &i0_g, &i1_g, &j0_g, &j1_g, &i0_f, 
				&i1_f, &j0_f, &j1_f, &dt, &Nthreads) ;
						
			LSM2_signedLinearExtrapolationBC(
                dataLS->phi_next, 
                &i0_g, &i1_g, &j0_g, &j1_g, &i0_f, &i1_f, &j0_f, &j1_f) ;
			//signedLinearExtrapolationBC(dataLS->phi_next, gridLS, 9) ;
		}       
			
	}
    
    return 0 ;
    
} /*** end of: lsm2_RK2() **********************************************/







/**********************************************************************
 * Evaluation of the kernel (Gaussian kernel; analytical solution)    *
 **********************************************************************/
double get_kernel_gauss_a(double r, double t) 
{
    unsigned int flag_original ;
    
    flag_original = 1 ;
    
    if (flag_original==1)
    {
		//if (h==0) return 6.366198e-03 ;
		if ( r < 3.0*h0*t )
			return exp(-(r*r)/(h0*h0*t)) / (3.14*h0*h0*t) ;
		else
			return 0.0 ;			
	}
	else
	{
		return exp(-(r*r)/(h0*h0*t)) / (M_PI*h0*h0*t) ;
	}
	
	return -1 ;
   
} /*** end of: get_kernel_gauss_a() ***********************************/



/**********************************************************************
 * Evaluation of the kernel (Gaussian kernel; numerical solution)     *
 **********************************************************************/
int get_kernel_gauss_n(double t_old, double t) 
{
    double last_grad, t_kernel, dt, D ;
    int i ;
    
    D = h0 * h0 / 4.0 ;
    
    // calculation of the coefficients 
    void coeff(int Nr, double dr)
    {
		int i ;
		
		for (i = 0; i < Nr; i++) 
			r[i] = i * dr ;
			
		for (i = 0; i < Nr-1; i++)
			Dr[i] = r[i+1] - r[i] ;
			
		for (i = 0; i < Nr-1; i++)
			ri[i] = 0.5 * (r[i+1] + r[i]) ;
			
		DR[0] = r[0] * Dr[0] ;
		for (i = 1; i < Nr-1; i++)
			DR[i] = r[i] * ( ri[i] - ri[i-1] ) ;
		DR[Nr-1] = r[Nr-1] * Dr[Nr-2] ; 
		
		for (i = 0; i < Nr-1; i++)
			kc[i] = ri[i] / Dr[i] / DR[i+1] ;
		
		kb[0] = 2.0 / (Dr[0]*Dr[0]) ;
		for (i = 1; i < Nr-1; i++)
			kb[i] =  ri[i] / Dr[i] / DR[i] ;
			
		ka[0] = -2.0 / (Dr[0]*Dr[0]) ;
		for (i = 1; i < Nr-1; i++)
			ka[i] = - ( kc[i-1] + kb[i] ) ;
		ka[Nr-1] = - ( 2.0*ri[Nr-2]/Dr[Nr-2] + 1.0 ) / DR[Nr-1] ; 
	
	}
	
	// solution of the tridiagonal system
	void thomas(int Nr)
	{
		int i ;
		
		th_m[0] = th_ka[0] ;
		for (i = 1; i < Nr; i++) th_m[i] = 0.0 ;
		
		th_y[0] = kernel_old[0] ;
		for (i = 1; i < Nr; i++) th_y[i] = 0.0 ;
		
		for (i = 1; i < Nr; i++)
		{
			double l = th_kc[i-1] / th_m[i-1] ;
			th_m[i] = th_ka[i] - l * th_kb[i-1] ;
			th_y[i] = kernel_old[i] - l * th_y[i-1] ;
		}
		
		kernel[Nr-1] = th_y[Nr-1] / th_m[Nr-1] ;
		for (i = Nr-2; i >= 0; i--)
			kernel[i] = ( th_y[i] - th_kb[i]*kernel[i+1] ) / th_m[i] ;
			
	}
	
	// initialization
	if (t == 0.0)
	{
		coeff(Nr, dr) ;
		kernel_old[0] = 1.0 / (M_PI*(r[1]/2.0)*(r[1]/2.0)) / 2.0 ;
		kernel[0] = kernel_old[0] ;
	}
	
	// time iteration 
	t_kernel = t_old ;
	
	while ( t_kernel < t )
	{
		t_kernel += deltat_kernel ;
		if (t_kernel > t)
		{
			dt = t - t_kernel + deltat_kernel ;
			t_kernel = t ;
		}
		else
		{
			dt = deltat_kernel ;
		}
		
		for (i = 0; i < Nr-1; i++)
		{
			th_ka[i] = 1.0 - dt * D * ka[i] ;
			th_kb[i] = - dt * D * kb[i] ;
			th_kc[i] = - dt * D * kc[i] ;
		}
		th_ka[Nr-1] = 1.0 - dt * D * ka[Nr-i] ;
			
		thomas(Nr) ;
		
		last_grad = fabs( (kernel[Nr-1]-kernel[Nr-2])/(r[Nr-1]-r[Nr-2]) ) ;
		
		if (last_grad > MAX_LAST_GRAD)
			printf("...WARNING: last_grad > MAX_LAST_GRAD (t_kernel=%f, t=%f)\n", 
				t_kernel, t) ;
		
		for (i = 0; i < Nr; i++)
			kernel_old[i] = kernel[i] ;
			
		printf(".") ;
				
	}

    return 0 ;
    
} /*** end of: get_kernel_gauss_n() ***********************************/



/**********************************************************************
 * Evaluation of phi                                                  *
 **********************************************************************/
 int update_phi()
{

	//#pragma omp parallel
    for (int idx = 0; idx < N_g; idx++) 
		dataLS->phi[idx] = dataLS->phi_next[idx] ;
	
	return 0 ;
	
} /*** end of: update_phi() **********************************************/	



/**********************************************************************
 * Evaluation of phieff (integral of K*phi)                           *
 **********************************************************************/
int get_phieff_LSM()
{   FILE *fphi_d, *fphi_c;
    fphi_d = fopen("phieff_slice.txt", "w+");
    fphi_c = fopen("phi_slice_continous.txt", "w+");
//    #pragma omp parallel

         for (int ja = j0_g; ja < j1_g + ic; ja++)
                 {
                for (int ia = i0_g; ia < i1_g + ic; ia++)
                        {
                                int idxa = ia + ja * Nx_g;
                                double xa = x_0_g + ia * dx + d_center_x ;
                                double ya = y_0_g + ja * dy + d_center_y ;
//                                if (dataLS->phi[idxa] <= 0.5 )phi_box[idxa]= 0.0;
//                                if (dataLS->phi[idxa] > 0.5) phi_box[idxa]= 1.0;
                                phieff[idxa] = dataLS->phi[idxa] ;
//                              phieff[idxa] = phi_box[idxa] ;
//                      if(ya*FT2M < 2505.0 && ya*FT2M  > 2495.0)fprintf(fphi_d,"%f %f %f %f\n",xa*FT2M,ya*FT2M,dataLS->phi[idxa],phi_box[idxa]);
//                        if(xa*FT2M < 3505.0 && xa*FT2M  > 3495.0)fprintf(fphi_d,"%f %f %f %f\n",xa*FT2M,ya*FT2M,dataLS->phi[idxa],phi_box[idxa]);
                        }
                }


        fclose(fphi_c);
        fclose(fphi_d);
        return 0 ;

} /*** end of: get_phieff_LSM() ***************************************/


/**********************************************************************
 * Evaluation of phieff (integral of K*phi)   for output from ForeFire                        *
 **********************************************************************/
int get_phieff_LSM_FF()
{       
    #pragma omp parallel
    for (int idx = 0; idx < N_g ; idx++)
    {
		//phieff[idx] = dataLS->phi[idx] ;
		phieff[idx] = phi_ff[idx] ;
	}
    
     FILE *fphi_d, *fphi_c;
    fphi_d = fopen("phieff_slice.txt", "w+");
//    #pragma omp parallel

         for (int ja = j0_g; ja < j1_g + ic; ja++)
                 {
                for (int ia = i0_g; ia < i1_g + ic; ia++)
                        {
                                int idxa = ia + ja * Nx_g;
                                double xa = x_0_g + ia * dx + d_center_x ;
                                double ya = y_0_g + ja * dy + d_center_y ;
                                phieff[idxa] = phi_ff[idxa] ;
                        //if(ya*FT2M < 2505.0 && ya*FT2M  > 2495.0)fprintf(fphi_d,"%f %f %f %f\n",xa*FT2M,ya*FT2M,phi_ff[idxa],phieff[idxa]);
                        if(xa*FT2M < 3505.0 && xa*FT2M  > 3495.0)fprintf(fphi_d,"%f %f %f %f\n",xa*FT2M,ya*FT2M,phi_ff[idxa],phieff[idxa]);
                        }
                }


        fclose(fphi_d);
	
	return 0 ;
    
} /*** end of: get_phieff_LSM(_FF) ***************************************/

/**********************************************************************
 * Evaluation of phieff (integral of K*phi)                           *
 **********************************************************************/
int get_phieff_Gauss(double t)
{       printf("Including the Turbulence and the time is = %f\n", t);
        //printf("h0,t,dx,dy %f %f %f %f\n", h0,t,dx,dy);
//        FILE *fp1, *fp2, *fp3, *fphi;
//        fp1 = fopen("phieff_0.2.txt", "w+");
//        fp2 = fopen("phieff_0.5.txt", "w+");
//        fp3 = fopen("phieff_0.9.txt", "w+");
//        fphi = fopen("phieff_slice.txt", "w+");
         for (int ja = j0_g; ja < j1_g + ic; ja++)
    {
                for (int ia = i0_g; ia < i1_g + ic; ia++)
        {
                int idxa = ia + ja * Nx_g ;
//                if (dataLS->phi[idxa] <= 0.5 ) phi_box[idxa]= 0.0;
//                if (dataLS->phi[idxa] > 0.5) phi_box[idxa]= 1.0;
        }
    }


    #pragma omp parallel for
    for (int ja = j0_g; ja < j1_g + ic; ja++)
    {
                for (int ia = i0_g; ia < i1_g + ic; ia++)
        {
            int idxa = ia + ja * Nx_g ;

            double xa = x_0_g + ia * dx + d_center_x ;
            double ya = y_0_g + ja * dy + d_center_y ;
            phieff[idxa] = 0.0 ;

                        double normalizzatore = 0.0 ;

                        int i0b = max(i0_g, ia-Nneighbour*(int)((h0*t)/dx)) ;
                        int i1b = min(i1_g, ia+Nneighbour*(int)((h0*t)/dx)) ;
                        int j0b = max(j0_g, ja-Nneighbour*(int)((h0*t)/dy)) ;
                        int j1b = min(j1_g, ja+Nneighbour*(int)((h0*t)/dy)) ;

            for (int jb = j0b; jb < j1b + ic; jb++)
            {
                for (int ib = i0b; ib < i1b + ic; ib++)
                {
                    int idxb = ib + jb * Nx_g ;

                    double xb = x_0_g + ib * dx + d_center_x ;
                    double yb = y_0_g + jb * dy + d_center_y ;
                    //double r = sqrt( (xa-xb)*(xa-xb)+(ya-yb)*(ya-yb) ) ;
                    //double r = sqrt( (ya-yb-4000)*(ya-yb-4000) ) ;
                                        double e = 0.0 ;

                                        if (flag_backward_compatibility)
                                        {
                                                if ( ((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya)) < (9.0*h0*h0*t) )
                                                        e = exp(-((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))/(h0*h0*t)) ;
                                                        //e = (exp(-((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))/(h0*h0*t)) + exp(-((xb-xa)*(xb-xa)+(yb-ya-500)*(yb-ya-500))/(h0*h0*t)));
                                                        //e = ( exp(-((xb-xa)*(xb-xa)+(yb-ya-500.)*(yb-ya-500.))/(h0*h0*t)));
                                        }
                                        else
                                        {
                                                e = exp(-((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))/(h0*h0*t)) ;
                                                //e =0.5*(exp(-((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))/(h0*h0*t))+exp(-((yb-ya+6640)*(yb-ya+6640))/(h0*h0*t))) ;
                                                //e =0.5*(exp(-((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))/(h0*h0*t))+exp(-((xb-xa)*(xb-xa)+(yb-ya+8440)*(yb-ya+8440))/(h0*h0*t))) ;
                                                //e =0.5*(exp(-((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))/(h0*h0*t))+exp(-((yb-ya+7440)*(yb-ya+7440))/(h0*h0*t))) ;
                                                //if (r>0)e =0.5*(exp(-((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))/(h0*h0*t)) + exp(-(log(r)-mu)*(log(r)-mu)/(2*si*si))/(sqrt(2*M_PI)*si*r)) ;
                                                //if (r>0)e = exp(-(log(r)-mu)*(log(r)-mu)/(2*si*si))/(sqrt(2*M_PI)*si*r) ;
                                        }
                                        phieff[idxa] += dataLS->phi[idxb] * e ;
                                        //phieff[idxa] += phi_box[idxb] * e ;   
                                        normalizzatore += e ;

                }
            }



            // NB: alpha = dx*dy/(M_PI*twosigma2) = dx*dy/(M_PI*h0*h0*t)
            if (flag_backward_compatibility)
            {
                                phieff[idxa] *= (dx*dy/(3.14*h0*h0*t)) ;
                        }
                        else
                        {
                //double alpha = dx * dy / (M_PI*h0*h0*t) ;

                                // normalizzatore *= alpha ;
                                // phieff[idxa] *= alpha / normalizzatore ;
                                //if (normalizzatore > 0) 
                                        phieff[idxa] /= normalizzatore ;
                                //else
                                //      phieff[idxa] *= alpha ;
                        }


//                                        if(xa*FT2M < 5005.0 && xa*FT2M  > 4995.0)fprintf(fphi,"%f %f %f %f %f\n",xa*FT2M,ya*FT2M,phi_box[idxa],phieff[idxa],dataLS->phi[idxa]);
//                                        if(phieff[idxa] <= 0.9)fprintf(fp3,"\n %f %f %f",xa*FT2M,ya*FT2M,phieff[idxa]);
//                                        if(phieff[idxa] <= 0.5)fprintf(fp2,"\n %f %f %f",xa*FT2M,ya*FT2M,phieff[idxa]);
//                                        if(phieff[idxa] <= 0.3)fprintf(fp1,"\n %f %f %f",xa*FT2M,ya*FT2M,phieff[idxa]);
        }
    }
//        fclose(fp1);
//        fclose(fp2);
//        fclose(fp3);
//        fclose(fphi);
    //printf("alpha: %f\n", dx*dy/(3.14*h0*h0*t));

    return 0 ;

} /*** end of: get_phieff_Gauss() *************************************/


/**********************************************************************
 * Evaluation of phieff for output from ForeFire (integral of K*phi)                           *
 **********************************************************************/
int get_phieff_Gauss_FF(double t)
{   printf("\n Including Turbulence at time = %f", t);
    #pragma omp parallel for
    for (int ja = j0_g; ja < j1_g + ic; ja++)
    {
		for (int ia = i0_g; ia < i1_g + ic; ia++)
        {
            int idxa = ia + ja * Nx_g ;
            
            double xa = x_0_g + ia * dx + d_center_x ;
            double ya = y_0_g + ja * dy + d_center_y ;
            
            phieff[idxa] = 0.0 ;
            
            double normalizzatore = 0.0 ;
			
			int i0b = max(i0_g, ia-Nneighbour*(int)((h0*t)/dx)) ;
			int i1b = min(i1_g, ia+Nneighbour*(int)((h0*t)/dx)) ;
			int j0b = max(j0_g, ja-Nneighbour*(int)((h0*t)/dy)) ;
			int j1b = min(j1_g, ja+Nneighbour*(int)((h0*t)/dy)) ;	
    
            for (int jb = j0b; jb < j1b + ic; jb++)  
            {                       
                for (int ib = i0b; ib < i1b + ic; ib++) 
                {
                    int idxb = ib + jb * Nx_g ;
                    
                    double xb = x_0_g + ib * dx + d_center_x ;
                    double yb = y_0_g + jb * dy + d_center_y ;
                    					
					double e = 0.0 ;
                    
					if (flag_backward_compatibility)
					{	printf ("with flag \n");			
						if ( ((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya)) < (9.0*h0*h0*t) )
							e = exp(-((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))/(h0*h0*t)) ;
					}
					else
					{
						e = exp(-((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))/(h0*h0*t)) ;
					}
					
//					phieff[idxa] += dataLS->phi[idxb] * e ;	
					phieff[idxa] += phi_ff[idxb] * e ;	
					normalizzatore += e ; 

                }
            }
            
            
            // NB: alpha = dx*dy/(M_PI*twosigma2) = dx*dy/(M_PI*h0*h0*t)
            if (flag_backward_compatibility)
            {
				phieff[idxa] *= (dx*dy/(3.14*h0*h0*t)) ;
			}
			else
			{
					phieff[idxa] /= normalizzatore ;
			}

        }
    }

    return 0 ;
    
} /*** end of: get_phieff_Gauss_FF() *************************************/


/**********************************************************************
 * If Flag_only_lognorm ==1, to use only the firespotting algorithm
 * with the Firefront output
 * ******************************************************************* */

int read_firefront(){
	int iq, i,  ii, count = 0 ;
	FILE *fp_x1, *fp_y1, *fp_fh,*fc, *f_phiff, *fphi ;
	char c, fn_firecontour[60], fn_phi_ff[60]; 
	sprintf(filename_xarray,"../data/FF_LS/xyarray.txt");
	sprintf(fn_phi_ff,"phiff.txt");

	fp_x1 = fopen(filename_xarray,"r");
	f_phiff = fopen(fn_phi_ff, "w");
//-------------------------------------------------------------------

    // Check if file exists
        if (fp_x1 == NULL)
    {
        printf("Could not open file %s", fp_x1);
        return 0;
    }

    // Extract characters from file and store in character c
        for (c = getc(fp_x1); c != EOF; c = getc(fp_x1))
        if (c == '\n') // Increment count if this character is newline
        count = count + 1;

    // Close the file
        //printf("The file %s has %d lines\n ", fp_x1, count);
        ntotal = count;
        fclose(fp_x1);

	fp_x1 = fopen(filename_xarray,"r");
//-----------------------------------------------------------------------	
	xf = (double*)calloc(ntotal, sizeof(double));
	yf = (double*)calloc(ntotal, sizeof(double));
	fhistory = (double*)calloc(ntotal, sizeof(double));

	for (int i = 0; i < ntotal ; i++)
	{
		fscanf(fp_x1,"%lg %lg %lg",&xf[i],&yf[i],&fhistory[i]);
	
	}
	fclose (fp_x1);

	int total_index = (int)fhistory[ntotal-1];
	//printf("total index = %d\n", total_index);
	for (int index = 1 ; index < (total_index)+1; index = index++)
	{
	nvert = 0;
	for (int i = 0; i < ntotal; i = i++)
		{ 
		if (fhistory[i] == index) nvert = nvert + 1;
		}
	//printf ("\n nvert %d ", nvert);

	xff = (double*)calloc(nvert, sizeof(double));
        yff = (double*)calloc(nvert, sizeof(double));
	ii = 0;
	for (int i = 0; i < ntotal; i = i++)
        {
		if(fhistory[i] == index) 
		{ 
		xff[ii] = xf[i];
        	yff[ii] = yf[i];
		ii=ii+1 ;       	
		}

        }
/****************************************************************
define phi for firefront                                        *
****************************************************************/
	
//	printf ("calculating through pnpoly");
	for (int i = i0_g; i < i1_g+ic ; i = i++)
	{
		for (int j = j0_g; j < j1_g+ic ; j = j++)
		{
			int idxa = i + j*Nx_g;
			double xa = (x_0_g + i * dx + d_center_x)*FT2M ;
	        	double ya = (y_0_g + j * dy + d_center_y)*FT2M ;
			
			if(phi_ff[idxa] < 1.0)phi_ff[idxa] = (double) pnpoly(nvert, xff , yff , xa, ya);
			
		}
	}
	free (xff);
	free (yff);
	}

	for (int i = i0_g; i < i1_g+ic ; i = i++)
        {
                for (int j = j0_g; j < j1_g+ic ; j = j++)
                {
			int idxa = i + j*Nx_g;
                        double xa = (x_0_g + i * dx + d_center_x)*FT2M ;
                        double ya = (y_0_g + j * dy + d_center_y)*FT2M ;
			if (fuelMap[idxa] < 1)phi_ff[idxa] = 0.0;
			if (phi_ff[idxa] == 1.)fprintf (f_phiff," %f %f %f\n",xa,ya,phi_ff[idxa]);
		}
	}
	
	
	fclose(f_phiff);
	return 0;
}

/*************************************************************** 
*find the points inside the polygon			        *
****************************************************************/
int pnpoly(int nvert, double *xff, double *yff, double testx, double testy)
{
//  int i, j, c = 1; //c = 1 , to ensure 0 inside and 1 outside;
  int i, j, c = 0; //c = 0 , to ensure 1 inside and 0 outside;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((yff[i]>testy) != (yff[j]>testy)) &&
	 (testx < (xff[j]-xff[i]) * (testy-yff[i]) / (yff[j]-yff[i]) + xff[i]) )
       c = !c;
  }
/*  for (i = 0; i < nvert; i= i++)
	{
	if(fabs(xff[i]-testx) < 2.0 && fabs(yff[i]-testy) < 2.0)
	//printf("\n %f %f", xff[i],yff[i]);
	 c = 1;
	}*/
  return c;
}

/**********************************************************************
 * Evaluation of phieff (integral of K*phi)                           *
 **********************************************************************/
int get_phieff_Lognorm(double t)

{   FILE *fphi;
        fphi = fopen("phieff_slice.txt", "w+");
    #pragma omp parallel for 
    for (int ja = j0_g; ja < j1_g + ic; ja++)
    {
                for (int ia = i0_g; ia < i1_g + ic; ia++)
        {
            int idxa = ia + ja * Nx_g ;

            double xa = x_0_g + ia * dx + d_center_x ;
            double ya = y_0_g + ja * dy + d_center_y ;
             //if (dataLS->phi[idxa] <= 0.5 ) phi_box[idxa]= 0.0;
             //if (dataLS->phi[idxa] > 0.5) phi_box[idxa]= 1.0;

            double wspd, wdir ;

            if (flag_wind_const)
            {
                wspd = WindSpd ;
                wdir = WindDir ;

            }
            else
            {
                wspd = fun_WindSpd(t, xa, ya, WindSpd) ;
                wdir = fun_WindDir(t, xa, ya) ;
            }

            double sint = deg2sin(wdir) ;
            double cost = deg2cos(wdir) ;
            double Ut = wspd * FTMIN2MS ;

            //get_lognorm_par(If, Ut) ; // THIS IS NOT THREAD-SAFE !!!

            double lognorm_Fr, lognorm_mu, lognorm_s, lognorm_mean_l ;
            double Iff = If / 1e3 ; // [kW/m] -> [MW/m]

            // Calculate Fr
            if (flag_kernel == 60)
            {
                double Lc = pow( If/(rho_amb*c_pg*T_amb*pow(g_acc, 0.5)), 2.0/3.0) ; // here we need If, not Iff
                lognorm_Fr = Ut / pow(g_acc*Lc, 0.5) ;
            }
            else if (flag_kernel == 61)
            {
                lognorm_Fr = 0.0 ;
            }
            else if (flag_kernel == 62)
 {
                lognorm_Fr = INFINITY ;
            }

                //printf("%f\n", lognorm_Fr);
            // Calculate mu, s
            if (lognorm_Fr < 1)
            {
                lognorm_mu = 1.47 * pow(Iff, 0.54) * pow(Ut, -0.55) + 1.14 ;
                lognorm_s  = 0.86 * pow(Iff, -0.21) * pow(Ut, 0.44) + 0.19 ;
            }
            else // lognorm_Fr >= 1
            {
                lognorm_mu = 1.32 * pow(Iff, 0.26) * pow(Ut, 0.11) - 0.02 ;
                lognorm_s  = 4.95 * pow(Iff, -0.01) * pow(Ut, -0.02) - 3.48 ;
                //lognorm_mu = 6.59;
               //lognorm_s = 0.125;
            }

            // see note in get_lognorm_par()

            lognorm_mu += LN_M2FT ;

            lognorm_mean_l = exp(lognorm_mu + 0.5*lognorm_s*lognorm_s) ;
//            lognorm_mu = 15.0;
//            lognorm_s = 5.5;
            //printf("%f %f \n", lognorm_mu,lognorm_s);      

            phieff[idxa] = 0.0 ;

            double normalizzatore = 0.0 ;

                        int i0b = max(i0_g, ia-Nneighbour) ;//*(int)(h0*t/dx)) ;
                        //int i0b = max(i0_g, ia-Nneighbour *(int)(h0*t/dx)) ;
                        int i1b = min(i1_g, ia+Nneighbour) ;//*(int)(h0*t/dx)) ;
                        //int i1b = min(i1_g, ia+Nneighbour *(int)(h0*t/dx)) ;
                        int j0b = max(j0_g, ja-Nneighbour) ;//*(int)(h0*t/dy)) ;
                        //int j0b = max(j0_g, ja-Nneighbour *(int)(h0*t/dy)) ;
                        int j1b = min(j1_g, ja+Nneighbour) ;//*(int)(h0*t/dy)) ;
                        //int j1b = min(j1_g, ja+Nneighbour *(int)(h0*t/dy)) ;



                        for (int jb = j0b; jb < j1b + ic; jb++)
                        {
                                for (int ib = i0b; ib < i1b + ic; ib++)
                                {
                                        int idxb = ib + jb * Nx_g ;

                                        //if (dataLS->phi[idxb] > PHI_THRESHOLD) continue ;

                                        double xb = x_0_g + ib * dx + d_center_x ;
                                        double yb = y_0_g + jb * dy + d_center_y ;

                                        double r = sqrt( (ya-yb)*(ya-yb) + (xa-xb)*(xa-xb) ) ;
                                        double cosu = (ya-yb) / r ;
                                        double sinu = (xa-xb) / r ;

                                        double nU_dot_n = nx[idxb]*sint + ny[idxb]*cost ;
                    double e ;
                    double twosigma2 = h0*h0*t ;

                    // punti sottovento e con n*nU>0
                    if ( (r>0.0) && fabs(cosu-cost) < 0.001 && fabs(sinu-sint) < 0.001 && nU_dot_n > 0.0 )
                    //if ( (r>0.0) && fabs(cosu-cost) < 0.001 && fabs(sinu-sint) < 0.001 && nU_dot_n > 0.0 && ya*FT2M > 2500.0)
                    //if ( (r>0.0) && fabs(cosu-cost) < 0.001 && fabs(sinu-sint) < 0.001  && ya*FT2M > 2500.0)
                    //if ( (r>0.0) && fabs(cosu-cost) < 0.001 && fabs(sinu-sint) < 0.001  && ya*FT2M > 1500.0)
                                        {
                                                double sum_quadratura = 0.0 ;
                                                for (int i = 0; i < Nquad; i++)
                                                {
//                                                      e = w_q[i] * exp( -( r*r + exp(2.0*sqrt(2.0)*lognorm_s*x_q[i]) - 
//                                                              2.0*r*exp(sqrt(2.0)*lognorm_s*x_q[i]) ) / twosigma2 + 
//                                                              sqrt(2.0)*lognorm_mu*x_q[i]/lognorm_s - 0.5*pow(lognorm_mu/lognorm_s,2) ) ;
                                                         e = w_q[i] * exp( -( r*r + exp(2.0*sqrt(2.0)*lognorm_s*x_q[i]) -
                                                                2.0*r*exp(sqrt(2.0)*lognorm_s*x_q[i]) ) / twosigma2 +
                                                                sqrt(2.0)*(fabs(lognorm_mu))*x_q[i]/lognorm_s - 0.5*pow(fabs(lognorm_mu)/lognorm_s,2) ) ;

                                        //printf("%f %f %f\n",xa-xb,r,ya-yb);
                                                        if ( e<INFINITY ) sum_quadratura += e ;
                                                }
                                                e = sum_quadratura / (pow(M_PI, 1.5)*twosigma2) ;

                        }
                        else
                                        {
                                                e = exp(-r*r/twosigma2) / (M_PI*twosigma2) ;
//                                              printf("we are here\n");
                        }
                        phieff[idxa] += dataLS->phi[idxb] * e ;
//                      phieff[idxa] += phi_box[idxb] * e ;
                        normalizzatore += e ;
                                }
            }

            //phieff[idxa] *= dx * dy ;      
            phieff[idxa] /= normalizzatore  ;
            //if(xa*FT2M < 5005.0 && xa*FT2M  > 4995.0)fprintf(fphi,"%f %f %f %f %f\n",xa*FT2M,ya*FT2M,phi_box[idxa],phieff[idxa],dataLS->phi[idxa]);

        }
    }
    fclose(fphi);
    return 0 ;

} /*** end of: get_phieff_Lognorm() *****************************************/



/**********************************************************************
 * Evaluation of phieff (integral of K*phi)  for output from ForeFire                         *
 **********************************************************************/
int get_phieff_Lognorm_FF(double t)
{	printf ("Including Fire-spotting");
//	FILE *fphi;
//	fphi = fopen("phieff_slice.txt", "w+");
	
    #pragma omp parallel for 
    for (int ia = i0_g; ia < i1_g + ic; ia++)
    {
                for (int ja = j0_g; ja < j1_g + ic; ja++)
        {
            int idxa = ia + ja * Nx_g ;
            double xa = x_0_g + ia * dx + d_center_x ;
            double ya = y_0_g + ja * dy + d_center_y ;

            double wspd, wdir ;

            if (flag_wind_const)
            {
                wspd = WindSpd ;
                wdir = WindDir ;

            }
            else
            {
                wspd = fun_WindSpd(t, xa, ya, WindSpd) ;
                wdir = fun_WindDir(t, xa, ya) ;
            }

            double sint = deg2sin(wdir) ;
            double cost = deg2cos(wdir) ;
            double Ut = wspd * FTMIN2MS ;

            //get_lognorm_par(If, Ut) ; // THIS IS NOT THREAD-SAFE !!!

            double lognorm_Fr, lognorm_mu, lognorm_s, lognorm_mean_l ;
            double Iff = If / 1e3 ; // [kW/m] -> [MW/m]

            // Calculate Fr
            if (flag_kernel == 60)
            {
                double Lc = pow( If/(rho_amb*c_pg*T_amb*pow(g_acc, 0.5)), 2.0/3.0) ; // here we need If, not Iff
                lognorm_Fr = Ut / pow(g_acc*Lc, 0.5) ;
            }
            else if (flag_kernel == 61)
            {
                lognorm_Fr = 0.0 ;
 }
            else if (flag_kernel == 62)
            {
                lognorm_Fr = INFINITY ;
            }

            // Calculate mu, s
            if (lognorm_Fr < 1)
            {
                lognorm_mu = 1.47 * pow(Iff, 0.54) * pow(Ut, -0.55) + 1.14 ;
                lognorm_s  = 0.86 * pow(Iff, -0.21) * pow(Ut, 0.44) + 0.19 ;
            }
            else // lognorm_Fr >= 1
            {
                lognorm_mu = 1.32 * pow(Iff, 0.26) * pow(Ut, 0.11) - 0.02 ;
                lognorm_s  = 4.95 * pow(Iff, -0.01) * pow(Ut, -0.02) - 3.48 ;
            }

            // see note in get_lognorm_par()

            lognorm_mu += LN_M2FT ;

            lognorm_mean_l = exp(lognorm_mu + 0.5*lognorm_s*lognorm_s) ;



            phieff[idxa] = 0.0 ;

            double normalizzatore = 0.0 ;

                        int i0b = max(i0_g, ia-Nneighbour) ;//*(int)(h0*t/dx)) ;
                        int i1b = min(i1_g, ia+Nneighbour) ;//*(int)(h0*t/dx)) ;
                        int j0b = max(j0_g, ja-Nneighbour) ;//*(int)(h0*t/dy)) ;
                        int j1b = min(j1_g, ja+Nneighbour) ;//*(int)(h0*t/dy)) ;

                        for (int jb = j0b; jb < j1b + ic; jb++)
                        {
                                for (int ib = i0b; ib < i1b + ic; ib++)
                                {
                                        int idxb = ib + jb * Nx_g ;

                                        //if (dataLS->phi[idxb] > PHI_THRESHOLD) continue ;

                                        double xb = x_0_g + ib * dx + d_center_x ;
                                        double yb = y_0_g + jb * dy + d_center_y ;

                                        double r = sqrt( (ya-yb)*(ya-yb) + (xa-xb)*(xa-xb) ) ;
                                        double cosu = (ya-yb) / r ;
                                        double sinu = (xa-xb) / r ;
                                        //double nU_dot_n = nx[idxb]*sint + ny[idxb]*cost ;

                    double e ;
                    double twosigma2 = h0*h0*t ;

                    // punti sottovento e con n*nU>0
                    if ( (r>0.0) && fabs(cosu-cost) < 0.001 && fabs(sinu-sint) < 0.001 && ya*FT2M > 2500.0)
                    //if ( (r>0.0) && fabs(cosu-cost) < 0.001 && fabs(sinu-sint) < 0.001 && nU_dot_n > 0.0 )
                                        {
                                                double sum_quadratura = 0.0 ;
                                                for (int i = 0; i < Nquad; i++)
                                                {
                                                        e = w_q[i] * exp( -( r*r + exp(2.0*sqrt(2.0)*lognorm_s*x_q[i]) -
                                                                2.0*r*exp(sqrt(2.0)*lognorm_s*x_q[i]) ) / twosigma2 +
                                                                sqrt(2.0)*lognorm_mu*x_q[i]/lognorm_s - 0.5*pow(lognorm_mu/lognorm_s,2) ) ;
                                                        if ( e<INFINITY ) sum_quadratura += e ;
                                                }
                                                e = sum_quadratura / (pow(M_PI, 1.5)*twosigma2) ;
                                                //sum_quad = sum_quadratura ;
                                                //printf ("\n %f %f", sum_quadratura, e);       
                        }
                        else
                                        {
                                                e = exp(-r*r/twosigma2) / (M_PI*twosigma2) ;
                        }
                        phieff[idxa] += phi_ff[idxb] * e ;
                        normalizzatore += e ;
                                }
            }

            phieff[idxa] /= normalizzatore  ;
//            if(ya*FT2M < 2505.0 && ya*FT2M  > 2495.0)fprintf(fphi,"%f %f %f %f %f\n",xa*FT2M,ya*FT2M,phi_ff[idxa],phieff[idxa],dataLS->phi[idxa]);
    }
    }
//	fclose(fphi);
    return 0 ;

} /*** end of: get_phieff_Lognorm_FF() *****************************************/


int get_phieff_Lognorm1(double t)
{    
//double dx = 30.0 ; double dy =30.0 ;
double sum_quadratura;
//    #pragma omp parallel for
    for (int ja = j0_g; ja < j1_g + ic; ja++)
    {
		for (int ia = i0_g; ia < i1_g + ic; ia++)
        {
            int idxa = ia + ja * Nx_g ;
            double xa = x_0_g + ia * dx + d_center_x ;
            double ya = y_0_g + ja * dy + d_center_y ;
            
            double wspd, wdir ;
            
            if (flag_wind_const)
            {
                wspd = WindSpd ;
                wdir = WindDir ;
                
            }
            else 
            {
                wspd = fun_WindSpd(t, xa, ya, WindSpd) ;
                wdir = fun_WindDir(t, xa, ya) ;
            }
            
            double sint = deg2sin(wdir) ;
            double cost = deg2cos(wdir) ;
            double Ut = wspd * FTMIN2MS ;
            
            //get_lognorm_par(If, Ut) ; // THIS IS NOT THREAD-SAFE !!!
            
            double lognorm_Fr, lognorm_mu, lognorm_s, lognorm_mean_l ;
            double Iff = If / 1e3 ; // [kW/m] -> [MW/m]

            // Calculate Fr
            if (flag_kernel == 30)
            {
                double Lc = pow( If/(rho_amb*c_pg*T_amb*pow(g_acc, 0.5)), 2.0/3.0) ; // here we need If, not Iff
                lognorm_Fr = Ut / pow(g_acc*Lc, 0.5) ;
            }
            else if (flag_kernel == 31)
            {
                lognorm_Fr = 0.0 ;
            }
            else if (flag_kernel == 32)
            {
                lognorm_Fr = INFINITY ;
            }
            
            // Calculate mu, s
            if (lognorm_Fr < 1)
            {
                lognorm_mu = 1.47 * pow(Iff, 0.54) * pow(Ut, -0.55) + 1.14 ;
                lognorm_s  = 0.86 * pow(Iff, -0.21) * pow(Ut, 0.44) + 0.19 ;
            }
            else // lognorm_Fr >= 1
            {
                lognorm_mu = 1.32 * pow(Iff, 0.26) * pow(Ut, 0.11) - 0.02 ;
                lognorm_s  = 4.95 * pow(Iff, -0.01) * pow(Ut, -0.02) - 3.48 ;
            }
            
            // see note in get_lognorm_par()
            
            lognorm_mu += LN_M2FT ;
            
            lognorm_mean_l = exp(lognorm_mu + 0.5*lognorm_s*lognorm_s) ;
            
                  
//          printf("lognorm_mu, lognorm_s");
//          printf("%7.2f %7.2f %9d\n",lognorm_mu, lognorm_s, idxa);   
            phieff[idxa] = 0.0 ;
		
            double normalizzatore = 0.0 ;
			
	/*		int i0b = max(i0_g, ia-Nneighbour) ;//*(int)(h0*t/dx)) ;
			int i1b = min(i1_g, ia+Nneighbour) ;//*(int)(h0*t/dx)) ;
			int j0b = max(j0_g, ja-Nneighbour) ;//*(int)(h0*t/dy)) ;
			int j1b = min(j1_g, ja+Nneighbour) ;//*(int)(h0*t/dy)) ;
			
			for (int jb = j0b; jb < j1b + ic; jb++)  
			{                       
				for (int ib = i0b; ib < i1b + ic; ib++) 
				{
					int idxb = ib + jb * Nx_g ;
                   
					//if (dataLS->phi[idxb] > PHI_THRESHOLD) continue ;
		*/
	    	
	              for (ib = 0; ib < ntotal ; ib++)
			{
/*					double xb = x_0_g + ib * dx + d_center_x ;
					double yb = y_0_g + jb * dy + d_center_y ;
*/			 
            		                
					double xb = xf[ib] ;
					double yb = yf[ib] ;
										

					double r = sqrt( (ya-yb)*(ya-yb) + (xa-xb)*(xa-xb) ) * (100.0/30.0) ; //M2FT 16 Jan 2015
					double r1 = sqrt( (ya-yb)*(ya-yb) + (xa-xb)*(xa-xb) )  ; 
					double cosu = (ya-yb) / r1 ;
					double sinu = (xa-xb) / r1 ;						
					//printf ("%f %f %f %f %f\n", r,cosu,sinu,cost,sint);	
					//double nU_dot_n = nx[idxb]*sint + ny[idxb]*cost ;
                    
                    double e ;
                    double twosigma2 = h0*h0*t ;
                    //printf ("%9.2f\n",twosigma2);
                    // punti sottovento e con n*nU>0
                    //if ( (r>0.0) && fabs(cosu-cost) < 0.001 && fabs(sinu-sint) < 0.001 && nU_dot_n > 0.0 )
                    if ( (r>0.0) && fabs(cosu-cost) < 0.001 && fabs(sinu-sint) < 0.001 && fhistory[ib] > 0.0 )
						
						{											
                                                double sum_quadratura = 0.0 ;
						for (int i = 0; i < Nquad; i++)
						{	//printf ("ok");
							e = w_q[i] * exp( -( r*r + exp(2.0*sqrt(2.0)*lognorm_s*x_q[i]) - 
								2.0*r*exp(sqrt(2.0)*lognorm_s*x_q[i]) ) / twosigma2 + 
								sqrt(2.0)*lognorm_mu*x_q[i]/lognorm_s - 0.5*pow(lognorm_mu/lognorm_s,2) ) ;
							if ( e<INFINITY ) sum_quadratura += e ;
						
						}

		    //include contribution to the integral from the burnt points:			
						for (j = 0; j < ntotal; j++)
						{
						if (fhistory[j] < 1.0)
						{    	        
							 
                                        		double xb = xf[j] ;
                                        		double yb = yf[j] ;
                                        		double r = sqrt( (ya-yb)*(ya-yb) + (xa-xb)*(xa-xb) ) * (100.0/30.0) ; //M2FT 16 Jan 2015

							//double sum_quadratura  ;
							for (int i = 0; i < Nquad; i++)
						{       //printf ("ok");
                                                        e = w_q[i] * exp( -( r*r + exp(2.0*sqrt(2.0)*lognorm_s*x_q[i]) -
                                                                2.0*r*exp(sqrt(2.0)*lognorm_s*x_q[i]) ) / twosigma2 +
                                                                sqrt(2.0)*lognorm_mu*x_q[i]/lognorm_s - 0.5*pow(lognorm_mu/lognorm_s,2) ) ;
                                                        if ( e<INFINITY ) sum_quadratura += e ;
							//printf ("\n new %d %d %f",j,i,e);
                                                }
                                                //printf ("%f %f %f \n",M_PI,twosigma2,(pow(M_PI, 1.5)*twosigma2)) ;
                        
						sum_quad = sum_quadratura ;	
						e = sum_quadratura / (pow(M_PI, 1.5)*twosigma2) ;
						//printf ("%f %f %f \n",M_PI,twosigma2,(pow(M_PI, 1.5)*twosigma2)) ;	
	                
//	                else
//					{
//						e = exp(-r*r/twosigma2) / (M_PI*twosigma2) ;
	                }
			
			}
			
			//printf ("\n sum_quad %f %f %f ", sum_quadratura,sum_quad,  e) ;
			
			
	                //phieff[idxa] += dataLS->phi[idxb] * e ;
			//printf ("\n %f ip ", sum_quad);
	                //phieff[idxa] += sum_quad ;
			//printf ("\n %f sum", sum_quadratura);
	                }
			else
                                      {
                                              e = exp(-r*r/twosigma2) / (M_PI*twosigma2) ;
                        }
			
	                phieff[idxa] += e ;
			}
	                //sum_quad[idxa] += sum_quadratura ;

                        
	                //phieff[idxa] += 1.0;
	                //normalizzatore += e ;
			
            			
            
            //phieff[idxa] *= dx * dy ;      
            printf("\n %f",phieff[idxa]);
	    //phieff[idxa] /= normalizzatore  ;
	       
        
    }
    }

    return 0 ;
    
} /*** end of: get_phieff_Lognorm1() *****************************************/



/**********************************************************************
 * Evaluation of phieff (integral of K*phi)                           *
 **********************************************************************/
int get_phieff_Mainardi(double t)
{
   
    // h0 = 2*sqrt(D) ; sigma^2=2*D*t   =>   twosigma2 = 4*D*t = h0*h0*t
    // NB: alpha = dx*dy/(M_PI*twosigma2) = dx*dy/(M_PI*h0*h0*t) = dx*dy/(4*M_PI*D*t)
    
    //alpha = dx * dy / (M_PI*h0*h0*t) ;
    
    double f = pow(Mbeta_D, -0.5) * pow(t, -Mbeta_beta) ;
    
    #pragma omp parallel for             
    for (int ja = j0_g; ja < j1_g + ic; ja++)                              
    {
		for (int ia = i0_g; ia < i1_g + ic; ia++)
        {
            int idxa = ia + ja * Nx_g ;
            
            double xa = x_0_g + ia * dx + d_center_x ;
            double ya = y_0_g + ja * dy + d_center_y ;
            
            phieff[idxa] = 0.0 ;
            
            double normalizzatore = 0.0 ;
			
			int i0b = max(i0_g, ia-Nneighbour*(int)((h0*t)/dx)) ;
			int i1b = min(i1_g, ia+Nneighbour*(int)((h0*t)/dx)) ;
			int j0b = max(j0_g, ja-Nneighbour*(int)((h0*t)/dy)) ;
			int j1b = min(j1_g, ja+Nneighbour*(int)((h0*t)/dy)) ;	
            
            for (int jb = j0b; jb < j1b + ic; jb++)  
            {                       
                for (int ib = i0b; ib < i1b + ic; ib++) 
                {
                    int idxb = ib + jb * Nx_g ;
                    
                    double xb = x_0_g + (ib-i0_g+0.5) * dx ;
                    double yb = y_0_g + (jb-j0_g+0.5) * dy ;
                    					
					double r = sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya)) ;
					double e = M_interp(f*r, 0) ;
										
					phieff[idxa] += dataLS->phi[idxb] * e ;	
					normalizzatore += e ; 

                }
            }
            
			phieff[idxa] /= normalizzatore ;
			            
        }
    }

    return 0 ;
    
} /*** end of: get_phieff_Mainardi() **********************************/




/**********************************************************************
 * Evaluation of phieff (integral of K*phi)                           *
 **********************************************************************/
int get_phieff(int Niter, double t_old, double t)
{
	//printf ("flag_kernel_lognorm:\n",flag_kernel_Lognorm);
	//get_phieff_lsm() ; // per tutti!
				
	if ( flag_kernel_LSM || Niter == 1 )
		get_phieff_LSM();
//	if ( flag_kernel_LSM  )
		
	else if (flag_kernel_Gauss) 
		get_phieff_Gauss(t) ;
		
	else if (flag_kernel_Lognorm) 
		get_phieff_Lognorm(t) ;
	
        else if (flag_kernel_Gauss_FF)
		{
		if ( t == 0.0 ) {
			get_phieff_LSM_FF() ;
				}
		else {
                	get_phieff_Gauss_FF(t) ;
			}
		}

        else if (flag_kernel_Lognorm_FF)
		{
		if ( t == 0.0 ) {
                        get_phieff_LSM_FF() ;
                                }
                else {
                get_phieff_Lognorm_FF(t) ;
		}
		}
//	else if (flag_kernel_Mainardi)
//		get_phieff_Mainardi(t) ;
	
	
	return 0 ;
}

/**********************************************************************
 * Evaluation of psi and fuelMap                                      *
 **********************************************************************/
int get_psi(LSMLIB_REAL dt, LSMLIB_REAL t)
{   
    for (int i = i0_g; i < i1_g + ic; i++)
    {
//              for (int i = i0_g; i < i1_g + ic; i++)
                for (int j = j0_g; j < j1_g + ic; j++)
        {
            int idx = i + j * Nx_g ;

            double wspd, wdir ;
             double x = x_0_g + i * dx + d_center_x ;
             double y = y_0_g + j * dy + d_center_y ;

            if (flag_wind_const)
            {
                wspd = WindSpd ;
                wdir = WindDir ;
            }
            else
            {
                double x = x_0_g + i * dx + d_center_x ;
                double y = y_0_g + j * dy + d_center_y ;

                wspd = fun_WindSpd(t, x, y, WindSpd) ;
                wdir = fun_WindDir(t, x, y) ;
            }

            double sint = deg2sin(wdir) ;
            double cost = deg2cos(wdir) ;

            double nU_dot_n = (nx[idx]*sint + ny[idx]*cost) ;

            double tc ;

            if ( nU_dot_n > 0.0 )
                tc = tc_dw ;
            else
                tc = tc_up ;
            psi[idx] += (1.0-phieff[idx]) * dt / tc ;
//          psi[idx] += (1.0-phieff[idx])  ;
            if (phieff[idx] <= PHIEFF_THRESHOLD)
                psi[idx] = PSI_THRESHOLD ;

            // if psi>=1 and there is fuel  =>  it's burnt, so phi=0
            if (psi[idx] >= PSI_THRESHOLD)
            {
                psi[idx] = PSI_THRESHOLD ;

                int idxMap_fuel ;

                if (flag_fireMap_const && !flag_fire_obstacle)
                    idxMap_fuel = 0 ;
                else
                    idxMap_fuel = idx ;

                if ( (!flag_kernel_LSM) && flag_hotfront && fuelMap[idxMap_fuel] > 0  )
                   {
                         dataLS->phi[idx] = 0.0 ;
                   }
            }
        }
    }
        return 0 ;

} /*** end of: get_psi() **********************************************/


/**********************************************************************
 * Evaluation of psi and fuelMap for output from ForeFire                                   *
 **********************************************************************/
int get_psi_FF(LSMLIB_REAL dt, LSMLIB_REAL t)
{	FILE *f_new_ign, *f_psi, *f_psi1;
	size_t size;
	
	FILE *fp_psi;
	fp_psi =  fopen("../data/FF_LS/psi.txt","a+");
//        printf("reading psi data\n");
	fseek (fp_psi, 0, SEEK_END);
	size = ftell(fp_psi);
//	printf(" size = %d\n", size);
	fseek (fp_psi, 0, SEEK_SET);
	//rewind(fp_psi);
	if( size > 0) 
	{
	for (int i = i0_g; i < i1_g + ic; i++)
    		{
                for (int j = j0_g; j < j1_g + ic; j++)
		{ 	int idx = i + j * Nx_g; 
			fscanf(fp_psi, "%lg", & psi[idx]);
		}
	}
	}
	else
	{for(int i = 0 ; i < N_g; i++)
        { psi[i] = 0.0;}
        }
	fclose(fp_psi);

	f_new_ign = fopen("../data/FF_LS/new_ignition.txt", "w+");
	f_psi = fopen("../data/FF_LS/psi.txt", "w+");
	f_psi1 = fopen("../data/FF_LS/psi1.txt", "w+");
//    #pragma omp parallel for
    for (int i = i0_g; i < i1_g + ic; i++)
    {
		for (int j = j0_g; j < j1_g + ic; j++)
        {
            int idx = i + j * Nx_g ;
                double x = x_0_g + i * dx + d_center_x ;
                double y = y_0_g + j * dy + d_center_y ;
           
		 double wspd, wdir ;
            
            if (flag_wind_const)
            {
                wspd = WindSpd ;
                wdir = WindDir ;
            }
            else 
            {
                double x = x_0_g + i * dx + d_center_x ;
                double y = y_0_g + j * dy + d_center_y ;
                
                wspd = fun_WindSpd(t, x, y, WindSpd) ;
                wdir = fun_WindDir(t, x, y) ;
            }
            
            double sint = deg2sin(wdir) ;
            double cost = deg2cos(wdir) ;
            
            
            double tc ;
            
//            if ( nU_dot_n > 0.0 )
              tc = tc_dw ;
//            else
	if (y < 2500*M2FT) tc = tc_up;
	//if (x < 2500*M2FT) tc = tc_up;
	else tc = tc_dw ;
                
            psi[idx] += (phieff[idx]) * dt / tc ;
            if (phieff[idx] >= PHIEFF_THRESHOLD)
		{// printf("\n %f %f %f", phieff[idx], psi[idx], phi_ff[idx]);
              		psi[idx] = PSI_THRESHOLD;
		}
		
            if (psi[idx] >= PSI_THRESHOLD)
                {
		psi[idx] = PSI_THRESHOLD ;
                fprintf(f_psi1," %f %f %f\n",x*FT2M,y*FT2M,psi[idx]);
		int idxMap_fuel ;
                
		if (flag_fireMap_const && !flag_fire_obstacle)
                    idxMap_fuel = 0 ;
                else
                    idxMap_fuel = idx ;
                    
                if ( fuelMap[idxMap_fuel] > 0 && phi_ff[idx] < 1.0 )
		{
                	phi_ff[idx] = 1.0 ;
			fprintf(f_new_ign," %f %f %f %d %f\n",x*FT2M,y*FT2M,psi[idx], fuelMap[idx], phi_ff[idx]); 
		}
		}
                fprintf(f_psi, " %f\n",psi[idx]);
        }
    }
        

	fclose(f_psi);
	fclose(f_psi1);
	fclose(f_new_ign);
	return 0 ;
		
} /*** end of: get_psi_FF() **********************************************/	



/**********************************************************************
 * Evaluation of time of burning                                      *
 **********************************************************************/
int get_tburn(double t)
{	
	if (flag_save_tburn_phi)
	{
		for (int idx = 0; idx < N_g; idx++)
			if ( dataLS->phi[idx] >= PHI_THRESHOLD ) 
				tburn_phi[idx] = t ; // ex "time"
	}
	
	if (flag_save_tburn_phieff)
	{
		for (int idx = 0; idx < N_g; idx++)
			if ( phieff[idx] >= PHIEFF_THRESHOLD ) 
				tburn_phieff[idx] = t ; // ex "gptime"
	}
	
	if (flag_save_tburn_psi)
	{
		for (int idx = 0; idx < N_g; idx++)
			if ( psi[idx] >= PSI_THRESHOLD && tburn_psi[idx] == T_NO_BURN ) 
				tburn_psi[idx] = t ;
	}
	
	return 0 ;
		
} /*** end of: get_tburn() ********************************************/



/**********************************************************************
 * Check if the fire has reached the border                           *
 **********************************************************************/
unsigned int check_fire_border()
{
	unsigned int flag_fire_at_border = 0 ;
	
	for (int i = i0_i; i < i1_i + ic ; i++)
	{
		// south and north
		//if (flag_backward_compatibility)
		//{
		//	if (dataLS->phi[i+j0_f*Nx_g] < PHI_THRESHOLD) flag_fire_at_border = 1 ;
		//	if (dataLS->phi[i+j1_f*Nx_g] < PHI_THRESHOLD) flag_fire_at_border = 1 ;
		//}
		//else
		//{
			if (psi[i+j0_f*Nx_g] >= PSI_THRESHOLD) flag_fire_at_border = 1 ;
			if (psi[i+j1_f*Nx_g] >= PSI_THRESHOLD) flag_fire_at_border = 1 ;
		//}
	}
	
	for (int j = j0_i; j < j1_i + ic ; j++)
	{
		// west and east
		//if (flag_backward_compatibility)
		//{
		//	if (dataLS->phi[i0_f+j*Nx_g] < PHI_THRESHOLD) flag_fire_at_border = 1 ;
		//	if (dataLS->phi[i1_f+j*Nx_g] < PHI_THRESHOLD) flag_fire_at_border = 1 ;
		//}
		//else
		//{
			if (psi[i0_f+j*Nx_g] >= PSI_THRESHOLD) flag_fire_at_border = 1 ;
			if (psi[i1_f+j*Nx_g] >= PSI_THRESHOLD) flag_fire_at_border = 1 ;
		//}
	}
		
	return flag_fire_at_border ;
		
} /*** end of: check_fire_border() ************************************/
		
		
		
/**********************************************************************
 * Evaluation of the time step according to the CFL condition         *
 **********************************************************************/
double get_dt(double cfl) 
{  
	double dtmin = INFINITY ;
	
    for (int j = j0_i; j < j1_i + ic; j++)
    {
        for (int i = i0_i; i < i1_i + ic; i++) 
        {
            int idx = i + j * Nx_g ;
            
			double dt = min(dx, dy) / ros[idx] ;
			
			if (dt < dtmin) dtmin = dt ;
        }
    }
    
    return (cfl*dtmin) ;
    
} /*** end of: get_dt() ***********************************************/



/**********************************************************************
 * Free all the allocated memory                                      *
 **********************************************************************/
int free_all()
{
    destroyGrid(gridLS) ;
    destroyLSMDataArrays(dataLS) ;
    
    free(fuelMap) ;
    free(ignMap) ;
    //free(flameMap) ;
    free(slpMap) ;
    free(aspMap) ;
    free(wspdMap) ;
    free(wdirMap) ;
    free(m1Map) ;
    free(m10Map) ;
    free(m100Map) ;
    free(mherbMap) ;
    free(mwoodMap) ;
    
    free(ros) ;
    
    return 0 ;
    
} /*** end of: free_all() *********************************************/




/**********************************************************************
 * Get the filename of the data to be saved                           *
 **********************************************************************/
 char* get_filename_full(char* filename, int Niter, char* filename_extension)
{
	//char filename_full[MAX_CHAR_NFILE] ;
	
	if (Niter >= 0)
		sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename,
			Niter, filename_extension) ;
	else
		sprintf(filename_full, "%s_%s_end.%s", filename_ws, filename,
			filename_extension) ;
			
	return filename_full ;
	
} /*** end of: get_filename_full() ************************************/			
				



/**********************************************************************
 * Save grid information to txt file                                  *
 **********************************************************************/
int save_datapre_to_file(char* filename)   
{
	FILE *fp ;
	char filename_full[MAX_CHAR_NFILE] ;
    
	sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename_data, 0, filename_extension_txt) ;
    
    fp = fopen(filename_full, "w") ;
    
    fprintf(fp, "Nthreads = %d; Nthreads_tot = %d; Nthreads_max = %d;\n", 
                Nthreads, Nthreads_tot, Nthreads_max) ;
    fprintf(fp, "flag_sanitycheck = %d; flag_backward_compatibility = %d; ic=%d;\n", 
                flag_sanitycheck, flag_backward_compatibility, ic) ;
    fprintf(fp, "\n") ;
    
    fprintf(fp, "x_0 = %f; x_1 = %f;\n", x_0, x_1) ;
    fprintf(fp, "y_0 = %f; y_1 = %f;\n", y_0, y_1) ;
    fprintf(fp, "dx = %f; dy = %f;\n\n", dx, dy) ;
    fprintf(fp, "Lx = %f; Ly = %f; A = %f;\n", Lx, Ly, A) ;
    fprintf(fp, "\n") ;
    
    fprintf(fp, "Nx = %d; Ny = %d; N = %d;\n", Nx, Ny, N) ;
    fprintf(fp, "\n") ;
    
    fprintf(fp, "t_max = %f; deltat = %f;\n", t_max, deltat) ;
    fprintf(fp, "flag_deltat_const = %d; flag_t_end = %d;\n", 
                flag_deltat_const, flag_t_end) ;
    fprintf(fp, "\n") ;
    
    fprintf(fp, "x_0_g = %f; x_1_g = %f;\n", x_0_g, x_1_g) ;
    fprintf(fp, "y_0_g = %f; y_1_g = %f;\n", y_0_g, y_1_g) ;
    fprintf(fp, "flag_d_center = %d; d_center_x = %f; d_center_y = %f;\n", 
                 flag_d_center, d_center_x, d_center_y ) ;
    fprintf(fp, "\n") ;
    
    fprintf(fp, "Nx_g = %d; Ny_g = %d; N_g = %d;\n", Nx_g, Ny_g, N_g) ;
    fprintf(fp, "\n") ;
    
    fprintf(fp, "i0_g = %d; i1_g = %d; j0_g = %d; j1_g = %d;\n", i0_g, i1_g, j0_g, j1_g) ;
    fprintf(fp, "i0_f = %d; i1_f = %d; j0_f = %d; j1_f = %d;\n", i0_f, i1_f, j0_f, j1_f) ;
    fprintf(fp, "i0_i = %d; i1_i = %d; j0_i = %d; j1_i = %d;\n", i0_i, i1_i, j0_i, j1_i) ;
	fprintf(fp, "\n") ;
	
	fprintf(fp, "flag_kernel = %d;\n", flag_kernel) ;	
	fprintf(fp, "D0 = %f; h0 = %f;\n", D0, h0) ;
	fprintf(fp, "weibull_lambda = %f; weibull_h = %f;\n", weibull_lambda, weibull_h) ;
	fprintf(fp, "lognorm_mu = %f; lognorm_s = %f; lognorm_Fr = %f;\n", lognorm_mu, lognorm_s, lognorm_Fr) ;
	fprintf(fp, "Mbeta_beta = %f; Mbeta_N = %d; Mbeta_D = %f; Mbeta_u = %f;\n", Mbeta_beta, Mbeta_N, Mbeta_D, Mbeta_u) ;
	fprintf(fp, "rho_amb = %f; T_amb = %f; c_pg = %f; g_acc = %f;\n", rho_amb, T_amb, c_pg, g_acc) ;
	fprintf(fp, "Io = %f; If = %f;\n", Io, If) ;
	fprintf(fp, "\n") ;
	
	//fprintf(fp, "dCsi = %f; NCsi = %d;\n", dCsi, NCsi) ;
	fprintf(fp, "Nr = %d; dr = %f; deltat_kernel = %f; Nquad = %d;\n", Nr, dr, deltat_kernel, Nquad) ;
	fprintf(fp, "\n") ;
	
	fprintf(fp, "flag_hotfront = %d; Nneighbour = %d;\n", flag_hotfront, Nneighbour) ;
	fprintf(fp, "tc_h_up = %f; tc_f_up = %f; tc_up = %f;\n", tc_h_up, tc_f_up, tc_up) ;
	fprintf(fp, "tc_h_dw = %f; tc_f_dw = %f; tc_dw = %f;\n", tc_h_dw, tc_f_dw, tc_dw) ;
	fprintf(fp, "\n") ;
	
	fprintf(fp, "flag_circle_init = %d; init_circle_N = %d; init_circle_inside = %d;\n", 
            flag_init_circle, init_circle_N, init_circle_inside) ;
	fprintf(fp, "init_circle_xc = %f; init_circle_xc = %f; init_circle_r = %f;\n", 
            init_circle_xc, init_circle_yc, init_circle_r) ;
	fprintf(fp, "flag_rectangle_init = %d; init_rectangle_N = %d; init_rectangle_inside = %d;\n", 
            flag_init_rectangle, init_rectangle_N, init_rectangle_inside) ;
	fprintf(fp, "init_rectangle_x0 = %f; init_rectangle_y0 = %f; init_rectangle_Lx = %f; init_rectangle_Ly = %f;\n", 
            init_rectangle_x0, init_rectangle_y0, init_rectangle_Lx, init_rectangle_Ly) ;
	fprintf(fp, "\n") ;
	
	fprintf(fp, "flag_fire_obstacle   = %d;\n", flag_fire_obstacle) ;
	fprintf(fp, "flag_fire_obstacle_1 = %d;\n", flag_fire_obstacle_1) ;
	fprintf(fp, "x0_fire_obstacle_1 = %f; x1_fire_obstacle_1 = %f;\n", x0_fire_obstacle_1, x1_fire_obstacle_1) ;
	fprintf(fp, "y0_fire_obstacle_1 = %f; y1_fire_obstacle_1 = %f;\n", y0_fire_obstacle_1, y1_fire_obstacle_1) ;
	fprintf(fp, "flag_fire_obstacle_2 = %d;\n", flag_fire_obstacle_2) ;
	fprintf(fp, "x0_fire_obstacle_2 = %f; x1_fire_obstacle_2 = %f;\n", x0_fire_obstacle_2, x1_fire_obstacle_2) ;
	fprintf(fp, "y0_fire_obstacle_2 = %f; y1_fire_obstacle_2 = %f;\n", y0_fire_obstacle_2, y1_fire_obstacle_2) ;
	fprintf(fp, "\n") ;
	
	fprintf(fp, "flag_saveformat_binary = %d;\n", flag_saveformat_binary) ;
	fprintf(fp, "flag_saveformat_text   = %d;\n", flag_saveformat_text) ;
	fprintf(fp, "flag_saveformat_ghostboxes = %d;\n", flag_saveformat_ghostboxes) ;
	fprintf(fp, "\n") ;
	
	fprintf(fp, "filename_ws = '%s'; save_each_iteration = %d;\n", filename_ws, save_each_iteration) ;
	fprintf(fp, "\n") ;
	
	fprintf(fp, "flag_save_data         = %d; filename_data         = '%s';\n", flag_save_data, filename_data) ;
	fprintf(fp, "flag_save_phi          = %d; filename_phi          = '%s';\n", flag_save_phi, filename_phi) ;
	fprintf(fp, "flag_save_phieff       = %d; filename_phieff       = '%s';\n", flag_save_phieff, filename_phieff) ;
	fprintf(fp, "flag_save_psi          = %d; filename_psi          = '%s';\n", flag_save_psi, filename_psi) ;
	fprintf(fp, "flag_save_tburn_phi    = %d; filename_tburn_phi    = '%s';\n", flag_save_tburn_phi, filename_tburn_phi) ;
	fprintf(fp, "flag_save_tburn_phieff = %d; filename_tburn_phieff = '%s';\n", flag_save_tburn_phieff, filename_tburn_phieff) ;
	fprintf(fp, "flag_save_tburn_psi    = %d; filename_tburn_psi    = '%s';\n", flag_save_tburn_psi, filename_tburn_psi) ;
	fprintf(fp, "flag_save_nx           = %d; filename_nx           = '%s';\n", flag_save_nx , filename_nx ) ;
	fprintf(fp, "flag_save_ny           = %d; filename_ny           = '%s';\n", flag_save_ny, filename_ny) ;
	fprintf(fp, "flag_save_kernel       = %d; filename_kernel       = '%s';\n", flag_save_kernel, filename_kernel) ;
	fprintf(fp, "flag_save_ros          = %d; filename_ros          = '%s';\n", flag_save_ros, filename_ros) ;
	fprintf(fp, "flag_save_fuelMap      = %d; filename_fuelMap      = '%s';\n", flag_save_fuelMap, filename_fuelMap) ;
	fprintf(fp, "flag_save_windMap      = %d; filename_windMap      = '%s';\n", flag_save_windMap, filename_windMap) ;
    fprintf(fp, "\n") ;
    
    fprintf(fp, "filename_extension_txt = '%s';\n",	filename_extension_txt) ;
    fprintf(fp, "filename_extension_bin = '%s';\n",	filename_extension_bin) ;
    fprintf(fp, "\n") ;
    
    fprintf(fp, "flag_byram = %d;\n", flag_byram) ;
    fprintf(fp, "LowHeatComb = %f; omega0 = %f; ros_byram = %f;\n", LowHeatComb, omega0, ros_byram) ;
	fprintf(fp, "\n") ;
	
	fprintf(fp, "flag_fireMap_const = %d; flag_ros_const = %d;\n", flag_fireMap_const, flag_ros_const) ;
	fprintf(fp, "Model = %zd;\n", Model) ;
	fprintf(fp, "flag_wind_const = %d; WindSpd = %f; Ut = %f; WindDir = %f;\n", 
             flag_wind_const, WindSpd, WindSpd/MS2FTMIN, WindDir) ;
	fprintf(fp, "Slope = %f; Aspect = %f;\n", Slope, Aspect) ;
	fprintf(fp, "M1 = %f; M10 = %f; M100 = %f;\n", M1, M10, M100) ;
	fprintf(fp, "Mherb = %f; Mwood = %f;", Mherb, Mwood) ;
    
    fclose(fp) ;
        
    return 0 ;
    
} /*** end of: save_datapre_to_file() *********************************/


    
    

/**********************************************************************
 * save array Kernel to file                                          *
 **********************************************************************/
int save_array_k_to_file(char* filename, int Niter)   
{
	char filename_full[MAX_CHAR_NFILE] ;
	FILE *fp ;
        
    // save data to binary file
    
    if (flag_saveformat_binary)
    {			
		if (Niter >= 0)
			sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename, 
				Niter, filename_extension_bin) ;
		else
			sprintf(filename_full, "%s_%s_end.%s", filename_ws, filename, 
				filename_extension_bin) ;
    
		fp = fopen(filename_full, "wb+") ;
		
		fwrite(nm_kernel, sizeof(int), 2, fp ) ;
		
		for (int i = 0; i < Nr; i++) 
			fwrite(r+i, sizeof(LSMLIB_REAL), 1, fp) ;
		
		for (int i = 0; i < Nr; i++) 
			fwrite(kernel+i, sizeof(LSMLIB_REAL), 1, fp) ;
		
		fclose(fp) ;
	}
	
	
	// save data to text file
	
	if (flag_saveformat_text)
	{
		if (Niter >= 0)
			sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename, 
				Niter, filename_extension_txt) ;
		else
			sprintf(filename_full, "%s_%s_end.%s", filename_ws, filename, 
				filename_extension_txt) ;
		
		fp = fopen(filename_full, "w") ;
		
		fprintf(fp, "%d %d\n", nm_kernel[0], nm_kernel[1] ) ;		
			
		for (int i = 0; i < nm_kernel[1]; i++) 
			fprintf(fp, "%e ", r[i]) ;
			
		fprintf(fp, "\n");
		
		for (int i = 0; i < nm_kernel[1]; i++) 
			fprintf(fp, "%e ", kernel[i]) ;

		fclose(fp) ;
	}
    
    return 0 ;
    
} /*** end of: save_array_k_to_file() *********************************/



/**********************************************************************
 * save array 2D to file                                              *
 **********************************************************************/
int save_array2d_to_file(char* filename, LSMLIB_REAL* data, int Niter, char type)   
{
	char filename_full[MAX_CHAR_NFILE] ;
	FILE *fp ;
    int flag_transpose ;
	
	// type "d": integer
	// type "f": double
	// type "t": size_t
        
    // save data to binary file
    
    if (flag_saveformat_binary)
    {			
		if (Niter >= 0)
			sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename, 
				Niter, filename_extension_bin) ;
		else
			sprintf(filename_full, "%s_%s_end.%s", filename_ws, filename, 
				filename_extension_bin) ;
    
		fp = fopen(filename_full, "wb+") ;
		
		fwrite(nm, sizeof(int), 2, fp ) ;
		
		for (int j = j0_save; j <= j1_save; j++) 
		{
			for (int i = i0_save; i <= i1_save; i++) 
			{
				int idx = i + j * Nx_g ;
				if (type=='d')
					fwrite(data+idx, sizeof(int), 1, fp) ;
				else if (type=='t')
					fwrite(data+idx, sizeof(size_t), 1, fp) ;
				else if (type=='f')
					fwrite(data+idx, sizeof(LSMLIB_REAL), 1, fp) ;
			}
		} 
		
		fclose(fp) ;
	}
	
	
	// save data to text file
	
	if (flag_saveformat_text)
	{
		if (Niter >= 0)
			sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename, 
				Niter, filename_extension_txt) ;
		else
			sprintf(filename_full, "%s_%s_end.%s", filename_ws, filename, 
				filename_extension_txt) ;
		
		fp = fopen(filename_full, "w") ;
		
		fprintf(fp, "%d %d\n", j1_save-j0_save+1, i1_save-i0_save+1 ) ;
		
		flag_transpose = 0 ;
        
        if (flag_transpose) {
            for (int i = i0_save; i < i1_save+ic; i++) 
            {
                for (int j = j0_save; j < j1_save+ic; j++) 
                {
                    int idx = i + j * Nx_g ;
                    
                    if (type=='d')
                        fprintf(fp, "%d ", (int)data[idx]) ;
                    else if (type=='t')
                        fprintf(fp, "%zu ", (size_t)data[idx]) ;
                    else if (type=='f')
                        fprintf(fp, "%13.6f ", data[idx]) ;
                }
                fprintf(fp, "\n") ;
            }
        }
        else {
            for (int j = j0_save; j <= j1_save; j++) 
            {
                for (int i = i0_save; i <= i1_save; i++) 
                {
                    int idx = i + j * Nx_g ;
                    if (type=='d')
                        fprintf(fp, "%d ", (int)data[idx]) ;
                    else if (type=='t')
                        fprintf(fp, "%zu ", (size_t)data[idx]) ;
                    else if (type=='f')
                        fprintf(fp, "%13.6f ", data[idx]) ;
                }
                fprintf(fp, "\n") ;
            }
        }
		
		fclose(fp) ;
	}
    
    return 0 ;
    
} /*** end of: save_array2d_to_file() *********************************/



/**********************************************************************
 * save data                                                          *
 **********************************************************************/
int save_dataiter_to_file(int Niter, LSMLIB_REAL t)
{
	char filename_full[MAX_CHAR_NFILE] ;
	int absNiter = abs(Niter) ;
	FILE *fp ;
	
	if (flag_saveformat_text)
	{
		//char *tmp = get_filename_full(filename_data, Niter, filename_extension_txt) ;
		//sprintf(filename_full, "%s", tmp) ;
		sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename_data, Niter, filename_extension_txt) ;
		fp = fopen(filename_full, "w") ;
		fprintf(fp, "%d %f\n", absNiter, t ) ;
		fclose(fp) ;
	}
	
	if (flag_saveformat_binary)
	{
		//char *tmp2 = get_filename_full(filename_data, Niter, filename_extension_bin) ;
		//sprintf(filename_full, "%s", tmp2) ;
		sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename_data, Niter, filename_extension_bin) ;
		fp = fopen(filename_full, "wb+") ;
		fwrite(&absNiter, sizeof(int), 1, fp ) ;
		fwrite(&t, sizeof(LSMLIB_REAL), 1, fp ) ;
		fclose(fp) ;
	}
		
	if (flag_save_phi) 
		save_array2d_to_file(filename_phi, dataLS->phi, Niter, 'f') ;
    
    if (flag_save_phieff) 
		save_array2d_to_file(filename_phieff, phieff, Niter, 'f') ;
    
    if (flag_save_psi) 
		save_array2d_to_file(filename_psi, psi, Niter, 'f') ;
		
	if (flag_save_tburn_phi) 
		save_array2d_to_file(filename_tburn_phi, tburn_phi, Niter, 'f') ;
		
	if (flag_save_tburn_phieff) 
		save_array2d_to_file(filename_tburn_phieff, tburn_phieff, Niter, 'f') ;
		
	if (flag_save_tburn_psi) 
		save_array2d_to_file(filename_tburn_psi, tburn_psi, Niter, 'f') ;
		
	if ((flag_save_ros && (!flag_ros_const)) || Niter==1) 
		save_array2d_to_file(filename_ros, ros, Niter, 'f') ;
    
    if (flag_save_nx) 
		save_array2d_to_file(filename_nx, nx, Niter, 'f') ;
    
    if (flag_save_ny) 
		save_array2d_to_file(filename_ny, ny, Niter, 'f') ;
    
    if (flag_save_kernel) 
		save_array_k_to_file(filename_kernel, Niter) ;
    
    return 0 ;
    
} /*** end of: save_dataiter_to_file() ********************************/



/**********************************************************************
 * Save last stuff to file                                            *
 **********************************************************************/
int save_datapost_to_file(int Niter, LSMLIB_REAL t)   
{
	FILE *fp ;
	char filename_full[MAX_CHAR_NFILE] ;
    
	sprintf(filename_full, "%s_%s_end.%s", filename_ws, filename_data, 
		filename_extension_txt) ;
    
    fp = fopen(filename_full, "w") ;
    
    fprintf(fp, "Niter = %d; t1 = %f;\n", Niter, t) ;
    
    fclose(fp) ;
        
    return 0 ;
    
} /*** end of: save_datapost_to_file() ********************************/



/**********************************************************************
 * get elapsed time                                                   *
 **********************************************************************/
int* get_hhmmss(double dt)
{
	int* hhmmss = (int*) calloc(4, sizeof(int)) ;
		
	int hh = (int)(dt/3600) ;
	int mm = (int)((dt-hh*3600)/60) ;
	int ss = (int)(dt-hh*3600-mm*60) ;
	int ms = (int)((dt-(double)(hh*3600+mm*60+ss))*1000.0) ;
	
	hhmmss[0] = hh ;
	hhmmss[1] = mm ;
	hhmmss[2] = ss ;
	hhmmss[3] = ms ;
	
	return hhmmss ;
	
} /*** end of: get_hhmmss() *******************************************/






/**********************************************************************
 * main                                                               *
 **********************************************************************/
int main(int argc, char **argv)
{
    int Niter ;
    unsigned int flag_iterate, flag_fire_at_border ;
    int *hhmmss, *whhmmss ;
    clock_t tic, toc ;    
    time_t  wtic, wtoc ;
    LSMLIB_REAL t, dt, dt_cfl, t_old ;
    //int c ;
    FILE *fp ;
    
    int Nprova, cprova, Nthreadsprova, i1, i2 ;
    
    hhmmss = (int*) calloc(4, sizeof(int)) ;
    
    // watch the clock
    tic = clock() ;
    wtic = time(NULL);
        

    // set here the kind of problem
    flag_adiff = 0 ;
    
	if (flag_adiff==0)
	{
	 //printf ("calling to read the input data\n") ;
         get_input_fire(argc, argv) ; printf ("input data read successfully\n"); }
    	else
        get_input_adiff(argc, argv) ;
	
        
    // preprocessing of input data
    
    preprocessing() ;

    // initializing OpenMP
    #ifdef _OPENMP
		#pragma omp parallel
		{
			if (omp_get_thread_num() == 0)
			{
				Nthreads_tot = omp_get_num_threads() ;
				
				if (Nthreads_max != 0) 
					Nthreads = Nthreads_max ;
				else
					Nthreads = Nthreads_tot ;
					
				printf("*** Executing LSfire+ on %d threads with OpenMP (%d threads available) ***\n", 
					Nthreads, Nthreads_tot) ;
			}
		}
		omp_set_num_threads(Nthreads) ;
	#else
		Nthreads_tot = 0 ;
		Nthreads = 0 ;
		printf("*** Executing LSfire+ on a single thread ***\n") ;
	#endif
	
	//Nprova = 80000 ; cprova = 0 ; i1 = 0 ; i2 = 0 ; Nthreadsprova = Nthreads ;
	//lsm2_test(&Nprova, &cprova, &i1, &i2, &Nthreadsprova) ;
    //printf("*** lsm2_test: %d, %d (%d, %d) %d ***\n", Nprova, cprova, i1, i2, Nthreadsprova) ;
	
	printf("*** Saving data to file '%s' ***\n", filename_ws) ;

    // initialization
    init_grid() ;

    init_lsm() ;

    init_fire() ;
    
    if (flag_kernel_Weibull || flag_kernel_Lognorm || flag_kernel_Lognorm_FF) 
		init_quad() ;
	
	if (flag_adiff)
		init_adiff() ;
			
	if (flag_kernel_Mainardi)
		init_anomalous() ;
	else
		Mbeta_beta = 0 ;
		
		
    if (flag_save_data) 
		save_datapre_to_file(filename_data) ;
		
	
	if (flag_save_fuelMap) 
	{
		if (flag_saveformat_text)
		{
			sprintf(filename_full, "%s_%s_000.%s", filename_ws, filename_fuelMap, 
				filename_extension_txt) ;
			
			fp = fopen(filename_full, "w") ;	
			
			fprintf(fp, "%d %d\n", j1_save-j0_save+1, i1_save-i0_save+1 ) ;
		
			
			for (int j = j0_save; j <= j1_save; j++) 
			{
				for (int i = i0_save; i <= i1_save; i++) 
				{
					int idx = i + j * Nx_g ;
					fprintf(fp, "%d ", (int)fuelMap[idx]) ;
				}
				fprintf(fp, "\n") ;
			}
		
			
			fclose(fp) ;
		}
		//save_array2d_to_file(filename_fuelMap, (LSMLIB_REAL*)fuelMap, 0, 'd') ;
		//printf("fuel: %zu, %d, %f\n", fuelMap[100], (int)fuelMap[100], (double)fuelMap[100]);
	}
		
		
    // time integration loop
    
    printf("Start time integration until t=%f\n", t_max) ;
    
    t  = 0.0 ; 

    dt = 0.0 ;
    flag_iterate = 1 ;
    Niter = 0 ;
    
    t_old = 0.0 ;
    
    if (flag_save_data) {
        save_dataiter_to_file(Niter, t) ;
        save_datapre_to_file(filename_data) ;
    }
    
	if (flag_kernel_Lognorm_FF || flag_kernel_Gauss_FF)
	{
		//printf ("\n calling read_firefront \n");
		read_firefront() ;
        	get_phieff(Niter, t_old, t_ff) ;
		//printf (" phieff calculated\n");
 		get_RoS(t);
        	dt_cfl = get_dt(0.5) ;
        	if (flag_deltat_const)
        	dt = min(deltat, dt_cfl) ;
        	else
        	dt = dt_cfl ;
       		//printf (" the cfl dt is = %f ",dt );
        	if (!flag_adiff) get_psi_FF(dt, t_ff) ;

	}
	else 
	{ 
	while (flag_iterate) 
    		{
        Niter++ ;

         //integration
        lsm2_RK2(dt, t) ;

        // update phi and calculate phieff, psi, tburn
        update_phi() ;
        //printf ("%f\n",If);
        get_phieff(Niter, t_old, t) ;

        if (!flag_adiff) get_psi(dt, t) ;
               
		if (flag_save_tburn) get_tburn(t) ;
		
		
		// watch the clock
        toc = clock() ;
        wtoc = time(NULL) ;
        hhmmss = get_hhmmss( (double)(toc-tic) / CLOCKS_PER_SEC ) ;
        whhmmss = get_hhmmss( (double)(wtoc-wtic) ) ;
        #ifdef _OPENMP 
			printf("[%02d:%02d:%02d|%d|%02d:%02d:%02d.%03d] #%03d, time: %f", 
				whhmmss[0], whhmmss[1], whhmmss[2], Nthreads,
				hhmmss[0], hhmmss[1], hhmmss[2], hhmmss[3], Niter, t) ;
		#else
			printf("[%02d:%02d:%02d.%03d] #%03d, time: %f", hhmmss[0], 
				hhmmss[1], hhmmss[2], hhmmss[3], Niter, t) ;
        #endif
        
		        
		// save data to file
		if (flag_save_data && ( Niter==1 || ( (save_each_iteration > 0) && 
			Niter%save_each_iteration==0 ) ) ) 
		{
			save_dataiter_to_file(Niter, t) ;
			save_datapost_to_file(Niter, t) ;
			printf(" [saved]") ;
		}
		printf("\n") ;
		
		
		// check if the fire has reached the border
		flag_fire_at_border = check_fire_border() ;
		
		// calculate the new time step
		dt_cfl = get_dt(0.5) ;
		if (flag_deltat_const)
			dt = min(deltat, dt_cfl) ;
		else
			dt = dt_cfl ;
				
		// decide if another iteration is in order
		if ( flag_fire_at_border )
		{
			printf("The fire has reached the border of the domain!\n") ;
			flag_iterate = 0 ;
		}
		if ( flag_t_end == -1 && t+dt > t_max )
		{
			flag_iterate = 0 ;
		}
		else if ( flag_t_end == 0 && t < t_max && t+dt > t_max )
		{
			t_old = t ;
			dt = t_max - t ;
			t = t_max ;
		}
		else if ( flag_t_end == 0 && t==t_max )
		{
			flag_iterate = 0 ;
		}
		else if ( flag_t_end == 1 && t > t_max )
		{
			flag_iterate = 0 ;
		}
		else
		{
			t_old = t ;
			t += dt ;
		}	
		
		
    }
    
    // save data to file
    if (flag_save_data && (save_each_iteration==0 || Niter%save_each_iteration))
		save_dataiter_to_file(Niter, t) ;
		
	if (flag_save_tburn)
	{       
		printf("...removing unnecessary files...") ;
		for (int n = 0; n < Niter; n++)
		{
			if (flag_saveformat_binary && flag_save_tburn_phi)
			{
				sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename_tburn_phi, 
					n, filename_extension_bin) ;
				remove(filename_full) ;
			}
			if (flag_saveformat_binary && flag_save_tburn_phieff)
			{
				sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename_tburn_phieff, 
					n, filename_extension_bin) ;
				remove(filename_full) ;
			}
			if (flag_saveformat_binary && flag_save_tburn_psi)
			{
				sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename_tburn_psi, 
					n, filename_extension_bin) ;
				remove(filename_full) ;
			}
			if (flag_saveformat_text && flag_save_tburn_phi)
			{
				sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename_tburn_phi, 
					n, filename_extension_txt) ;
				remove(filename_full) ;
			}
			if (flag_saveformat_text && flag_save_tburn_phieff)
			{
				sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename_tburn_phieff, 
					n, filename_extension_txt) ;
				remove(filename_full) ;
			}
			if (flag_saveformat_text && flag_save_tburn_psi)
			{
				sprintf(filename_full, "%s_%s_%03d.%s", filename_ws, filename_tburn_psi, 
					n, filename_extension_txt) ;
				remove(filename_full) ;
			}
		}
		//printf("done!\n") ;
	}
			
	save_datapost_to_file(Niter, t) ;
    
    if (flag_sanitycheck) 
		sanity_check(flag_sanitycheck, ros, nx, ny, nx, ny) ;
}
    free_all() ;
    
    return 0 ;
    
} /*** end of: main() *************************************************/
