#include <unistd.h>
#include "LSfire+.h"

/**********************************************************************
 * read input data (fire model)                                       *
 **********************************************************************/
int get_input_fire(int argc, char **argv)
{
	
	int scaling ;
	double Ut ;
	
	int Mbeta_Nx, Mbeta_Ny, Mbeta_Nneighbour ;
	double Mbeta_dx ;
	
	/* 
     * Some of the input data can be provided  via command line options
     * as follows:
     * 
     * ./LSfire+.exe -k10 -s2 -D225 -T300 -t60 -I10000 -U6.7 -o0 -w1 -e10 -fLSout_k3_s2_D225_tf60_I10_U0670_o0 -nv1645 -t1.85
     * 
     * The previous command will run the code setting the following
     * values for some of the input data:
     * 
     * - flag_kernel = 10 (option -k10)
     * - scaling = 2 (option -s2)
     * - D0 = 225 (option -D225)
     * - t_max = 300 (option -T300)
     * - tc_dw = 60 (option -t60)
     * - Io = 10000 (option -I10000)
     * - Ut = 6.7 (option -U6.7)
     * - flag_fire_obstacle = 0 (option -o0)
     * - flag_wind_const = 1 (option -w1)
     * - save_each_iteration = 10 (option -e10)
     * - filename_ws = "LSout_k3_s2_D225_tf60_I10_U0670_o0" (option -fLSout_k3_s2_D225_tf60_I10_U0670_o0)
     * - number_of_vertices_from_FF = 1645 (option -n1645 )
     * - dt_time_step = 1.85 (option -d1.85)
     * When some (of all) the command line arguments are not provided,
     * the default values of the corresponding variables are used
     */ 
     
     
	/*******************************************************************
     * parsing of input arguments                                      *
     ******************************************************************/
	
    int i_flag_kernel = -1 ;
    int i_scaling     = -1 ;
    double i_D0       = -1 ;
    double i_t_max    = -1 ;
    double i_tc_dw    = -1 ;
    double i_Io       = -1 ;
    double i_Ut       = -1 ;
    double i_Mbeta_beta = -1 ;
    int i_flag_fire_obstacle = -1 ;
    int i_flag_wind_const = -1 ;
    int i_save_each_iteration = -9 ;
    char *i_filename_ws = NULL ;
    double i_dt		= -1;

    int index, cgetopt, optind, optopt ;  
    //opterr = 0 ;
    while ( (cgetopt = getopt(argc, argv, "k:s:D:T:t:I:U:B:o:w:e:f:x:")) != -1)
		switch (cgetopt)
		{
			case 'k': /* flag_kernel */
				i_flag_kernel = (int) strtol(optarg, (char**)NULL, 10) ;
				break ;
                
			case 's': /* scaling */
				i_scaling = (int) strtol(optarg, (char**)NULL, 10) ;
				break ;
			
            case 'D': /* D0 */
				i_D0 = (double) strtod(optarg, (char**)NULL) ;
				break ;
                
            case 'T': /* t_max */
				i_t_max = (double) strtod(optarg, (char**)NULL) ;
				break ;
                
			case 't': /* tc_dw */
				i_tc_dw = (double) strtod(optarg, (char**)NULL) ;
				break ;
                
			case 'I': /* Io */
				i_Io = (double) strtod(optarg, (char**)NULL) ;
				break ;
			
            case 'U': /* Ut */
				i_Ut = (double) strtod(optarg, (char**)NULL) ;
				break ;
                
			case 'B': /* Mbeta_beta */
				i_Mbeta_beta = (double) strtod(optarg, (char**)NULL) ;
				break ;
                
            case 'o': /* flag_fire_obstacle */
				i_flag_fire_obstacle = (int) strtol(optarg, (char**)NULL, 10) ;
				break ;
                
			case 'w': /* flag_wind_const */
				i_flag_wind_const = (int) strtol(optarg, (char**)NULL, 10) ;
				break ;
                
            case 'e': /* save_each_iteration */
				i_save_each_iteration = (int) strtol(optarg, (char**)NULL, 10) ;
				break ;
                
            case 'f': /* filename_ws */
				i_filename_ws = optarg ;
				break ;
	   case 'x': /* dt */
                                i_dt = (double) strtod(optarg, (char**)NULL) ;
                               break ;
                
            case '?':
				if (optopt == 'k')
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint(optopt))
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
				return 1 ;
                
			default:
				abort () ;
		}
			
		//for (index = optind; index < argc; index++)
		//	printf ("Non-option argument %s\n", argv[index]) ;        
       
	/*******************************************************************
     * END of parsing input arguments                                  *
     ******************************************************************/
	
	flag_sanitycheck = 0 ;
	flag_backward_compatibility = 0 ;
	
	flag_d_center = 0 ;
    
    //flag_adiff = 0 ;
	
	if (flag_adiff)
	{
		filename_ws = "Mout" ;
		Mbeta_beta = 2./4. ;
		Mbeta_u = 1.0 ;
		Mbeta_D = 225.0 ;
		
		if (i_Mbeta_beta>=0) printf("i_Mbeta_beta     : %f\n", Mbeta_beta = i_Mbeta_beta) ;
		
		scaling = 2 ;
		Mbeta_Nx = 100 ;
		Mbeta_Ny = 100 ;
		Mbeta_dx = 100.0 ;
		Mbeta_Nneighbour = 50 ;
	}
	else
	{
		scaling = 1 ;
		filename_ws = "LSout" ;
		Mbeta_Nx = Mbeta_Ny = Mbeta_dx = Mbeta_Nneighbour = 0 ;
	}
		
	if (i_filename_ws!=NULL) printf("i_filename_ws     : %s\n", filename_ws = i_filename_ws) ;
		
	Nthreads_max = 0 ; // number of threads (OpenMP)
	
    
	/****************************   space and time discretization   ***/
	
	if (i_scaling>=0) printf("i_scaling         : %d\n", scaling = i_scaling) ;
	
	if (flag_adiff)
	{
		Nx = Mbeta_Nx*scaling ; Ny = Mbeta_Ny*scaling ; // number of cells
		dx = Mbeta_dx/scaling ; dy = dx ; // size of the grid cell [ft]
	}
	else
	{
		//Nx = 200*scaling ; Ny = 100*scaling ; // number of cells
		Nx = 250*scaling ; Ny = 250*scaling ; // number of cells
		dx = 20.0*M2FT/scaling ; dy = dx ; // size of the grid cell [ft]
	}
	
	x_0 = 0.0 ; y_0 = 0.0 ; // off-set of the grid [ft]
		
	flag_t_end = 0 ; // [-1]: non supera t_max; [0]: finisce esattamente a t_max; [+1]: supera t_max
	t_max = 1000 ;  // max simulation time [min]
	if (i_t_max>=0) printf("i_t_max           : %f\n", t_max = i_t_max) ;
	
	flag_deltat_const = 0 ; // compute the solution each deltat
	
	deltat = 50 ; // (max) time step [min]
	
    
	/************************************   kernel discretization   ***/
	
	flag_kernel = 40 ;  // 0: LSM puro
						// 10: Gaussian (analytical)
						// 11: Gaussian (numerical)
						// 20: convolution, q: Weibull
						// 30: convolution, q: lognorm 
						// 31: convolution, q: lognorm (Fr<1)
						// 32: convolution, q: lognorm (Fr>1)
						// 40: anomalous diffusion (Mv)
	
	if (i_flag_kernel>=0) printf("i_flag_kernel     : %d\n", flag_kernel = i_flag_kernel) ;
	                   
                       
	/******************************   randomized-model parameters   ***/
	
	Nquad = 20 ;
	dr = 10.0 ; deltat_kernel = 1 ;
	
	// parameters for Gaussian distribution [sigma^2=2*D*t]
	D0 = 100.0 ; // diffusion coefficient [ft^2/min]
	//D0 = 75.0 ; // diffusion coefficient [ft^2/min]
	//D0 = 2.55 ; // diffusion coefficient [ft^2/min]
	if (flag_adiff) D0 = Mbeta_D ;
	if (i_D0>=0) printf("i_D0              : %f\n", D0 = i_D0) ;
			 
	h0 = 2.0*sqrt(D0) ; //15.0; 30.0 ; // smoothing length [ft/min^0.5]
	Nneighbour = 4*scaling ;
	//Nneighbour = 4*scaling ;
	//Nneighbour = 6*scaling ;
	if (flag_adiff) Nneighbour = Mbeta_Nneighbour*scaling ;
	
	// parameters for Weibull distribution 
	weibull_lambda = 100.0 ;
	weibull_h = 2.0 ;
	
	// parameter for Mainardi function (anomalous diffusion)
	if (!flag_adiff) Mbeta_beta = 1./2. ;
	
	// surface fireline intensity (for Byram formula and for If)
	Io = 10000.0 ; // [kW/m]	
	if (i_Io>=0) printf("i_Io              : %f\n", Io = i_Io) ;
	
	// parameter for Log-normal distribution
	If = Io + 0.015 ; // fire intensity [kW/m]
	
	// parameters for calculation of the Freude number (Fr)	
	rho_amb = 1.1    ; // ambient gas density [kg/m^3]
    	T_amb   = 300.0  ; // ambient temperature [K] 
    	c_pg    = 1121.0 ; // mean cp of the gas [kJ/(kg*K)]
    	g_acc   = 9.81   ; // gravity acceleration [m/s^2]
	
	
	// parameters from FireFront
	//ntotal = i_nv; // number of burnt points of forefire contour
	t_ff = i_dt * S2MIN; // time step in forefire [min]
	/******************   heating-before-burning model parameters   ***/
	
	
	flag_hotfront = 1 ; // activate (1) or deactivate (0)
	// sopravento (upwind)
	tc_h_up = 20.0 ; // ignition delay due to hot air [min]
	tc_f_up = 20.0 ; // ignition delay due to firebrands [min]
	// sottovento (downwind)
	tc_h_dw = 20.0 ; // ignition delay due to hot air [min]
	tc_f_dw = 20.0 ; // ignition delay due to firebrands [min] ***
	
	tc_up = 10.0 ; // SOPRAVENTO //tc_h_up*tc_f_up/(tc_h_up+tc_f_up) ; // time scale fuel [min]
	tc_dw =  1.0 ; // SOTTOVENTO //tc_h_dw*tc_f_dw/(tc_h_dw+tc_f_dw) ; // time scale fuel [min]
	
	if (i_tc_dw>=0) printf("i_tc_dw           : %f\n", tc_dw = i_tc_dw) ;
		
        
	/*********************************************   initial data   ***/
	
	flag_init_circle = 1 ;
	init_circle_N = 1 ;
	//init_circle_xc = x_0+0.5*Nx*dx ; init_circle_yc = y_0+0.5*Ny*dy ; init_circle_r = 500.0 ; 
	init_circle_xc = 2500.0*M2FT ; init_circle_yc = 2500.0*M2FT ; init_circle_r = 300.0*M2FT ; 
	init_circle_inside = -1 ; // +1 (-1) -> phi > 0 inside (outside) the circle
	
	flag_init_rectangle = 0 ;
	init_rectangle_N = 1 ;
	init_rectangle_Lx = 1000.0 ; init_rectangle_Ly = 1000.0 ;
	init_rectangle_x0 = (x_0+0.5*Nx*dx) - 0.5*init_rectangle_Lx ;
	init_rectangle_y0 = (y_0+0.5*Ny*dy) - 0.5*init_rectangle_Ly ;
	init_rectangle_inside = -1 ;
	
    
	/*********************************   presence of the obstacle   ***/
	
	flag_fire_obstacle_1 = 1 ;
	flag_fire_obstacle_2 = 1 ;
	flag_fire_obstacle = 2 ;
	
	if (i_flag_fire_obstacle>=0)
	{
		if (i_flag_fire_obstacle == 0)
		{	flag_fire_obstacle_1 = flag_fire_obstacle_2 = 0 ; }
		else if (i_flag_fire_obstacle == 1)
		{	flag_fire_obstacle_1 = 1 ; flag_fire_obstacle_2 = 0 ; }
		else if (i_flag_fire_obstacle == 2)
		{	flag_fire_obstacle_1 = flag_fire_obstacle_2 = 1 ; }
		printf("i_flag_ob         : %d, %d [%d]\n", flag_fire_obstacle_1, 
            flag_fire_obstacle_2, flag_fire_obstacle = i_flag_fire_obstacle) ;
	}
		
	//x0_fire_obstacle_1 = 3100.0*M2FT      ; x1_fire_obstacle_1 = x0_fire_obstacle_1 + 60.0*M2FT ;
	x0_fire_obstacle_1 = -1.0*INFINITY      ; x1_fire_obstacle_1 = INFINITY ;
    //x0_fire_obstacle_1 = 3100.0        ; x1_fire_obstacle_1 = x0_fire_obstacle_1 + 200.0 ;
	//y0_fire_obstacle_1 = -1.0*INFINITY ; y1_fire_obstacle_1 = INFINITY ;
	//y0_fire_obstacle_1 = 3100.0*M2FT ; y1_fire_obstacle_1 = y0_fire_obstacle_1 + 40.0 * M2FT ;
	y0_fire_obstacle_1 = 3100.0*M2FT ; y1_fire_obstacle_1 = y0_fire_obstacle_1 + 60.0 * M2FT ;
		
	//x0_fire_obstacle_2 = 6000.0        ; x1_fire_obstacle_2 = 6300.0 ;
	x0_fire_obstacle_2 = -1.0*INFINITY      ; x1_fire_obstacle_2 = INFINITY ;
    //x0_fire_obstacle_2 = 2000.0+6000.0        ; x1_fire_obstacle_2 = x0_fire_obstacle_2 + 350.0 ;
	//y0_fire_obstacle_2 = -1.0*INFINITY ; y1_fire_obstacle_2 = INFINITY ;
	y0_fire_obstacle_2 = 3600.0*M2FT ; y1_fire_obstacle_2 = y0_fire_obstacle_2 + 80.0 * M2FT ;
	
	
	/***************************************   FireLib input data   ***/
	
	flag_byram = 1 ; // 1: set ros_max (i.e. the ros when u*n=0) equal to ros_Byram
	LowHeatComb = 22000.0 ; // fuel low heat of combustion [kJ/kg]
	omega0 = 2.243 ; // oven-dry mass of fuel consumed per unit area [kg/m^2] 
	ros_byram = Io / (LowHeatComb*omega0) * MS2FTMIN ; // ros [ft/min]
	
	flag_fireMap_const = 1 ;
	
	Ut = 6.7 ;
	if (i_Ut>=0) printf("i_Ut              : %f\n", Ut = i_Ut) ;
    
    flag_wind_const = 1 ;
    if (i_flag_wind_const>=0) printf("i_flag_wind_const : %d\n", flag_wind_const = i_flag_wind_const) ;

	Model    = 9  ;  // NFFL model number
	WindSpd  = Ut * MS2FTMIN ; // wind speed [ft/min]
	WindDir  = 0.0000001 ; // wind angle, clockwise from north [degrees]
	Slope    = 0.0 ; // fraction rise/reach [-]
	Aspect   = 0.0 ; // clockwise from north [degrees] 
	M1       = 0.1 ; // 1-hr dead fuel moisture 
	M10      = 0.0 ; // 10-hr dead fuel moisture
	M100     = 0.0 ; // 100-hr dead fuel moisture
	Mherb    = 0.0 ; // live herbaceous fuel moisture 
	Mwood    = 0.0 ; // live woody fuel moisture
   
	
	/**********************************************   data saving   ***/
		
	flag_saveformat_binary = 1 ; // save data in binary format
	flag_saveformat_text   = 0 ; // save data in text format
	flag_saveformat_ghostboxes = 1 ;  // save arrays including ghostboxes
	
	save_each_iteration = 1 ; // 0: salva solo alla fine
	
	if (i_save_each_iteration!=-9) printf("i_save_each_iteration: %d\n", 
        save_each_iteration = i_save_each_iteration) ;
   	
    printf ("ok\n"); 
    flag_save_data         = 1 ; filename_data         = "data" ;    
    flag_save_phi          = 1 ; filename_phi          = "phi" ;
    flag_save_phieff       = 1 ; filename_phieff       = "phieff" ;
    flag_save_psi          = 1 ; filename_psi          = "psi" ; 
    flag_save_tburn_phi    = 1 ; filename_tburn_phi    = "tburn_phi" ;
    flag_save_tburn_phieff = 1 ; filename_tburn_phieff = "tburn_phieff" ;
    flag_save_tburn_psi    = 1 ; filename_tburn_psi    = "tburn_psi" ;
    flag_save_nx           = 0 ; filename_nx           = "nx" ;
    flag_save_ny           = 0 ; filename_ny           = "ny" ;
    flag_save_kernel       = 0 ; filename_kernel       = "kernel" ;
    flag_save_ros          = 0 ; filename_ros          = "ros" ;
    flag_save_fuelMap      = 0 ; filename_fuelMap      = "fuel" ;
    flag_save_windMap      = 0 ; filename_windMap      = "wind" ;
    
    filename_extension_txt = "txt" ;
	filename_extension_bin = "dat" ;
	
	return 0 ;
	
} /*** end of: get_input_fire() ***************************************/



/***********************************************************************
 * Wind speed and direction as functions of time and position          *
 * (only used when flag_wind_const==1)                                 *
 **********************************************************************/ 

double t1 = 90.0 ;
double t2 = 300.0 ;
double z_start = 90.0 ;
double z_final = 30.0 ;



/* Wind direction                                                     *
 * This is the angle of the wind (in degrees, clockwise from north)   */
double fun_WindDir(LSMLIB_REAL t, LSMLIB_REAL x, LSMLIB_REAL y)
{
    double z ; /* [degrees], clockwise from north */
    
    if (t <= t1 )
        z = z_start ;
    else if (t <= t2 )
        z = z_start + (z_final-z_start)*(t-t1)/(t2-t1) ;
    else
        z = z_final ;
            
    return z ;
}

/* Derivative of the wind direction                                   *
 * This is the derivative of the angle of the wind (in degrees,       *
 * clockwise from nort) with respect to time (in minutes)             */
double fun_WindDir_der(LSMLIB_REAL t, LSMLIB_REAL x, LSMLIB_REAL y)
{    
    double dz_dt ; /* [degrees/min] */
        
    if (t <= t1 )
        dz_dt = 0.0 ;
    else if (t <= t2 )
        dz_dt = (z_final-z_start)/(t2-t1) ;
    else
        dz_dt = 0.0 ;
 
    return dz_dt ;
}

/* Wind speed                                                         *
 * This is the wind speed (in ft/min)                                 */ 
double fun_WindSpd(LSMLIB_REAL t, LSMLIB_REAL x, LSMLIB_REAL y, LSMLIB_REAL WindSpd_max)
{    
    double WindSpd ; /* [ft/min] */
    
    WindSpd = WindSpd_max ; 
 
    return WindSpd ;
}
