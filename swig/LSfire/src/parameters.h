/***********************************************************************
 * Some parameters and conversion factors                              *
 **********************************************************************/
 
#define MAX_CHAR_NFILE (60)  // max length of filename

#define NMAX_CATALOG (13) // max number of fuel items

#define MAX_LAST_GRAD (1.0e-12) // parameter for the Gauss_n subroutine

#define T_NO_BURN        (-1.0) // time of burning for unburnt region
#define PHI_THRESHOLD    (0.5)  // burnt region threshold
#define PHIEFF_THRESHOLD (0.5)  // burnt region threshold
#define PSI_THRESHOLD    (1.0)  // burnt region threshold

#define DEG2RAD (M_PI/180.0)     // degrees->radiants
	
/* I suppose that the formulas by (Sardoy et al., 2008) implemented
 * in the subroutine 'get_lognorm_par' provide mu and s for 
 * q(u,t) = 1/(sqrt(2*pi)*s*u) exp(-(ln(u)-mu)^2/(2*s^2))     [*]
 * when u is given in meters. So, if u is given in ft (as in the 
 * present case), the formula [*] must be modified as follows:
 * - in the argument of the exponential:
 *     ln(u*ft2m)-mu = ln(u)+ln(ft2m)-mu = ln(u)-(mu-ln(ft2m))
 * - in the denominator:
 *     s*(u*ft2m) = (s*ft2m)*u
 *   being:
 *     ft2m = 1/3.2808399 = 0.304799999536704
 *     ln(ft2m) = ln(1/3.2808399) = -1.1880994566896459
 * NB: When integrating, it is now necessary to multiply by ft2m*u
 *     instead of u, so the term ft2m in the denominator is canceled
 * The only modification is thus:
 *     mu -> mu+ln(m2ft) = mu+1.1880994566896459
 */
 
#define M2FT    (3.2808399)       // 1 m = 3.2808399 ft
#define MPH2MS  (0.44704)         // 1 mph = 0.44704 m/s 
#define MS2MPH  (2.2369362920544) // 1 m/s = 2.2369362920544 mph
#define MS2FTMIN (196.850393701)  // 1 m/s = 196.850394 ft/min
#define FTMIN2MS (0.00508)        // 1 ft/min = 0.00508 m/s
 
#define LN_M2FT (1.1880994566896459)  // log(m2ft)=ln(3.2808399)	
#define FT2M (0.3048) // 1 ft = 0.3048 m
#define S2MIN (0.0166667) // 1 sec = 0.0166667min 
