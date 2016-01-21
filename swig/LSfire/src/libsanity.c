#include "LSfire+.h"

/**********************************************************************
 * sanity check (temporary function for development)                  *
 **********************************************************************/
int sanity_check(unsigned int flag_sanity_check,  LSMLIB_REAL* rosMap,
	LSMLIB_REAL* nxMap, LSMLIB_REAL* nyMap, LSMLIB_REAL* timeMap, 
	LSMLIB_REAL* gptimeMap)
{
	int i, j, k ;
    LSMLIB_REAL phi_min, rosMap_min, timeMap_min, gptimeMap_min ;
    LSMLIB_REAL phieff_min, psi_min, nxMap_min, nyMap_min ;
	LSMLIB_REAL phi_max, rosMap_max, timeMap_max, gptimeMap_max ;
	LSMLIB_REAL phieff_max, psi_max, nxMap_max, nyMap_max ;
	LSMLIB_REAL phi_tot, rosMap_tot, timeMap_tot, gptimeMap_tot ;
	LSMLIB_REAL phieff_tot, psi_tot, nxMap_tot, nyMap_tot ;
	
	//N_gb  = gridLS->num_gridpts ;
		
    phi_min = INFINITY ; phi_max = -1.0*INFINITY ; phi_tot = 0.0 ;
    rosMap_min = INFINITY ; rosMap_max = -1.0*INFINITY ; rosMap_tot = 0.0 ;
    timeMap_min = INFINITY ; timeMap_max = -1.0*INFINITY ; timeMap_tot = 0.0 ;
    gptimeMap_min = INFINITY ; gptimeMap_max = -1.0*INFINITY ; gptimeMap_tot = 0.0 ;
    phieff_min = INFINITY ; phieff_max = -1.0*INFINITY ; phieff_tot = 0.0 ;
    psi_min = INFINITY ; psi_max = -1.0*INFINITY ; psi_tot = 0.0 ;
    nxMap_min = INFINITY ; nxMap_max = -1.0*INFINITY ; nxMap_tot = 0.0 ;
    nyMap_min = INFINITY ; nyMap_max = -1.0*INFINITY ; nyMap_tot = 0.0 ;
    
    for (j = gridLS->jlo_gb; j <= gridLS->jhi_gb; j++) 
    {
        for (i = gridLS->ilo_gb; i <= gridLS->ihi_gb; i++) 
        {
			k = i + j * gridLS->grid_dims_ghostbox[0];
			
			if (dataLS->phi[k] < phi_min) phi_min = dataLS->phi[k] ;
			if (rosMap[k] < rosMap_min) rosMap_min = rosMap[k] ;
			if (timeMap[k] < timeMap_min) timeMap_min = timeMap[k] ;
			if (gptimeMap[k] < gptimeMap_min) gptimeMap_min = gptimeMap[k] ;
			if (phieff[k] < phieff_min) phieff_min = phieff[k] ;
			if (psi[k] < psi_min) psi_min = psi[k] ;
			if (nxMap[k] < nxMap_min) nxMap_min = nxMap[k] ;
			if (nyMap[k] < nyMap_min) nyMap_min = nyMap[k] ;
			
			if (dataLS->phi[k] > phi_max) phi_max = dataLS->phi[k] ;
			if (rosMap[k] > rosMap_max) rosMap_max = rosMap[k] ;
			if (timeMap[k] > timeMap_max) timeMap_max = timeMap[k] ;
			if (gptimeMap[k] > gptimeMap_max) gptimeMap_max = gptimeMap[k] ;
			if (phieff[k] > phieff_max) phieff_max = phieff[k] ;
			if (psi[k] > psi_max) psi_max = psi[k] ;
			if (nxMap[k] > nxMap_max) nxMap_max = nxMap[k] ;
			if (nyMap[k] > nyMap_max) nyMap_max = nyMap[k] ;
			
			phi_tot       += dataLS->phi[k] ;
			rosMap_tot    += rosMap[k] ;
			timeMap_tot   += timeMap[k] ;
			gptimeMap_tot += gptimeMap[k] ;
			phieff_tot    += phieff[k] ;
			psi_tot       += psi[k] ;
			nxMap_tot     += nxMap[k] ;
			nyMap_tot     += nyMap[k] ;
		}
	}	    
	
	if (flag_sanity_check == 2)
    {
		printf("Sanity Check...\n") ;
		printf("...indexes: %d, %d, %d, %d, %d, %d\n", gridLS->jlo_gb, gridLS->jhi_gb, 
			gridLS->ilo_gb, gridLS->ihi_gb, gridLS->grid_dims_ghostbox[0], k) ;
		printf("...phi (min, max, tot): %f, %f, %f\n", phi_min, phi_max, phi_tot) ;
		printf("...rosMap  (min, max, tot): %f, %f, %f\n", rosMap_min, rosMap_max, rosMap_tot) ;
		printf("...timeMap (min, max, tot): %f, %f, %f\n", timeMap_min, timeMap_max, timeMap_tot) ;
		printf("...gptimeMap (min, max, tot): %f, %f, %f\n", gptimeMap_min, gptimeMap_max, gptimeMap_tot) ;
		printf("...phieff (min, max, tot): %f, %f, %f\n", phieff_min, phieff_max, phieff_tot) ;
		printf("...psi (min, max, tot): %f, %f, %f\n", psi_min, psi_max, psi_tot) ;
		printf("...nxMap (min, max, tot): %f, %f, %f\n", nxMap_min, nxMap_max, nxMap_tot) ;
		printf("...nyMap (min, max, tot): %f, %f, %f\n", nyMap_min, nyMap_max, nyMap_tot) ;
	}
	
	if (flag_sanity_check > 0)
	{   
		printf("Sanity check grand total: %f\n", phi_min+phi_max+phi_tot+rosMap_min+
			rosMap_max+rosMap_tot+timeMap_min+timeMap_max+timeMap_tot+
			gptimeMap_min+gptimeMap_max+gptimeMap_tot+phieff_min+
			phieff_max+phieff_tot+psi_min+psi_max+psi_tot+nxMap_min+
			nxMap_max+nxMap_tot+nyMap_min+nyMap_max+nyMap_tot) ;
	}
    
    return 0 ; 
    
} /*** end of: sanity_check() *****************************************/
