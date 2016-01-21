#include <stdlib.h>
#include <stdio.h>
#include <math.h>


double fun(double x, double csi, double eta, double sigma, double s, double U, double theta)
{
	double e1, e2, e3 ;
	
	e1= -(pow(csi,2)+pow(eta,2)+exp(2*pow(2,0.5)*s*x)-2.0*exp(pow(2,0.5)*s*x)*(csi*sin(theta)+eta*cos(theta)))/(2.0*pow(sigma,2)) ;
	e2 = pow(2,0.5)*U*x/s ;
	e3 =-pow(U,2)/(2*pow(s,2)) ;
	
	return exp(e1+e2+e3) ;
}
 
int main(int argc, char **argv)
{
    int i, N ;
	char *filename_x, *filename_w ;
	double *x, *w ;
	double sigma, s, csi, eta, theta, U, Int ;
	
	FILE *fp_x, *fp_w ;
	
	csi = 2.5 ;
	eta = 1.5 ;
	sigma = 10.8 ;
	s = 10.4 ;
	U = 0.7 ;
	theta = 0.9 ;
	
	N = 100 ;
	x = (double*) calloc(N, sizeof(double)) ;
	w = (double*) calloc(N, sizeof(double)) ;
	
	filename_x = "quad_hermite_100_x.txt" ;
	filename_w = "quad_hermite_100_w.txt" ;
	
	fp_x = fopen(filename_x, "r") ;
	fp_w = fopen(filename_w, "r") ;
	
	Int = 0.0 ;
	for (i=0; i<N; i++)
	{
		fscanf(fp_x, "%lf", &x[i]) ;
		fscanf(fp_w, "%lf", &w[i]) ;
		//printf("%03d: %e, %e\n", i+1, x[i], w[i]) ;
		
		Int += w[i]*fun(x[i], csi, eta, sigma, s, U, theta) ;
	}
	
	Int /= (2.0*pow(M_PI,1.5)*pow(sigma,2)) ;
	
	printf("integrale: %e\n", Int) ;
	
	return 0 ;

}
