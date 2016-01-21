#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char **argv)
{
    int i, N ;
	char *filename ;
	double *x ;
	
	FILE *fp ;
	
	N = 100 ;
	x = (double*) calloc(N, sizeof(double)) ;
	
	filename = "glglib/gl_alpha1_A0_B1_N100_w.txt" ;
	
	fp = fopen(filename, "r") ;
	
	for (i=0; i<N; i++)
	{
		fscanf(fp, "%lf", &x[i]) ;
		printf("%02d: %e\n", i+1, x[i]) ;
	}
	
	return 0 ;

}	
