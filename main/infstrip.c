#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include "random.h"


double pi = 4.0*atan(1.);

int main (int argc, char *argv[]){

	int i;
	double x, y;
	double z;
	double B[2];
	char *out_name;
		out_name = malloc(100*sizeof(char));
	FILE *out_file;
	
	double L = 10;
	double N = 10000;
	double dt = .001;
	double D = 1.;
	double b = 1.;
	
	x = 0.;
	y = 0.;
	z = sin(pi*x/b);
	
	system("mkdir -p output");
	sprintf(out_name, "output/infstrip_D=%.3lf_b=%.3lf_dt=%.3lf_L=%.1lf",
						D, b, dt, L);
	out_file = fopen(out_name, "w");
	
	/* Initialization of the randomizer */
	srand(time(NULL));
	rlxd_init(1,rand());
	
	for (i=0; (x<2*L && x>-2*L) ; i++){
		gauss_dble(B,2);
		x += 2*D*pi/b*tanh(pi/b*x)*dt + sqrt(2*D*dt)*B[0];
		z += -3*D*pow(pi/b,2)*z*dt + sqrt(2*D*dt*(1-z*z))*pi/b*B[1];
		y = b/pi*asin(z);
		fprintf(out_file, "%lf\t%.4lf\t%.4lf\t%.4lf\n", i*dt, x, y, z);
	}
	fclose(out_file);
	
	
	exit(EXIT_SUCCESS);
}
