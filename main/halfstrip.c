#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include "random.h"


double pi = 4.0*atan(1.);

double L = 10;
double N = 10000;
double dt = .001;
double D = .1;
double b = 1.;



double complex inverse (double complex z){
	return b/pi*cacosh((1.+z)/(1.-z));
}

double complex transf (double complex w){
	return ((ccosh(pi/b*w) - 1.)/(ccosh(pi/b*w) + 1.));
}

double complex jacob (double complex z){
	return pi/b*(1.-z)*csqrt(z);
}

int main (int argc, char *argv[]){
	
	int i;
	double complex w, z;	/* values of the processes */
	double complex jac;		/* jacobian of the transformation */
	double dtau;			/* rescaled time step */
	double B[2];			/* random normal 2D vector (Re and Im parts) */

	char *out_name;
		out_name = malloc(100*sizeof(char));
	FILE *out_file;	

	if (argc < 3) exit(EXIT_FAILURE); 
		w = atof(argv[1]) + I*atof(argv[2]);
/*		z = transf(w);
	printf("w = %.2lf + I %.2lf,\t arg(w) = %.3lf\n", creal(w), cimag(w), carg(w)/pi);
	printf("z = %.2lf + I %.2lf,\t arg(z) = %.3lf\n", creal(z), cimag(z), carg(z)/pi);
	
	TEST TRASFORMAZIONE	
	printf("transf :\tz = %.3lf + I %.3lf\n",creal(z), cimag(z));
	printf("inverse:\tw = %.3lf + I %.3lf\n",creal(inverse(z)), cimag(inverse(z)));
	printf("jacob  :\tz'= %.3lf + I %.3lf\n",creal(jacob(z)), cimag(jacob(z)));
*/
	
	/* initial value of the process in the half strip geometry */
/*	w = 1. + .5*I; */
	z = transf(w);
	printf("w = %.2lf + I %.2lf,\t arg(w) = %.3lf\n", creal(w), cimag(w), carg(w)/pi);
	printf("z = %.2lf + I %.2lf,\t arg(z) = %.3lf\n", creal(z), cimag(z), carg(z)/pi);
	printf("jac(z) = %.3lf + I %.3lf\n", jacob(z));
	printf("arg(w)= %.3lf\targ(z) = %.3lf\n", carg(w)/pi, carg(z)/pi);
	
	system("mkdir -p output");
	sprintf(out_name, "output/halfstrip_D=%.3lf_b=%.3lf_dt=%.3lf_L=%.1lf",
						D, b, dt, L);
	out_file = fopen(out_name, "w");
	fprintf(out_file, "#\n# time\t\tRe & Im in UHP\t\tRe & Im in half strip\n#\n");
	
	/* Initialization of the randomizer */
	srand(time(NULL));
	rlxd_init(1,rand());
	
	
	for (i=0; i<N; i++){
	
		gauss_dble(B,2);
		fprintf(out_file, "%.4lf\t%.5e\t%.5e\t%.5e\t%.5e\t\t%lf + I %lf\n", \
				i*dt, creal(z), cimag(z), creal(w), cimag(w), creal(-I/(2.*carg(z)*z)), cimag(-I/(2.*carg(z)*z)));
		jac = jacob(z);
		dtau = cabs(jac)*cabs(jac)*dt; /* con o senza cabs(jac)*cabs(jac) */
		z = z + 4*D*(-I*1.0/(2.*carg(z)*(-creal(z)+I*cimag(z))))*dtau	 \
								+ sqrt(2*D*dtau)*(B[0]+I*B[1]);
		if(cimag(z)<0) break;
		
		w = inverse(z);
	}

	fclose(out_file);

	
	exit(EXIT_SUCCESS);
}
