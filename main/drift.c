#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "extras.h"
#include "random.h"
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_bessel.h>
#include <mpi.h>


double pi = 4.0*atan(1.);

double R = 1.;
double m, lambda, srml; /* srml = sqrt(m^2 - lambda^2) */
double lambdaextr[2];	/* lower and upper bound for lambda */
double rhoextr[2];		/* lower and upper bound for rho */
int maxrec = 100;		/* maximum recursion depth */
double s;


/*
 * Recursive auxiliary function for adaptiveSimpsons() function below
 */                                                                                                 
double adaptiveSimpsonsAux(double (*f)(double), double a, double b, double epsilon,                 
                         double S, double fa, double fb, double fc, int bottom) {                 
  double c = (a + b)/2, h = b - a;                                                                  
  double d = (a + c)/2, e = (c + b)/2;                                                              
  double fd = f(d), fe = f(e);                                                                      
  double Sleft = (h/12)*(fa + 4*fd + fc);                                                           
  double Sright = (h/12)*(fc + 4*fe + fb);                                                          
  double S2 = Sleft + Sright;                                                                       
  if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)   /* 15 from error analysis */
    return S2 + (S2 - S)/15;                                                                        
  return adaptiveSimpsonsAux(f, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +                    
         adaptiveSimpsonsAux(f, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);                     
}


/*
 * Adaptive Simpson's Rule
 */
double adaptiveSimpsons(double (*f)(double),   		/* ptr to function */
                           double a, double b,  	/* interval [a,b] */
                           double epsilon,			/* error tolerance */
                           int maxRecursionDepth) {	/* recursion cap */
  double c = (a + b)/2, h = b - a;                                                                  
  double fa = f(a), fb = f(b), fc = f(c);                                                           
  double S = (h/6)*(fa + 4*fc + fb);                                                                
  return adaptiveSimpsonsAux(f, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);                   
}                                                                                                   


/*
 *	Function whose first zero is the value of lambda
 */
double funcforlambda (double var){
	double u = sqrt(m*m - var*var);
	return  var*gsl_sf_bessel_J1(var*R) - u * gsl_sf_bessel_J0(var*R) * gsl_sf_bessel_K1_scaled(R*u)/gsl_sf_bessel_K0_scaled(R*u);
}



double f (double rho){
	if (rho == 0.)
		return 0;

	double aux = sqrt(m*m - (s+lambda)*(s+lambda));
	return rho*rho*gsl_sf_bessel_J0(lambda*rho)*gsl_sf_bessel_K0_scaled(aux*rho)*( lambda*gsl_sf_bessel_J1(lambda*rho)*gsl_sf_bessel_I0_scaled(aux*rho) + aux*gsl_sf_bessel_J0(lambda*rho)*gsl_sf_bessel_I1_scaled(aux*rho) );
}



/*******************************************************************************
 *
 *	MAIN ROUTINE
 *
 ******************************************************************************/
 
int main (int argc, char *argv[]) {

	if(argc < 2){
		printf("\nEnter average jump for the free process.\n\n");
		exit(EXIT_FAILURE);
	}
	
	double avjmp;
	double smin, smax, ds;
	double pref, C;
	
	char *out_name;
	FILE *out;
	
	avjmp = atof(argv[1]);
	m = 2./avjmp;

	out_name = malloc(100*sizeof(char));
	system("mkdir -p output/genf_z");
	sprintf(out_name, "output/genf_z/genf_%.3e.dat", avjmp);
	out = fopen(out_name,"w");
	
	lambdaextr[0] = 0.000001*m;
	lambdaextr[1] = 2.41/R; 	/* first zero of J0, plus a bit */
	if (m < lambdaextr[1])
		lambdaextr[1] = m - 1.e-10;
	lambda = Zbisection(funcforlambda, lambdaextr, 1.e-6);

	pref = 4.*m*m/R/(gsl_sf_bessel_J0(lambda*R)*gsl_sf_bessel_J0(lambda*R) + gsl_sf_bessel_J1(lambda*R)*gsl_sf_bessel_J1(lambda*R));
	
	smax = 1.;
	if(m-lambda < smax)
		smax = .1*(m - lambda);
	smin = -smax;
	ds = (smax-smin)/1000.;
	
	s = smin;
	while(s<smax){
		C = pref/(m*m - s*s - 2.*s*lambda);
		fprintf(out, "%e\t%.10lf\n", s, log(C*adaptiveSimpsons(f, 0., R, 1.e-10, maxrec)));
		s += ds;
 	}

 	fclose(out);
	
	exit(EXIT_SUCCESS);
	
}
