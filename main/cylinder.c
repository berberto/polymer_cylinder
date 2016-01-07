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


double pi = 4.0*atan(1.);

int Njumps = 1000;
double R = 1.;
double m, lambda, srml; /* srml = sqrt(m^2 - lambda^2) */
double lambdaextr[2];	/* lower and upper bound for lambda */
double rho, eta;		/* coordinates in the transverse plane (cylindrical) */
double phi, xi;			/* azimutal and (cosine of) zenithal angles for x'-x */
int maxrec = 100;		/* maximum recursion depth */
double minPhi;			/* minimum value of phi */
double minXi;			/* minimum value of cos(theta) */



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
	return  var*gsl_sf_bessel_J1(var*R) * gsl_sf_bessel_K0(R*sqrt(m*m - var*var)) - sqrt(m*m - var*var) * gsl_sf_bessel_J0(var*R) * gsl_sf_bessel_K1( R*sqrt(m*m - var*var));
}


/*
 *	Maximum transverse elongation (as function of rho, eta and phi)
 */
double B (double phi){
	return (sqrt(R*R - rho*rho*pow(sin(phi - eta), 2)) - rho*cos(phi-eta));
}



/*
 *  Distribution for the azimutal angle, phi
 */
double pdfPhi (double phi){
	return (.5/pi*(1. - B(phi)*srml*gsl_sf_bessel_K1( B(phi)*srml))/(1. - R*srml*gsl_sf_bessel_I0( rho*srml)*gsl_sf_bessel_K1( R*srml)));
}


/*
 * Distribution for the (cosine of the) zenital angle, xi
 */
double pdfXi (double xi){
	double alpha, norm, arstar;
	norm = 2./(m*m - lambda*lambda)*(1. - B(phi)*srml*gsl_sf_bessel_K1( B(phi)*srml));
	alpha = m - lambda*xi;

	if((xi == 1.)||(xi == -1.)) return 1./(alpha*alpha)/norm;
	
	else{
		arstar = alpha*B(phi)/sqrt(1. - xi*xi);
		return (1. - (1. + arstar)*exp(-arstar))/(alpha*alpha)/norm;
	}
}


/*
 *	Inversion of the cumulative distribution function:
 *	it gives x such that cdf(x) = w.
 *	
 *	pdf		: pointer to the pdf function (of one real variable)
 *	minX 	: minimum value of the R.V. X
 *	w		: real between 0 and 1
 *	deltaX	: accuracy in the determination of x = X(w) = G^{-1}(w)
 *
 */
double cdfInversion (double (*pdf)(double), double minX, double w, double deltaX) {
	double cdf, cdfnew;
	double startX;
	int i, j;
	
	cdfnew = 0.;
	cdf = 0.;
	i=-1;
	while((cdfnew - w)*(cdf-w) > 0) {
		i++;
		cdf = cdfnew;
		cdfnew += adaptiveSimpsons(pdf, minX + i*deltaX*100, minX + (i+1)*deltaX*100., .00001, maxrec);
	}
	j=-1;
	startX = minX + i*deltaX;
	while((cdfnew - w)*(cdf-w)>0) {
		j++;
		cdf = cdfnew;
		cdfnew = cdf + adaptiveSimpsons(pdf, startX + j*deltaX, startX + (j+1)*deltaX, .00001, maxrec);
	}
	return minX + deltaX*(i*100. + (j - .5));
}


/*
 *	Angle of a vector (x,y) with respect to the x axis
 */
double argument (double x, double y){
	if (x*x + y*y == 0)	return 0;
	
	else return carg(x + I*y);
}


void printlambert (double a, double b, int npoints){
	double x, delta;
	delta = (b - a)/npoints;
	for(x=a; x<= b; x +=delta)
		printf("%lf\t%lf\n", x, -gsl_sf_lambert_Wm1(-x));
}


int main (int argc, char *argv[]){
	if(argc < 2){
		printf("Set average jump length (in units of R)\n");
		exit(EXIT_SUCCESS);
	}
	
	char *out_name;
	FILE *out_traj, *out_phis, *out_xis, *out_rhos;
	
	int counter, i;
	double x0, y0, z0;
	double x, y, z;
	double xnew, ynew, znew, rhonew;
	double alpha, beta, arstar, r;

	double u[4]; /* 4 random numbers uniformly distributed in (0,1) */
	
	/* Initialization of the randomizer */
	srand(time(NULL));
	rlxd_init(1,rand());
	
	/*
	 *	Set constants
	 */
	m = 2./R/atof(argv[1]);
	lambdaextr[0] = 0.000001*m;
	lambdaextr[1] = 2.41/R; 	/* first zero of J0, plus a bit */
	if (m < lambdaextr[1]) lambdaextr[1] = m-.000001;
	
	lambda = Zbisection(funcforlambda, lambdaextr, .0001);
	
	srml = sqrt(m*m - lambda*lambda);
	
	printf("\nm=%lf\t\tlambda = %lf\n\n", m, lambda);
	
	/* Set extremes for phi and xi = cos(theta) */
	minPhi = -pi;
	minXi  = -1.;
	
	/*
	 *	Open output files
	 */
	out_name = malloc(100*sizeof(char));
	sprintf(out_name, "output/cylinder_%2.3lf.dat", atof(argv[1]));

	out_rhos = fopen("output/rhos.dat", "w");
	out_traj = fopen(out_name,"w");
	
	/*
	 *	Set initial point of the trajectory
	 */
	x0 = 0.;
	y0 = 0.;
	z0 = 0.;
	counter=0;
	
	x = x0; y=y0; z=z0;

	/* UNCOMMENT TO TEST DISTRIBUTION FOR PHI */
	/* out_phis = fopen("output/phis.dat", "w");
	fprintf(out_phis, "%lf\t%lf\t%lf\n", m, lambda, sqrt(x*x + y*y)); */
	
	/* UNCOMMENT TO TEST DISTRIBUTION FOR XI */
	phi = pi/4.;
	out_xis  = fopen("output/xis.dat", "w");
	fprintf(out_xis, "%lf\t%lf\t%lf\t%lf\n", m, lambda, sqrt(x*x + y*y), phi);
	
	while(counter<Njumps){
	
		if(counter%100==0) printf("%d\n",counter);
		
		rho = sqrt(x*x + y*y);
		fprintf(out_rhos, "%d\t%lf\n", counter, rho);
		eta = argument(x,y);	
	
		ranlxd(u, 4);
		
		/*
		 *	Generation of phi and xi with the inversion method (numerically)
		 */
		/* phi	= cdfInversion(pdfPhi, minPhi, u[0], 1.0e-4); */
		/* UNCOMMENT TO TEST DISTRIBUTION FOR PHI */
		/* fprintf(out_phis, "%lf\n", phi);
		counter++;
		continue; */
		
		xi	= cdfInversion(pdfXi, minXi, u[1], 1.0e-5);
		/* UNCOMMENT TO TEST DISTRIBUTION FOR XI */
		fprintf(out_xis, "%lf\n", xi);
		fflush(out_xis);
		counter++;
		continue;
				
		/*
		 *	Generation of r with inversion method via Lambert W function
		 */
		alpha = m - lambda*xi;
		if (xi == 1. || xi == -1.)
			beta = 1.;
		else {
			arstar = alpha*B(phi)/sqrt(1. - xi*xi);
			beta = 1. - (1. + arstar)*exp(-arstar);
		}
		r = (-gsl_sf_lambert_Wm1((beta*u[2] - 1.)*exp(-1.)) - 1.)/alpha;
		
		/*
		 *	Definition of x'
		 */
		xnew = x + r*sqrt(1-xi*xi)*cos(phi);
		ynew = y + r*sqrt(1-xi*xi)*sin(phi);
		znew = z + r*xi;
		
		/*
		 *	Rejection sampling
		 */
		rhonew = sqrt(xnew*xnew + ynew*ynew);
		if(u[3] <  gsl_sf_bessel_J0(lambda*rhonew)){
			x=xnew; y=ynew; z=znew;
			counter++;
			fprintf(out_traj, "%d\t%.5e\t%.5e\t%.5e\n", counter, x, y, z);
		}
	}
	
	/*fclose(out_phis);*/
	fclose(out_xis);
	fclose(out_traj);
	
	exit(EXIT_SUCCESS);
}
