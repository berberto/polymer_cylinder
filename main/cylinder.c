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

double R = 1.;
double m, lambda, srml; /* srml = sqrt(m^2 - lambda^2) */
double lambdaextr[2];	/* lower and upper bound for lambda */
double rhoextr[2];		/* lower and upper bound for rho */
double rho, eta;		/* coordinates in the transverse plane (cylindrical) */
double phi, xi;			/* azimutal and (cosine of) zenithal angles for x'-x */
int maxrec = 100;		/* maximum recursion depth */
double minPhi;			/* minimum value of phi */
double minXi;			/* minimum value of cos(theta) */
double normPhi, normXi, normRho, b;	/* normalization constants */
double w;



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
	return (1. - B(phi)*srml*gsl_sf_bessel_K1(B(phi)*srml))*normPhi;
}


/*
 * Distribution for the (cosine of the) zenital angle, xi
 */
double pdfXi (double xi){
	double alpha, arstar;
	alpha = m - lambda*xi;

	if((xi >= 1.)||(xi <= -1.)) return 1./(alpha*alpha)/normXi;
	
	else{
		arstar = alpha*b/sqrt(1. - xi*xi);
		return (1. - (1. + arstar)*exp(-arstar))/(alpha*alpha)/normXi;
	}
}


/*
 *	Distriburion of the Tanh(xi)
 */
double pdfAtanhXi (double u) {
	double aux;
	aux = b*srml*cosh(u);
	return b*b*(1. - (1. + aux)*exp(-aux))/(aux*aux)/normXi;
}


/*
 *	Stationary distribution for rho
 */
double pdfRho (double r) {
	return lambda * r * gsl_sf_bessel_J0(lambda*r) / normRho;
}


/*
 *	Function whose zero is the the rho sampled according its stationary distribution
 */
double funcforrho (double r) {
	return (r * gsl_sf_bessel_J1(lambda*r)/normRho - w);
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
double cdfInversion (double (*pdf)(double), double minX, double w, double deltaX_rough, double deltaX_fine) {
	double cdf, cdfnew;
	double startX;
	int i, j;
	
	cdfnew = 0.;
	cdf = 0.;
	i=-1;
	while((cdfnew - w)*(cdf-w) > 0) {
		i++;
		cdf = cdfnew;
		cdfnew += adaptiveSimpsons(pdf, minX + i*deltaX_rough, minX + (i+1)*deltaX_rough, 1.e-8, maxrec);
	}
	j=-1;
	startX = minX + i*deltaX_rough;
	while((cdfnew - w)*(cdf-w)>0) {
		j++;
		cdf = cdfnew;
		cdfnew = cdf + adaptiveSimpsons(pdf, startX + j*deltaX_fine, startX + (j+1)*deltaX_fine, 1.e-8, maxrec);
	}
	return minX + i*deltaX_rough + (j-.5)*deltaX_fine;
}


/*
 *	Angle of a vector (x,y) with respect to the x axis (from -pi to pi)
 */
double argument (double x, double y){
	if (x*x + y*y == 0)	return 0;
	
	else return carg(x + I*y);
}



/*******************************************************************************
 *
 *	MAIN ROUTINE
 *
 ******************************************************************************/
 
int main (int argc, char *argv[]) {

  printf("%s --> start\n",argv[3]);

	if(argc < 5){
		printf("Set average jump length (in units of R), number of jumps and realization counter\n");
		exit(EXIT_SUCCESS);
	}
	
	char *out_name, *createdir, *dir;
	FILE *out_traj;
	
	int seed;
	int Njumps;
	int counter;
	double x0, y0, z0, rho0, eta0;
	double x, y, z;
	float **pts;
	double xnew, ynew, znew, rhonew;
	double alpha, beta, arstar, r;

	double u[5]; /* 5 random numbers uniformly distributed in (0,1) */	
	
	/* Initialization of the randomizer */
	seed = atoi(argv[4]);
	srand(seed);
	rlxd_init(1,rand());

	/*
	 *	Set constants
	 */
	Njumps = atoi(argv[2]);
	m = 2./R/atof(argv[1]);
	lambdaextr[0] = 0.000001*m;
	lambdaextr[1] = 2.41/R; 	/* first zero of J0, plus a bit */
	if (m < lambdaextr[1])
		lambdaextr[1] = m - 1.e-10;
	lambda = Zbisection(funcforlambda, lambdaextr, 1.e-6);
	srml = sqrt(m*m - lambda*lambda);
	return 0;
	pts = malloc((Njumps)*sizeof(long int));
	for(counter=0; counter<Njumps; counter++)
		pts[counter] = malloc(3*sizeof(float));
	
	/* printf("\nm=%lf\t\tlambda = %lf\n\n", m, lambda); */	
	
	/* Set extremes for phi and xi = cos(theta) */
	minPhi = -pi;
	minXi  = -1.;
	
	/*
	 *	Open output files
	 */
	out_name = malloc(100*sizeof(char));
	dir = malloc(100*sizeof(char));
	createdir = malloc(100*sizeof(char));
	sprintf(dir, "output/avjmp_%.3lf", atof(argv[1]));
	sprintf(createdir, "mkdir -p %s", dir);
	sprintf(out_name, "%s/rep_%d.dat", dir, atoi(argv[3]));

	system(createdir);
	out_traj = fopen(out_name,"w");
	
	/*
	 *	Set initial point of the trajectory
	 */
	ranlxd(u,2);
	w = u[0];
	rhoextr[0] = 0.;
	rhoextr[1] = R;

	rho0 = Zbisection(funcforrho, rhoextr, 1.e-6);
	eta0 = 2.*pi*(u[1]-.5);
	x0 = rho0*cos(eta0);
	y0 = rho0*sin(eta0);
	z0 = 0.;

	x=x0; y=y0; z=z0;	
	
	counter=0;
	while(counter<Njumps) {
	
		rho = sqrt(x*x + y*y);
		eta = argument(x,y);

		ranlxd(u, 5);	/* 4 random real numbers ~ U(0,1) */
	
		/*
		 *	Generation of phi and xi with the inversion method (numerically)
		 */
		normPhi = .5/pi/(1. - R*srml*gsl_sf_bessel_I0( rho*srml)*gsl_sf_bessel_K1( R*srml));
		phi	= cdfInversion(pdfPhi, -pi, u[0], 1.e-3, 1.e-6);
		
		b = B(phi);
		normXi = 2./(m*m - lambda*lambda)*(1. - b*srml*gsl_sf_bessel_K1(b*srml));
		
		if(u[4]<.5)
			xi	= tanh(atanh(lambda/m) + cdfInversion(pdfAtanhXi, 0., .5*u[1], 1.e-3, 1.e-8));
		else
			xi	= tanh(atanh(lambda/m) - cdfInversion(pdfAtanhXi, 0., .5*u[1], 1.e-3, 1.e-8));
	
		/*
		 *	Generation of r with inversion method via Lambert W function
		 */
		alpha = m - lambda*xi;
		if (xi >= 1.) {
			xi = 1.; 
			beta = 1.;
		} else if (xi <= -1.) {
			xi = -1.;
			beta = 1.;
		} else {
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
			pts[counter][0] = (float)x;
			pts[counter][1] = (float)y;
			pts[counter][2] = (float)z;
			counter++;
		}
	}
	fwrite(seed, sizeof(int), 1, out_traj);
	for(counter=0; counter<Njumps; counter++)
		fwrite(pts[counter], sizeof(float), 3, out_traj);
	
	fclose(out_traj);
	
	printf("%s --> end\n",argv[3]); fflush(stdout);
	exit(EXIT_SUCCESS);

	
}
