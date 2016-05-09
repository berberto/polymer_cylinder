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
double normRho;			/* normalization constants */
double w;


/*
 *	Function whose first zero is the value of lambda
 */
double funcforlambda (double var){
	double u = sqrt(m*m - var*var);
	return  var*gsl_sf_bessel_J1(var*R) - u * gsl_sf_bessel_J0(var*R) * gsl_sf_bessel_K1_scaled(R*u)/gsl_sf_bessel_K0_scaled(R*u);
}



/*
 *	Stationary distribution for rho
 */
double pdfRho (double r) {
	return r * gsl_sf_bessel_J0(lambda*r) * gsl_sf_bessel_J0(lambda*r) / normRho;
}


/*
 *	Function whose zero is the the rho sampled according its stationary distribution
 */
double funcforrho (double rho) {
	return (.5*rho*rho)*( gsl_sf_bessel_J0(lambda*rho)*gsl_sf_bessel_J0(lambda*rho) + gsl_sf_bessel_J1(lambda*rho)*gsl_sf_bessel_J1(lambda*rho) )/normRho - w;
}



/*******************************************************************************
 *
 *	MAIN ROUTINE
 *
 ******************************************************************************/
 
int main (int argc, char *argv[]) {


/*	if(argc < 3){
		printf("Please set number of parameters and points to sample.\n");
		exit(EXIT_SUCCESS);
	} */
	
	/* Seed of the randomizer */
	int seed;
	FILE *file_seeds;
	file_seeds=fopen("/dev/urandom","rb");
	fread(&seed, sizeof(int), 1, file_seeds);
	seed = abs(seed);


	/*
	 * Declaration of variables
	 */
	FILE *input, *output;
	int i, j, counter, par;
	int Npar; 				/* number of parameters avjmp */
	int Npoints, Npairs;	/* points and pair of points */
	double avjmp;
	double *rho, *eta;
	double *rtrans, rtransav, rtranssd;

	double u[2]; /* 2 random numbers uniformly distributed in (0,1) */	
	
	/* Initialization of the randomizer */
	srand(seed);
	rlxd_init(1,rand());


	/*
	 *	Set parameters
	 */
	output	= fopen("transverse_asy.dat", "w");
	scanf("%d", &Npar);
	Npoints	= atoi(argv[1]);
	Npairs	= (int)(Npoints*(Npoints-1)/2);
	
	rho		= malloc(Npoints*sizeof(double));
	eta		= malloc(Npoints*sizeof(double));
	rtrans	= malloc(Npairs*sizeof(double));

	printf("%d\n", Npar);


	for(par=0; par<Npar; par++){
		
		scanf("%lf", &avjmp);
		printf("%lf\t", avjmp);
		m = 2./R/avjmp;
		lambdaextr[0] = 0.000001*m;
		lambdaextr[1] = 2.41/R; 	/* first zero of J0, plus a bit */
		if (m < lambdaextr[1])
			lambdaextr[1] = m - 1.e-10;
		lambda = Zbisection(funcforlambda, lambdaextr, 1.e-6);
		srml = sqrt(m*m - lambda*lambda);
	
		normRho = (.5*R*R)*(gsl_sf_bessel_J0( lambda*R ) * gsl_sf_bessel_J0(lambda*R) + gsl_sf_bessel_J1(lambda*R) * gsl_sf_bessel_J1(lambda*R));
	
		/*
		 *	Generate points with the stationary distribution in the circle
		 */
		for(counter=0; counter < Npoints; counter++){
			ranlxd(u,2);
			w = u[0];
			rhoextr[0] = 0.;
			rhoextr[1] = R;
			rho[counter] = Zbisection(funcforrho, rhoextr, 1.e-6);
			eta[counter] = 2.*pi*(u[1]-.5);
		}
	
		/*
		 *	Calculate average and std dev of the distance between pairs of points
		 */
		counter = 0;
		rtransav = 0.;
		rtranssd = 0.;
		for (i=0; i<Npoints; i++){
			for(j=0; j<i; j++){
				rtrans[counter] = sqrt(rho[i]*rho[i] + rho[j]*rho[j] - 2.*rho[i]*rho[j]*cos(eta[i]*eta[j]));
				rtransav += rtrans[counter]/Npairs;
				rtranssd += rtrans[counter]*rtrans[counter]/Npairs;
				counter++;
			}
		}
		rtranssd = sqrt(rtranssd - rtransav*rtransav);

		fprintf(output, "%lf\t%lf\t%lf\n", avjmp, rtransav, rtranssd);
		printf("%lf\t%lf\n", rtransav, rtranssd);
	}
	
	fclose(output);
	
	exit(EXIT_SUCCESS);
	
}
