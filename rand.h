#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_linalg.h>
}

double myrand();
///
std::vector<unsigned long> sample( unsigned long n, unsigned long r );
///
double myrandRange( double Min, double Max );
///
void smyrand( long seed );
///
double gengam(double aa,double bb);
///
double genbet(double aa,double bb);
///
double gennor(double av,double sd);
//
double gennorLtail( double av, double sd, double a );
//
double gennorRtail( double av, double sd, double a );
//
void genMultNormal( gsl_matrix, gsl_matrix, double* );
///
double ks();
//
int genbinomial( int n, double p );
//
unsigned int genpoi( double );
//
std::vector<int> genmultinomial( int n, std::vector<double> p );
//
double MultinomialLikelihood( std::vector<int> r, std::vector<double> theta );
///
///Poisson generator
long ignpoi(double mu);
///
int SampleFromDiscrete( std::vector<double> *cdf );
///
int SampleFromDiscrete2( std::vector<double> *probs );
///
int SampleFromDiscrete3( double*, int numberofelements);
///
int SampleFromDiscrete3( std::vector<double>);
///
int SampleFromDiscrete4( double max, std::vector<double> LogProbs );
//
unsigned long SampleFromRange( unsigned long, unsigned long );
//
std::vector<double> gendirichlet( std::vector<double> );
///
void ddigam( double *, double * );
///
void trigam( double *, double * );
///
