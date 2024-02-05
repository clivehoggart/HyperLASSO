#include "rand.h"

using namespace std;

static gsl_rng *RandomNumberGenerator = gsl_rng_alloc( gsl_rng_taus );

double myrand()
{
   return( gsl_rng_uniform( RandomNumberGenerator ) );
}

double myrandRange( double Min, double Max )
{
   double Range = Max - Min;

   return( Min + Range * gsl_rng_uniform( RandomNumberGenerator ) );
}

vector<unsigned long> sample( unsigned long n, unsigned long r )
{
   vector<unsigned long> ret;
   unsigned long a[r], b[n];
   for(unsigned long i = 0; i < n; i++){
      b[i] = i;
   }
   gsl_ran_choose ( RandomNumberGenerator, a, r, b, n, sizeof (unsigned long));
   for(unsigned long i = 0; i < r; i++){
      ret.push_back( a[i] );
   }
   return ret;
}

void smyrand( long seed )
{
   gsl_rng_set(RandomNumberGenerator, static_cast< unsigned long int >( seed ) );
}

double gengam( double bb, double aa )
{
   return( gsl_ran_gamma( RandomNumberGenerator, aa, 1.0 / bb ) );
}

double genbet( double aa, double bb )
{
   return( gsl_ran_beta( RandomNumberGenerator, aa, bb ) );
}

double gennor( double av, double sd )
{
   return( av + gsl_ran_gaussian( RandomNumberGenerator, sd ) );
}

double gennorLtail( double av, double sd, double a )
{
   return( av - gsl_ran_gaussian_tail( RandomNumberGenerator, -(a-av), sd ) );
}

double gennorRtail( double av, double sd, double a )
{
   return( av + gsl_ran_gaussian_tail( RandomNumberGenerator, a-av, sd ) );
}

double ks()
{
  // See Devroye 1986 for explanation of algorithm
  double U, e0, e0new, e1, e1new, g=0.0, xs=0.0, q, pp, z, w, t_prime, e;
  int ok=0, n, accept=0;
  double t = 0.75, p = 0.37283296, pi=3.14159265358;
  t_prime = pi * pi / (8*t*t);
  double u = gsl_rng_uniform( RandomNumberGenerator );
  if(u<p)
  {
     ok=0;
     while(ok==0)
     {
        while(accept==0)
        {
           e0 = -1.0 * log(gsl_rng_uniform( RandomNumberGenerator ));
           e1 = -1.0 * log(gsl_rng_uniform( RandomNumberGenerator ));
           e0new = e0 / (1.0 - (1.0/(2.0 * t_prime)));
           e1new = 2.0 * e1;
           g = t_prime + e0new;
           if(((e0new*e0new) <= (t_prime*e1new*(g+t_prime))) or (((g/t_prime)-1-log(g/t_prime)) <=e1new)) accept=1;
           else accept=0;
        }
        xs = pi / sqrt(8.0*g);
        w=0;
        z = 1/(2*g);
        pp = exp(-g);
        n=1;
        q=1;
        U = gsl_rng_uniform( RandomNumberGenerator );
        while(U>=w && ok==0)
        {
           w += (z*q);
           if(U >=w) {ok=1;}
           n += 2;
           q = pow(pp,n*n-1);
           w -= n*n*q;
            }
     }
     return xs;
  }
  else
  {
     ok=0;
     while(ok==0)
     {
        e = -log(gsl_rng_uniform( RandomNumberGenerator ));
        u = gsl_rng_uniform( RandomNumberGenerator );
        xs = sqrt(t*t+ e/2);
        w=0;
        n=1;
        z=exp(-2*xs*xs);
        while (u>w && ok==0)
        {
           n++;
           w  += n*n*pow(z,n*n-1);
           if(u>=w) ok=1;
           n++;
           w -= n*n*pow(z,n*n-1);
        }
     }
     return xs;
  }
}

void genMultNormal( gsl_matrix* mean, gsl_matrix* A, double *beta )
{
  // A - precision
   unsigned int d = mean->size1;
   gsl_matrix *V;
   gsl_vector *S;
   gsl_vector *work;
   V = gsl_matrix_calloc(d,d);
   S = gsl_vector_calloc(d);
   work = gsl_vector_calloc(d);
   gsl_linalg_SV_decomp( A, V, S, work );

   double draw[d];
   for ( unsigned int i = 0; i < d; i++ )
      draw[i] = gennor( (double)0.0, (double)1.0 );

   for ( unsigned int i = 0; i < d; i++ ){
      beta[i] = 0.0;
      for( unsigned int j = 0; j < d; j++)
         beta[i] += gsl_matrix_get(A,i,j)*draw[j]/sqrt(gsl_vector_get(S,j));
      beta[i] += gsl_matrix_get(mean,i,0);
   }
   gsl_matrix_free( V );
   gsl_vector_free( S );
   gsl_vector_free( work );
}

int genbinomial( int n, double p )
{
   return( gsl_ran_binomial( RandomNumberGenerator, p, n ) );
}

unsigned int genpoi( double mu )
{
   return( gsl_ran_poisson( RandomNumberGenerator, mu ) );
}

vector<int> genmultinomial(int N, vector<double> theta)
{
  int K = theta.size();
  unsigned int n[ K ];
  double p[ K ];
  vector<int> sample( K );
  for(int i = 0; i < K; i++){
      p[i] = theta[i];
  }
  gsl_ran_multinomial(RandomNumberGenerator, K, N, p, n);
  for(int i = 0; i<K; i++){
    sample[i] = n[i];
  }
  return( sample );
}

double MultinomialLikelihood( vector<int> r, vector<int> theta )
{
   if( r.size() != theta.size() ){
      cout << "Length of vector arguments to MultinomialLikelihood not equal." << endl;
      exit(0);
   }
   int K = r.size();
   unsigned int n[ K ];
   double p[ K ];
   for( int i = 0; i < K; i++ ){
      p[i] = theta[i];
      n[i] = r[i];
   }
   return( gsl_ran_multinomial_pdf( K, p , n ) );
}

long ignpoi( double mu )
{
   return( gsl_ran_poisson( RandomNumberGenerator, mu ) );
}

int SampleFromDiscrete( vector<double> probs )
{
   unsigned long numberofelements=probs.size();
   vector<double> cdf;
   cdf.push_back( probs[0] );
   for( unsigned long i = 1; i < numberofelements; i++ )
     cdf.push_back( cdf[i-1] + probs[i] );
   double u = myrand();
   for( unsigned long i = 0; i < numberofelements; i++ )
      cdf[i] /= cdf[ numberofelements - 1 ];
   int k = 0;
   while( u > cdf[k] )
      k++;

   return k;
}

int SampleFromDiscrete( vector<double> *cdf )
{
   int numberofelements = (*cdf).size();
   double u = myrand();
   for( int i=0; i < numberofelements; i++ )
     (*cdf)[i] /= (*cdf)[ numberofelements - 1 ];
   int k = 0;
   while( u > (*cdf)[k] )
      k++;

   return(k);
}

int SampleFromDiscrete2( vector<double> *probs )
{
   int numberofelements = (*probs).size();
   vector<double> cdf( numberofelements );
   cdf[0] = (*probs)[0];
   for( int i = 1; i < numberofelements; i++ )
      cdf[i] = cdf[i-1] + (*probs)[i];
   double u = myrand();
   for( int i = 0; i < numberofelements; i++ )
     cdf[i] /= cdf[ numberofelements - 1 ];
   int k = 0;
   while( u > cdf[k] )
      k++;

   return(k);
}

int SampleFromDiscrete3( double* probs , int numberofelements )
{
   vector<double> cdf( numberofelements );
   cdf[0] = probs[0];
   for( int i = 1; i < numberofelements; i++ )
     cdf[i] = cdf[i-1] + probs[i];
   double u = myrand();
   for( int i = 0; i < numberofelements; i++ )
     cdf[i] /= cdf[ numberofelements - 1 ];
   int k = 0;
   while( u > cdf[k] )
      k++;

   return(k);
}

int SampleFromDiscrete3( vector<double> probs )
{
   unsigned long numberofelements=probs.size();
   vector<double> cdf( numberofelements );
   cdf[0] = probs[0];
   for( unsigned long i = 1; i < numberofelements; i++ )
     cdf[i] = cdf[i-1] + probs[i];
   double u = myrand();
   for( unsigned int i = 0; i < numberofelements; i++ )
     cdf[i] /= cdf[ numberofelements - 1 ];
   int k = 0;
   while( u > cdf[k] )
      k++;

   return(k);
}

int SampleFromDiscrete4( double max, vector<double> LogProbs )
{
   double *probs = new double[LogProbs.size()];
   for( unsigned short i = 0; i < LogProbs.size(); i++ )
      probs[i] = exp(LogProbs[i]-max);

   int k = SampleFromDiscrete3( probs, LogProbs.size() );
   delete [] probs;
   return k;
}

unsigned long SampleFromRange(unsigned long MIN, unsigned long MAX)
{	
   return (unsigned long)( myrand()*(MAX-MIN) + MIN );
}

vector<double> gendirichlet( vector<double> alpha )
{
   int d = alpha.size();
   double sum = 0;
   vector<double> theta( d );

   for( int i = 0; i < d; i++ )
   {
//      assert( (double)alpha(i) > 0 );
      if( alpha[i] > 0 )
         do{
            theta[i] = gengam( 1.0, (double)alpha[i] );
         }while( theta[i] == 0 );
      else
         theta[i] = 0.0;
      sum += theta[i]; 
   }

   assert( sum > 0.0 );
   for( int i = 0; i < d; i++ )
     theta[i] /= sum;

   return( theta );
}

void ddigam(  double *X, double *ddgam  )
{
//  FOLLOWS CLOSELY ALG AS 103 APPL.STATS.(1976) VOL 25
//  CALCS DIGAMMA(X)=D(LOG(GAMMA(X)))/DX
//
//  SET CONSTANTS.SN=NTH STIRLING COEFFICIENT,D1=DIGAMMA(1.)
//
   double S, C, S3, S4, S5, D1, Y, R;

   S = 1.0e-5;
   C = 8.5e0;
   S3 = 8.333333333e-2; 
   S4 = 8.333333333e-3;
   S5 = 3.968253968e-3;
   D1 = -0.5772156649;

//      DATA S,C,S3,S4,S5,D1/1.0D-5,8.5D0,8.333333333D-2,
//    1  8.333333333D-3,3.968253968D-3,-0.5772156649D0/


//  CHECK ARGUMENT IS POSITIVE

   *ddgam=0.0;
   Y=*X;
   if(Y < 0.0){
      std::cout << " ARG OF DIGAMMA FN = " << Y << ". ARG MUST BE POSITIVE " << std::endl;
   }
   assert(Y>=0.0);

//  USE APPROXIMATION IF ARGUMENT .LE.S

   if(Y > S){

//  REDUCE TO DIGAMMA(X+N),(X+N).GE.C

      while( Y < C ){
         *ddgam=*ddgam-(1.0/Y);
         Y=Y+1.0;}
   
//  USE STIRLING IF ARGUMENT .GE.C

      R=1.0/Y;
      *ddgam=*ddgam+log(Y)-0.5*R;
      R=R*R;
      *ddgam=*ddgam-R*(S3-R*(S4-R*S5));}
   else
      *ddgam=D1-1.0/Y;
}

void trigam( double *x, double *trgam )
{
/*
 * closely follows alg. as 121 appl.stats. (1978) 
 * vol 27, 97-99. (b.e. schneider)
 *
 * calculates trigamma(x)=d**2(log(gamma(x)))/dx**2
 */
   double a=1.0e-4,b=5.0,one=1.0,half=0.5,y,z,trigam1=0.0;
   double b2=0.1666666667,b4=-0.03333333333;
   double b6=0.02380952381,b8=-0.0333333333;
/*
 *  b2,b4,b6,b8 are bernoulli numbers
 *
 *  check that argument is positive
 *
 */ 
   z=*x;
/*
 *  use small value approximation if x.le.a
 */
   if(z<=a){ 
      trigam1=1.0/(z*z);
      *trgam=trigam1;}
   else{
/*
 *  increase argument to (x+i).ge.b
 */
      while(z<b){
         trigam1=trigam1+1.0/(z*z);
         z=z+1.0;
      }
/*
 *  apply asymptotic formula if argument.ge.b
 */
      y=1.0/(z*z);
      trigam1=trigam1+half*y+(one+y*(b2+y*(b4+y*(b6+y*b8))))/z;
      *trgam=trigam1;
   }
}
