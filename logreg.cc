#include "logreg.h"

using namespace std;

LogReg::LogReg()
{
}

LogReg::LogReg( unsigned long in, vector<short> iy )
{
   n = in;
   y = iy;
   eta = new double[n];
}

LogReg::~LogReg()
{
   delete [] eta;
}

void LogReg::SetData( unsigned long p, const vector< double >& X, const vector<bool>& missing )
{
   if( p != p_work ){
      p_work = p;
      delta.resize(p);
   }
   X_work = X;
   missing_work = missing;
}

void LogReg::SetResponse( unsigned long in, vector<short> iy )
{
   n = in;
   y = iy;
   eta = new double[n];
}

void LogReg::SetPriorParams( double ilambda )
{
   lambda = ilambda;
}

double LogReg::MargLikelihoodApprox( const double *beta )
{
   double LL=0;
   double *zeta = new double[n];
   for( unsigned int j = 1; j < p_work; j++ )
      LL -= beta[j]*beta[j];
   LL *= 0.5*lambda;
   LL += 0.5*(p_work-1)*log(lambda);
   cout << LL << " " << getLogLikelihood( beta ) << " ";
   LL += getLogLikelihood( beta );
   for( unsigned int i = 0; i < n; i++ ){
      zeta[i] = exp(eta[i]) / (( 1.0 + exp(eta[i]) )*( 1.0 + exp(eta[i]) ));
   }
   if( p_work > 1 ){
      double LogDet = LogDeterminantInfo( zeta );
      cout << LogDet << " ";
      LL -= 0.5*LogDet;
   }
   delete [] zeta;
   cout << LL << endl;
   return LL;
}

double LogReg::LogDeterminantInfo( const double *zeta )
{
   double element,det=0;
   gsl_matrix *_matrix = gsl_matrix_alloc(p_work-1,p_work-1);
   for( unsigned int i = 1; i < p_work; i++ ){
      element = 0.0;
      for( unsigned int k = 0; k < n; k++ ){
         element += X_work[k*p_work+i] * X_work[k*p_work+i] * zeta[k];
      }
      gsl_matrix_set(_matrix, i-1, i-1, element+lambda);
      for( unsigned int j = i+1; j < p_work; j++ ){
         element=0.0;
         for( unsigned int k = 0; k < n; k++ ){
            element += X_work[k*p_work+i] * X_work[k*p_work+j] * zeta[k];
         }
         gsl_matrix_set(_matrix, i-1, j-1, element);
         gsl_matrix_set(_matrix, j-1, i-1, element);
      }
   }
   det = LogDeterminant( _matrix );
//    for( unsigned int i = 0; i < p_work-1; i++ ){
//        for( unsigned int j = 0; j < p_work-1; j++ ){
//           cout << gsl_matrix_get(_matrix,i,j) << " ";
//        }
//        cout << endl;
//     }
//    cout << det << endl << endl;
   gsl_matrix_free(_matrix);
   return det;
}

double LogReg::LogDeterminant( gsl_matrix* _matrix )
{
   gsl_permutation *perm = gsl_permutation_alloc(p_work-1);
   int signum;
   gsl_linalg_LU_decomp (_matrix, perm, &signum);
   double det = gsl_linalg_LU_lndet( _matrix );
   gsl_permutation_free(perm);
   return det;
}

void LogReg::getEta( const double *beta )
{
   for( unsigned int i = 0; i < n; i++ ){
      eta[i] = beta[0]*y[i];
      for( unsigned int j = 1; j < p_work; j++ )
         if( !missing_work[i*p_work+j] )
            eta[i] += X_work[i*p_work+j]*beta[j]*y[i];
   }
}

double LogReg::getLogLikelihood( const double *beta )
{
   double l=0,theta;
   getEta( beta );
   for( unsigned int i = 0; i < n; i++ ){
      theta = 1.0 / ( 1.0 + exp(-eta[i]) );
      l += log(theta);
   }
   return l;
}
