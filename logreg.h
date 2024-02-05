// *-*-C++-*-*
#ifndef LogReg_H
#define LogReg_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_eigen.h>

class LogReg
{
  public:
   LogReg();
   LogReg( unsigned long, std::vector<short> );
   ~LogReg();
   void SetResponse( unsigned long, std::vector< short > );
   void SetData( unsigned long, const std::vector< double >&, const std::vector<bool>& missing );
   void SetPriorParams( double );
   double MargLikelihoodApprox( const double* );
   double getLogLikelihood( const double* );
   
  private:
   void getEta( const double* );
   double LogDeterminantInfo( const double* );
   double LogDeterminant( gsl_matrix* _matrix );
   unsigned long p_work;
   unsigned long n;
   std::vector<short> y;
   double *eta;
   double lambda;
   std::vector< double > delta;
   std::vector< double > X_work;
   std::vector< bool > missing_work;
};

#endif /* !defined LogReg_H */
