/*
 *  This file is part of the HyperLasso program.
 *  Copyright (c) Clive J Hoggart   <c.hoggart@imperial.ac.uk>
 *
 *  HyperLasso is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License, version 2, as
 *  published by the Free Software Foundation.
 *
 *  HyperLasso is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */ 
// *-*-C++-*-*
#ifndef CLG_H
#define CLG_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_eigen.h>

#include "rand.h"

class CLG
{
  public:
   CLG();
   CLG(short, short);
   CLG( std::vector<short>&, std::vector<unsigned long>&, short, short );
   CLG( std::vector<double>&, std::vector<unsigned long>&, short, short );
   void SetModel(short, short, double);
   void SetResponse( std::vector<short>&, std::vector<unsigned long>& );
   void SetResponse( std::vector<double>&, std::vector<unsigned long>& );
   unsigned short Maximise( double&, double*, short*, double*,
			    std::vector<unsigned long>&,
			    std::vector<unsigned long>&,
			    const std::vector<std::vector<double> >&,
			    const std::vector<std::vector<double> >&,
			    const std::vector< std::vector<bool> >&, const std::vector<bool>&,
			    const std::vector<double>&,
			    const std::vector<bool>&, double& );
   void SetPriorParams( std::vector<std::vector<double> >&, std::vector< std::vector<double> >&,
			std::vector< std::vector<double> >&, bool );
  void SetPriorParams_d( double, double, short );
  std::vector<double> getCorrection( std::vector<unsigned long>&,
				     const std::vector<short>&,
				     const std::vector<std::vector<double> >&,
				     const std::vector<std::vector<double> >&,
				     const std::vector<bool>& );
  double getCorrection();

private:
  bool BF_implement;
  unsigned short TotalObservations;
  double getXij( short, double, double, short, bool );

  /* double LogDeterminantInfo(double*, short*, double*,  */
  /* 			    const std::vector<short>&, const std::vector<bool>&, */
  /* 			    const std::vector<double>&, const std::vector<bool>&, */
  /* 			    const std::vector<std::vector<double> >&, */
  /* 			    const std::vector<std::vector<double> >&); */

  void getEta( double, const double*, const short*, const double*,
	       const std::vector< std::vector<bool> >&,
	       const std::vector<bool>&, unsigned long,
	       const std::vector<std::vector<double> >&,
	       const std::vector<std::vector<double> >&,
	       const std::vector<double>&, const std::vector<bool>&, unsigned long );

  double getLogLikelihood();

  void getEta_lin( double, const double*, const short*, const double*,
		   const std::vector< std::vector<bool> >&, const std::vector<bool>&, unsigned long,
		   const std::vector<std::vector<double> >&,
		   const std::vector<std::vector<double> >&,
		   const std::vector<double>&, const std::vector<bool>&, unsigned long );
  double getLogLikelihood_lin();

  void Move( bool &movement, double &delta_beta, double beta, double delta,
	     double, double, double df, double ddf,
	     double lambda_local, double gamma_local, short model_local );

  void MoveIntercept( double& beta0 );
  void MoveShortCovariates(double* beta,  short* beta_type,
			   const std::vector< std::vector<bool> >& X_work,
			   const std::vector<bool>& missing,
			   unsigned long,unsigned long,
			   const std::vector< std::vector<double> >& mean,
			   const std::vector< std::vector<double> >& sd,
			   std::vector< std::vector<double> >&,
			   std::vector< std::vector<double> >&,
			   std::vector<unsigned long>& ptr_j, bool );
  
void MoveDoubleCovariates(double* beta, 
			  const std::vector<double>& X_d,
			  const std::vector<bool>& missing,
			  unsigned long,
			  std::vector<double>& SumXpos_d,
			  std::vector<double>& SumXneg_d,
			  std::vector<unsigned long>& ptr_d,  bool );

  unsigned long n;
  double beta00;
  std::vector<unsigned long> observed_target;
  std::vector<short> y;
  std::vector<double> yy;
  std::vector<double> eta;
  std::vector<double> error;
  std::vector<std::vector<double> > lambda;
  std::vector<std::vector<double> > gamma;
  std::vector<std::vector<double> > Info;
  double lambda_d;
  double gamma_d;
  double tau;
  std::vector< double > delta;
  std::vector< double > delta_d;
//   std::vector< double > X_work;
  bool logistic;
  short model;
  short model_d;
  short disease_models;
  std::vector<std::vector<double> > DerivativeAtOrigin;
  std::vector<std::vector<double> > DDerivativeAtOrigin;
  double DerivativeAtOrigin_d;
  double DDerivativeAtOrigin_d;
  double EtaMin;
  double EtaMax;
};

#endif /* !defined CLG_H */
