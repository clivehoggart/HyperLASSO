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
#ifndef HLASSO_H
#define HLASSO_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <map>
#include <string>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cdf.h>

#include "rand.h"
#include "InputData.h"
#include "StringConvertor.h"
#include "CLG.h"
#include "logreg.h"

class CLG;

struct equal
{
   bool operator()( std::vector<unsigned long> pos1, std::vector<unsigned long> pos2 ) const
      {
         bool match = false;
         if( pos1.size() == pos2.size() && pos1.size() != 0 ){
            sort( pos1.begin(), pos1.end() );
            sort( pos2.begin(), pos2.end() );
            unsigned long j=0;
            do{
               if( pos1[j] < pos2[j] )
                  match = true;
               j++;
            }while( j < pos1.size() && pos1[j-1] == pos2[j-1] );
         }else if( pos1.size() < pos2.size() )
            match = true;
         return true;
	 //         return match;
      }
};

typedef std::map< std::vector<unsigned long>, double, equal > model_list;

class HLasso
{
public:
  HLasso( bool, char*, char*, char*, char* );
  void SetParams( unsigned short, double, double, double,
		  unsigned short, double, double, double,
		  bool, bool, short, double, double, double, char*, char* );
  void permute_outcome();
  void getLogPosterior();
  ~HLasso();
  void runCLG( unsigned short, bool, bool, bool );
  void OpenOutFiles(char*);
  void OpenSelectionProbFile(char*);
  void CloseOutFiles(unsigned int);
  void readFile2( const char*, const char*, bool, bool );
  void SetBF( const char *fname );
  void SetGW_BF( double BF, const char *fname );
  void SampleRows();
  void WriteSelectionProbs(unsigned short);

private:
   std::vector<unsigned long> SelectedRows;
   unsigned short prior;
   unsigned short prior_d;
   unsigned long n;
   unsigned long p;
   unsigned long p2;
   double precision;
   std::vector<std::string> Labels;
   std::vector<std::string> Labels_d;
   double beta0;
   std::vector< std::vector<bool> > X;
   std::vector<double> X_d;
   std::vector< std::vector<double> > mean;
   std::vector< std::vector<double> > info;
   std::vector< double > mean_d;
   std::vector< std::vector<double> > sd;
   std::vector< double > sd_d;
   std::vector<bool> missing;
   std::vector<bool> missing_d;
   std::vector<bool> missing_target;
   std::vector<unsigned long> observed_target;
   std::vector<short> y;
   std::vector<double> yy;
   std::vector<unsigned long> case_ptr;
   std::vector<unsigned long> control_ptr;
   double epsilon;
   double epsilon0;
   double epsilon1;
   std::vector< std::vector<double> > lambda;
   std::vector< std::vector<double> > gamma;
   double lambda_d;
   double gamma_d;
   double penalty_d;
   std::ofstream outfile;
   std::ofstream selectionprobfile;
   std::ofstream outfile2;
   bool logistic;
   short disease_models;
   model_list ModelList;
   unsigned long cv_sets;
   unsigned long cv_size;
   unsigned long NumberOfOutputModels;
   double *beta_shrinkage;
   double *beta_d;
   short *beta_type;
   CLG mle_fit;
   CLG shrinkage_fit;
   double data_mean;
   double data_sd;
   std::vector<double> BF;
   std::vector<unsigned short> SumSelected;
   std::vector<unsigned short> SumSelected_d;
};

#endif /* !defined HLasso_H */
