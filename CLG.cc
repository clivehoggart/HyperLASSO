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

#include "CLG.h"

extern "C" {
   extern void pbdv_(double*,double*,double*,double*);
}
double GetF( double d, double eta );
void getPenaltyHLasso( const double&, const double&, const double&, double&, double& );

using namespace std;

CLG::CLG()
{
}

void CLG::SetModel( short imodel, short idisease_models, double ibeta0 )
{
   disease_models=idisease_models;
   model = imodel;
   beta00 = ibeta0;
}

CLG::CLG( short imodel, short idisease_models )
{
   disease_models=idisease_models;
   model = imodel;
}

CLG::CLG( unsigned long in, vector<short> iy, vector<bool> imissing_target, short imodel, short idisease_models )
{
   disease_models=idisease_models;
   n = in;
   y = iy;
   missing_target = imissing_target;
   eta = new double[n];
   error = new double[n];
   logistic = true;
   model = imodel;
}

CLG::CLG( unsigned long in, vector<double> iy, vector<bool> imissing_target, short imodel, short idisease_models )
{
   disease_models=idisease_models;
   n = in;
   yy = iy;
   missing_target = imissing_target;
   eta = new double[n];
   error = new double[n];
   logistic = false;
   model = imodel;
}

void CLG::SetResponse( unsigned long in, vector<short> iy, vector<bool> imissing_target )
{
   n = in;
   y = iy;
   missing_target = imissing_target;
   eta = new double[n];
   error = new double[n];
   logistic = true;
}

void CLG::SetResponse( unsigned long in, vector<double> iy, vector<bool> imissing_target )
{
   n = in;
   yy = iy;
   missing_target = imissing_target;
   eta = new double[n];
   error = new double[n];
   logistic = false;
}

CLG::~CLG()
{
   delete [] eta;
   delete [] error;
}

void CLG::SetPriorParams( vector<double> ilambda, vector<double> igamma )
{
   lambda = ilambda;
   gamma = igamma;
   unsigned int p = lambda.size();
   DerivativeAtOrigin.resize(p);
   DDerivativeAtOrigin.resize(p);

   if( model == 0 )
     for( unsigned int i = 0; i < p; i++ ){
       getPenaltyHLasso( 0.0, lambda[i], gamma[i], DerivativeAtOrigin[i], DDerivativeAtOrigin[i] );
     }
   else if( model == 1 ){
     gamma.assign( p, 0 );
     for( unsigned int i = 0; i < p; i++ ){
       DerivativeAtOrigin[i] = lambda[i];
       DDerivativeAtOrigin[i] = 0.0;
     }
   }else
     gamma.assign( p, 0 );     
}

void CLG::SetPriorParams_d( double ilambda,
			    double igamma, short imodel )
{
   lambda_d = ilambda;
   gamma_d = igamma;
   model_d = imodel;

   if( model_d == 0 )
     getPenaltyHLasso( 0.0, lambda_d, gamma_d, DerivativeAtOrigin_d, DDerivativeAtOrigin_d );
   else{
     DerivativeAtOrigin_d = lambda_d;
     DDerivativeAtOrigin_d = 0.0;
   }
}

unsigned short CLG::Maximize( double &beta0, double* beta,
			      short* beta_type, double* beta_d,
			      vector<unsigned long>& ptr_j,
			      vector<unsigned long>& ptr_d,
			      const vector< vector<double> >& mean,
			      const vector< vector<double> >& sd,
			      const vector<short>& X_work,
			      const vector<bool>& missing,
			      const vector<double>& X_d,
			      const vector<bool>& missing_d, double& LL)
{
  bool change_in_status;
  unsigned short iter = 0;
  unsigned long i,j,jj;
  double sumdeltaeta, sumeta, rss=0;
  unsigned long p_work = ptr_j.size();
  unsigned long p = X_work.size() / n;
  unsigned long p_d = X_d.size() / n;

  tau = 1.0;
  //  double tau_alpha = 2.0, tau_beta = 1.0;

  delta.resize(p);
  for( jj = 0; jj < p; jj++ ){
    delta[jj] = 1.0;
  }
  
  delta_d.resize(p_d);
  for( jj = 0; jj < p_d; jj++ ){
    delta_d[jj] = 1.0;
  }
  
  vector< vector<double> > SumXpos(p), SumXneg(p);
  for( jj = 0; jj < p_work; jj++ ){
    j=ptr_j[jj];
    SumXpos[j].resize(disease_models);
    SumXneg[j].resize(disease_models);
    for( short type=0; type < disease_models; type++ ){
      SumXpos[j][type] = 0.0;
      SumXneg[j][type] = 0.0;
      for( i = 0; i < n; i++ ){
	double xx;
	if( type == 0 )
	  xx = (double)X_work[i*p+j];
	else if( type == 1 ){
	  if( X_work[i*p+j] == 2 )
	    xx = 1;
	  else
	    xx = 0;
	}else if( type == 2 ){
	  if( X_work[i*p+j] == 0 )
	    xx = 0;
	  else
	    xx = 1;
	}else{
	  if( X_work[i*p+j] == 1 )
	    xx = 1;
	  else
	    xx = 0;
	}
	xx = (xx-mean[j][type])/sd[j][type];
	if( logistic && !missing[i*p+j] && !missing_target[i] ){
	  if( y[i]*xx > 0 )
	    SumXpos[j][type] += y[i]*xx;
	  else
	    SumXneg[j][type] += y[i]*xx;
	}
      }
    }
  }
  
  vector<double> SumXpos_d(p_d), SumXneg_d(p_d);
  for( jj = 0; jj < p_d; jj++ ){
    j=ptr_d[jj];
    for( i = 0; i < n; i++ ){
      if( logistic && !missing_d[i*p_d+j] ){
	if( y[i] * X_d[i*p_d+j] > 0 && !missing_target[i] )
	  SumXpos_d[j] += y[i] * X_d[i*p_d+j];
	else
	  SumXneg_d[j] += y[i] * X_d[i*p_d+j];
      }
    }
  }
  

  if( logistic )
    getEta( beta0, beta, beta_type, beta_d, X_work, missing,
	    p, mean, sd, X_d, missing_d, p_d );
  else
    getEta_lin( beta0, beta, beta_type, beta_d, X_work, missing,
		p, mean, sd, X_d, missing_d, p_d );
//   if( !logistic ){
//     rss = 0;
//     for( i = 0; i < n; i++ )
//       if( !missing_target[i] )
// 	rss += yy[i]*yy[i]-eta[i]*yy[i];
//     tau = (0.5*n + tau_alpha-1.0) / (tau_beta + 0.5*rss);
//   }
   do{
      change_in_status = false;
      sumdeltaeta=0.0; sumeta=0.0;
      for( i = 0; i < n; i++ )
	error[i] = 0.0;
      MoveDoubleCovariates( beta_d, X_d, missing_d, p_d,
			    SumXpos_d, SumXneg_d, ptr_d,
			    change_in_status );
      MoveShortCovariates( beta, beta_type, X_work, missing,
			   p_work, p, mean, sd,
			   SumXpos, SumXneg, ptr_j, change_in_status );
      MoveIntercept( beta0 );
      for( i = 0; i < n; i++ ){
	sumeta += fabs(eta[i]);
	sumdeltaeta += fabs(error[i]);
      }
//       if( !logistic ){
// 	rss = 0;
// 	for( i = 0; i < n; i++ )
// 	  if( !missing_target[i] )
// 	    rss += yy[i]*yy[i]-eta[i]*yy[i];
// 	tau = (0.5*n + tau_alpha-1.0) / (tau_beta + 0.5*rss);
// 	cout << "tau: " << tau << endl;
//       }
      iter++;
   }while( sumdeltaeta/(1.0+sumeta) > 0.0005 || change_in_status == true );
   
   if( logistic )
     LL = getLogLikelihood();
   else
     LL = getLogLikelihood_lin();

//    double log2pi = 1.837877;
//    double logdet = log2pi*(double)p_work/2
//      -0.5*LogDeterminantInfo( beta, beta_type, beta_d, X_work,
// 			      missing, X_d, missing_d, mean, sd );
   return iter;
}



void CLG::Move( bool &movement, double &delta_beta, double beta, double delta,
		double penalty, double dpenalty, double df, double ddf,
		double lambda_local, double gamma_local, short model_local )
{
  if( model_local <= 1 ){
    if( beta == 0.0 ){
      delta_beta=(df-penalty)/(ddf+dpenalty);
      if( delta_beta < 0.0 ){
	delta_beta=(df+penalty)/(ddf+dpenalty);
	if( delta_beta > 0.0 ){
	  delta_beta = 0.0;
	}else
	  movement = true;
      }else
	movement = true;
    }else{
      movement = true;
      if( model_local == 0 ){
	getPenaltyHLasso( beta, lambda_local, gamma_local, penalty, dpenalty );
      }else{
	penalty = lambda_local;
	dpenalty = 0.0;
      }
      if( beta < 0.0 ){
	delta_beta=(df+penalty)/(ddf+dpenalty);
	if( delta_beta+beta > 0.0 )
	  delta_beta = -beta;
      }
      else{
	delta_beta=(df-penalty)/(ddf+dpenalty);
	if( delta_beta+beta < 0.0 )
	  delta_beta = -beta;
      }
    }
  }
  else{//ridge regression
    movement = true;
    penalty = beta*lambda_local;
    dpenalty = lambda_local;
    delta_beta = (df-penalty)/(ddf+dpenalty);
  }
  if( isnan(delta_beta) )
    delta_beta = 0.0;
  if( delta_beta < -delta )
    delta_beta = -delta;
  if( delta_beta > delta )
    delta_beta = delta;
}

void getPenaltyHLasso( const double &b, const double &lambda,
		       const double &gamma, double &penalty,
		       double &dpenalty )
{
   double x,v,temp,D1,D2,D3;
   x=fabs(b)/gamma;
   v = -2.0*lambda-1.0;
   pbdv_(&v,&x,&D1,&temp);
   v--;
   pbdv_(&v,&x,&D2,&temp);
   v--;
   pbdv_(&v,&x,&D3,&temp);
   if( D3 == 0 ){
     cout << "Beta = " << b << ", value out of range for specified prior.\n";
     exit(0);
   }
   penalty = 2.0*(lambda+0.5)*D2/(gamma*D1);
   dpenalty = 4.0*((lambda+0.5)*(lambda+1.0)*D3/D1
   		   -(lambda+0.5)*(lambda+0.5)*D2*D2/(D1*D1) )
     / (gamma*gamma);
}

double CLG::LogDeterminantInfo( double* beta, short* beta_type,
				double* beta_d,
				const vector<short>& X_work,
				const vector<bool>& missing,
				const vector<double>& X_d,
				const vector<bool>& missing_d,
				const vector<vector<double> >& mean,
				const vector<vector<double> >& sd )
{
  unsigned long p = X_work.size() / n;
  unsigned long p_d = X_d.size() / n;
  double element, Xi, Xj;
   vector<double> zeta(n);
   vector<short> selected_snps;
   vector<short> selected_vars;

   for( unsigned int i = 0; i < p; i++ )
     if( beta[i] != 0 )
       selected_snps.push_back(i);

   for( unsigned int i = 0; i < p_d; i++ )
     if( beta_d[i] != 0 )
       selected_vars.push_back(i);

   unsigned short p_work = selected_snps.size() + selected_vars.size();  
   gsl_matrix *_matrix = gsl_matrix_alloc(p_work,p_work);

   for( unsigned int i = 0; i < n; i++ ){
      zeta[i] = exp(eta[i]) / (( 1.0 + exp(eta[i]) )*( 1.0 + exp(eta[i]) ));
   }

   for( unsigned int i = 0; i < p_work; i++ ){
     element = 0.0;
     for( unsigned int k = 0; k < n; k++ ){
       if( i < selected_snps.size() )
	 Xi = getXij( beta_type, mean, sd, X_work, missing, k, i );
       else{
	 if( !missing_d[k*p_d+(i-selected_snps.size())] )
	   Xi = X_d[k*p_d+(i-selected_snps.size())];
	 else
	   Xi = 0.0;
       }
       element += Xi * Xi * zeta[k];
     }
     gsl_matrix_set(_matrix, i, i, element+lambda[0]);
     for( unsigned int j = i+1; j < p_work; j++ ){
       element=0.0;
       for( unsigned int k = 0; k < n; k++ ){
	 if( i < selected_snps.size() )
	   Xi = getXij( beta_type, mean, sd, X_work, missing, k, i );
	 else{
	   if( !missing_d[k*p_d+(i-selected_snps.size())] )
	     Xi = X_d[k*p_d+(i-selected_snps.size())];
	   else
	     Xi = 0.0;
	 }
	 if( j < selected_snps.size() )
	   Xj = getXij( beta_type, mean, sd, X_work, missing, k, j );
	 else{
	   if( !missing_d[k*p_d+(j-selected_snps.size())] )
	     Xj = X_d[k*p_d+(j-selected_snps.size())];
	   else
	     Xj = 0.0;
	 }

	 element += Xi * Xj * zeta[k];
       }
       gsl_matrix_set(_matrix, i, j, element);
       gsl_matrix_set(_matrix, j, i, element);
     }
   }
   
   gsl_permutation *perm = gsl_permutation_alloc(p_work);
   int signum;
   gsl_linalg_LU_decomp (_matrix, perm, &signum);
   double det = gsl_linalg_LU_lndet( _matrix );
   gsl_permutation_free(perm);
   gsl_matrix_free(_matrix);
   return det;
}
   
double CLG::getXij( const short* beta_type,
		    const vector< vector<double> >& mean,
		    const vector< vector<double> >& sd,
		    const vector<short>& X_work,
		    const vector<bool>& missing, long i, long j )
{  
  unsigned long p = X_work.size() / n;
  double X;
  if( !missing[i*p+j] ){
    if( beta_type[j] == 0 )
      X = ((double)X_work[i*p+j]-mean[j][0])/sd[j][0];
    else if( beta_type[j] == 1 )
      X = (2*(double)(X_work[i*p+j]/2)-mean[j][1])/sd[j][1];
    else
      X = (2*(double)((2-X_work[i*p+j])/2)-mean[j][2])/sd[j][2];
  }else
    X = 0;
  return X;
}

void CLG::getEta( double beta0, const double *beta,
		  const short *beta_type, const double *beta_d,
		  const vector<short>& X_work,
		  const vector<bool>& missing, unsigned long p,
                  const vector<vector<double> >& mean,
                  const vector<vector<double> >& sd,
		  const vector<double>& X_d,
		  const vector<bool>& missing_d,
		  unsigned long p_d )

{
  EtaMin =  999999.9;
  EtaMax = -999999.9;
   double X;
   for( unsigned int i = 0; i < n; i++ ){
     if( !missing_target[i] ){
       eta[i] = beta0*y[i];
       for( unsigned int j = 0; j < p; j++ ){
         if( !missing[i*p+j] ){
	   X = getXij( beta_type, mean, sd, X_work, missing, i, j );
	   eta[i] += X * beta[j] * y[i];
	   if( isnan(eta[i]) )
	     cout << "here\n";
         }
       }
       for( unsigned int j = 0; j < p_d; j++ ){
         if( !missing_d[i*p_d+j] ){
	   eta[i] += X_d[i*p_d+j] * beta_d[j] * y[i];
	   if( isnan(eta[i]) )
	     cout << "here\n";
         }
       }
       if( eta[i] > EtaMax )
	 EtaMax = eta[i];
       if( eta[i] < EtaMin )
	 EtaMin = eta[i];
     }
   }
}

double CLG::getLogLikelihood()
{
   double l=0,theta;
   for( unsigned int i = 0; i < n; i++ ){
      theta = 1.0 / ( 1.0 + exp(-eta[i]) );
      l += log(theta);
   }
   return l;
}

void CLG::getEta_lin(double beta0, const double *beta,
		     const short *beta_type, const double *beta_d,
		     const vector<short>& X_work,
		     const vector<bool>& missing,
		     unsigned long p,
		     const vector<vector<double> >& mean,
		     const vector<vector<double> >& sd,
		     const vector<double>& X_d,
		     const vector<bool>& missing_d,
		     unsigned long p_d )
{
   double X;
   for( unsigned int i = 0; i < n; i++ ){
      eta[i] = beta0;
      for( unsigned int j = 0; j < p; j++ ){
	if( !missing[i*p+j] ){
	  if( beta_type[j] == 0 )
	    X = ((double)X_work[i*p+j]-mean[j][0])/sd[j][0];
	  else if( beta_type[j] == 1 )
	    X = (2*(double)(X_work[i*p+j]/2)-mean[j][1])/sd[j][1];
	  else
	    X = (2*(double)((2-X_work[i*p+j])/2)-mean[j][2])/sd[j][2];
	  eta[i] += X*beta[j];
	  if( isnan(eta[i]) )
	    cout << "here\n";
	}
      }
      for( unsigned int j = 0; j < p_d; j++ ){
	if( !missing_d[i*p_d+j] ){
	  eta[i] += X_d[i*p_d+j]*beta_d[j];
	  if( isnan(eta[i]) )
	    cout << "here\n";
	}
      }
   }
}

double CLG::getLogLikelihood_lin()
{
   double l=0;
   for( unsigned int i = 0; i < n; i++ ){
     if( !missing_target[i] )
       l += -0.5 * ( yy[i] - eta[i] ) * ( yy[i] - eta[i] );
   }
   return l;
}

double GetF( double d, double eta )
{
  double F;
  if( eta > d )
    F =  1.0 / ( 2.0 + exp(fabs(eta)-d) + exp(d-fabs(eta)) );
  else
    F = 0.25;
  return F;
}

void CLG::MoveDoubleCovariates(double* beta_d,
			       const vector<double>& X_d,
			       const vector<bool>& missing,
			       unsigned long p_d,
			       vector<double>& SumXpos_d,
			       vector<double>& SumXneg_d,
			       std::vector<unsigned long>& ptr_d,
			       bool change_in_status )
{
  double F, df, ddf, change;
  for( unsigned long jj = 0; jj < p_d; jj++ ){
    bool movement = false;
    unsigned long j = ptr_d[jj];
    double delta_beta = 0.0;
    if( !logistic || beta_d[j] != 0 || model_d == 2 ||
	SumXpos_d[j]/(1+exp(EtaMin))+SumXneg_d[j]/(1+exp(EtaMax)) > DerivativeAtOrigin_d ||
	SumXneg_d[j]/(1+exp(EtaMin))+SumXpos_d[j]/(1+exp(EtaMax)) < -DerivativeAtOrigin_d ){
      df=0.0;
      ddf=0.0;
      for( unsigned long i = 0; i < n; i++ ){
	if( !missing_target[i] ){
	  if( isnan(eta[i]) )
	    cout << "here\n";
	  if( !missing[i*p_d+j] ){
	    if( logistic ){
	      F = GetF( fabs(delta_d[j]*X_d[i*p_d+j]), eta[i] );
	      df += X_d[i*p_d+j]*(double)y[i]/(1.0+exp(eta[i]));
	      ddf += X_d[i*p_d+j]*X_d[i*p_d+j]*F;
	    }else{
	      df += ( yy[i]-eta[i] ) * X_d[i*p_d+j] * tau;
	      ddf += X_d[i*p_d+j]*X_d[i*p_d+j] * tau;
	    }
	  }
	}
      }
      if( abs(df) > DerivativeAtOrigin_d || beta_d[j] != 0 || model_d == 2 ){
	Move( movement, delta_beta, beta_d[j], delta_d[j],
	      DerivativeAtOrigin_d, DDerivativeAtOrigin_d,
	      df, ddf, lambda_d, gamma_d, model_d );
	beta_d[j] += delta_beta;
	if( delta_beta != 0 && beta_d[j] == 0 )
	  change_in_status = true;
      }
    }
    if( movement ){
      for( unsigned long i = 0; i < n; i++ ){
	if( !missing_target[i] ){
	  change=0.0;
	  if( !missing[i*p_d+j] ){
	    if( logistic ){
	      change += delta_beta*X_d[i*p_d+j]*y[i];
	    }else{
	      change += delta_beta*X_d[i*p_d+j];
	    }
	    if( isnan(change) )
	      cout << "here\n";
	    eta[i] += change;
	    error[i] += change;
	    if( eta[i] > EtaMax )
	      EtaMax = eta[i];
	    if( eta[i] < EtaMin )
	    EtaMin = eta[i];
	  }
	}
      }
    }
    if( model_d < 2 ){
      delta_d[j] *= 0.5;
      if( fabs(2.0*delta_beta) > delta_d[j] )
	delta_d[j] = fabs(2.0*delta_beta);
    }
  }
}

void CLG::MoveIntercept( double &beta0 )
{
  double df=0.0, ddf=0.0, change, penalty, delta_beta;
  for( unsigned long i = 0; i < n; i++ ){
    if( !missing_target[i] ){
      if( isnan(eta[i]) )
	cout << "here\n";
      if( logistic ){
	ddf += exp(eta[i])/((1.0+exp(eta[i]))*(1.0+exp(eta[i])));
	df += (double)y[i]/(1.0+exp(eta[i]));
      }else{
	df += ( yy[i]-eta[i] ) * tau;
	ddf += tau;
      }
    }
  }
  //  cout << beta0 << " " << beta00 << endl;
  penalty = (beta0-beta00);
  delta_beta = (df-penalty)/(ddf+1.0);
  beta0 += delta_beta;
  //  cout << delta_beta << " " << df << " " << ddf << " " << penalty << endl;
  for( unsigned long i = 0; i < n; i++ ){
    if( !missing_target[i] ){
      change=0.0;
      if( logistic ){
	change += delta_beta*y[i];
      }else{
	change += delta_beta;
      }
      if( isnan(change) )
	cout << "here\n";
      eta[i] += change;
      error[i] += change;
      if( eta[i] > EtaMax )
	EtaMax = eta[i];
      if( eta[i] < EtaMin )
	EtaMin = eta[i];
    }
  }
}

void CLG::MoveShortCovariates(double* beta,  short* beta_type,
			      const vector<short>& X_work,
			      const vector<bool>& missing,
			      unsigned long p_work, unsigned long p,
			      const vector< vector<double> >& mean,
			      const vector< vector<double> >& sd,
			      vector< vector<double> >& SumXpos,
			      vector< vector<double> >& SumXneg,
			      vector<unsigned long>& ptr_j,
			      bool change_in_status )
{
  short type;
  vector<double> delta_beta(disease_models), df(disease_models),
    ddf(disease_models);
  double X, F, change;
  for( unsigned long jj = 0; jj < p_work; jj++ ){
    bool movement = false, test = false;
    unsigned long j=ptr_j[jj];
    if( logistic ){
      for( type=0; type < disease_models; type++ ){
	delta_beta[type] = 0.0;
	if( SumXpos[j][type]/(1+exp(EtaMin))+SumXneg[j][type]/(1+exp(EtaMax)) > DerivativeAtOrigin[j] ||
	    SumXneg[j][type]/(1+exp(EtaMin))+SumXpos[j][type]/(1+exp(EtaMax)) < -DerivativeAtOrigin[j] )
	  test = true;
      }
    }
    if( !logistic || beta[j] != 0 || model == 2 || test ){
      type = -1;
      double beta_current;
      do{
	if( beta[j] != 0 ){
	  type = beta_type[j];
	  beta_current=beta[j];
	}else{
	  type++;
	  beta_current=0.0;
	}
	df[type]=0.0;
	ddf[type]=0.0;
	for( unsigned long i = 0; i < n; i++ ){
	  if( !missing_target[i] ){
	    if( isnan(eta[i]) )
	      cout << "here\n";
	    if( !missing[i*p+j] ){
	      short xx;
	      if( type == 0 )
		xx = X_work[i*p+j];
	      else if( type == 1 ){
		if( X_work[i*p+j] == 2 )
		  xx = 1;
		else
		  xx = 0;
	      }else if( type == 2 ){
		if( X_work[i*p+j] == 0 )
		  xx = 0;
		else
		  xx = 1;
	      }else{
		if( X_work[i*p+j] == 1 )
		  xx = 1;
		else
		  xx = 0;
	      }
	      X = ((double)xx-mean[j][type])/sd[j][type];
	      if( logistic ){
		F = GetF( fabs(delta[j]*X), eta[i] );
		df[type] += X*(double)y[i]/(1.0+exp(eta[i]));
		ddf[type] += X*X*F;
	      }else{
		df[type] += ( yy[i]-eta[i] ) * X *tau;
		ddf[type] += X*X*tau;
	      }
	    }
	  }
	}
      }while( beta[j] == 0 && type < disease_models-1 );
      if( beta[j] != 0 )
	type = beta_type[j];
      else{
	type = 0;
	for( short i=1; i < disease_models; i++ ){
	  if( abs(df[type]) < abs(df[i]) ){
	    type = i;
	  }
	}
      }
      if( abs(df[type]) > DerivativeAtOrigin[j] || beta_current != 0 || model == 2 ){
	Move( movement, delta_beta[type], beta_current, delta[j],
	      DerivativeAtOrigin[j], DDerivativeAtOrigin[j],
	      df[type], ddf[type], lambda[j], gamma[j], model );
	beta[j] += delta_beta[type];

	if( beta[j] == 0 )
	  beta_type[j] = 0;
	
	if( delta_beta[type] != 0 && beta[j] == 0 )
	  change_in_status = true;
      }
    }
    if( movement ){
      if( delta_beta[type] != 0 )
	beta_type[j] = type;
      for( unsigned long i = 0; i < n; i++ ){
	if( !missing_target[i] ){
	  change=0.0;
	  if( !missing[i*p+j] ){
	    if( delta_beta[type] != 0 ){
	      short xx;
	      if( type == 0 )
		xx = X_work[i*p+j];
	      else if( type == 1 ){
		if( X_work[i*p+j] == 2 )
		  xx = 1;
		else
		  xx = 0;
	      }else if( type == 2 ){
		if( X_work[i*p+j] == 0 )
		  xx = 0;
		else
		  xx = 1;
	      }else{
		if( X_work[i*p+j] == 1 )
		  xx = 1;
		else
		  xx = 0;
	      }
	      X = ((double)xx-mean[j][type])/sd[j][type];
	      if( logistic ){
		change += delta_beta[type]*X*y[i];
	      }else{
		change += delta_beta[type]*X;
	      }
	    }
	    if( isnan(change) )
	      cout << "here\n";
	    eta[i] += change;
	    error[i] += change;
	    if( eta[i] > EtaMax )
	      EtaMax = eta[i];
	    if( eta[i] < EtaMin )
	      EtaMin = eta[i];
	  }
	}
      }
    }
    if( model < 2 ){
      delta[j] *= 0.5;
      for( short type1=0; type1<disease_models; type1++ )
	if( fabs(2.0*delta_beta[type1]) > delta[j] )
	  delta[j] = fabs(2.0*delta_beta[type1]);
    }
  }
}
