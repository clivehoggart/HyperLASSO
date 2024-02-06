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

CLG::CLG( vector<short>& iy, vector<unsigned long>& iobserved_target, short imodel, short idisease_models )
{
  disease_models=idisease_models;
  n = iobserved_target.size();
  y = iy;
  TotalObservations = y.size();
  observed_target = iobserved_target;
  eta.resize(TotalObservations);
  error.resize(TotalObservations);
  logistic = true;
  model = imodel;
}

CLG::CLG( vector<double>& iy, vector<unsigned long>& iobserved_target, short imodel, short idisease_models )
{
  disease_models=idisease_models;
  n = iobserved_target.size();
  yy = iy;
  TotalObservations = yy.size();
  observed_target = iobserved_target;
  eta.resize(TotalObservations);
  error.resize(TotalObservations);
  logistic = false;
  model = imodel;
}

void CLG::SetResponse( vector<short>& iy, vector<unsigned long>& iobserved_target )
{
  n = iobserved_target.size();
  y = iy;
  TotalObservations = y.size();
  observed_target = iobserved_target;
  eta.resize(TotalObservations);
  error.resize(TotalObservations);
  logistic = true;
}

void CLG::SetResponse( vector<double>& iy, vector<unsigned long>& iobserved_target )
{
  n = iobserved_target.size();
  yy = iy;
  TotalObservations = yy.size();
  observed_target = iobserved_target;
  eta.resize(TotalObservations);
  error.resize(TotalObservations);
  logistic = false;
}

void CLG::SetPriorParams( vector<vector<double> >& ilambda, vector<vector<double> >& igamma,
			  vector<vector<double> >& iInfo, bool iBF_implement )
{
  BF_implement = iBF_implement;
  lambda = ilambda;
  gamma = igamma;
  Info = iInfo;
  unsigned int p = lambda.size();
  DerivativeAtOrigin.resize(p);
  DDerivativeAtOrigin.resize(p);

   if( model == 0 ){
     for( unsigned int i = 0; i < p; i++ ){
       DerivativeAtOrigin[i].resize(disease_models);
       DDerivativeAtOrigin[i].resize(disease_models);
       for( short type=0; type < disease_models; type++ ){
	 getPenaltyHLasso( 0.0, lambda[i][type], gamma[i][type],
			   DerivativeAtOrigin[i][type], DDerivativeAtOrigin[i][type] );
       }
     }
   }else if( model == 1 ){
     for( unsigned int i = 0; i < p; i++ ){
       DerivativeAtOrigin[i].resize(disease_models);
       DDerivativeAtOrigin[i].resize(disease_models);
       for( short type=0; type < disease_models; type++ ){
	 DerivativeAtOrigin[i][type] = lambda[i][type];
	 DDerivativeAtOrigin[i][type] = 0.0;
       }
     }
   }else{
     for( unsigned int i = 0; i < p; i++ ){
       DerivativeAtOrigin[i].resize(disease_models);
       DDerivativeAtOrigin[i].resize(disease_models);
       for( short type=0; type < disease_models; type++ ){
	 DerivativeAtOrigin[i][type] = 0.0;
	 DDerivativeAtOrigin[i][type] = 0.0;
       }
     }
   }
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

unsigned short CLG::Maximise( double &beta0, double* beta,
			      short* beta_type, double* beta_d,
			      vector<unsigned long>& ptr_j,
			      vector<unsigned long>& ptr_d,
			      const vector< vector<double> >& mean,
			      const vector< vector<double> >& sd,
			      const vector< vector<bool> >& X_work,
			      const vector<bool>& missing,
			      const vector<double>& X_d,
			      const vector<bool>& missing_d, double& LL)
{
  bool change_in_status;
  unsigned short iter = 0;
  unsigned long i,ii,j,jj;
  double sumdeltaeta, sumeta, rss=0;
  unsigned long p_work = ptr_j.size();
  unsigned long p = X_work[0].size() / TotalObservations;
  unsigned long p_d = X_d.size() / TotalObservations;

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
  
  vector< vector<double> > SumXpos, SumXneg;
  SumXpos.resize(p);
  SumXneg.resize(p);
  for( jj = 0; jj < p_work; jj++ ){
    j=ptr_j[jj];
    SumXpos[j].resize(disease_models);
    SumXneg[j].resize(disease_models);
    for( short type=0; type < disease_models; type++ ){
      SumXpos[j][type] = 0.0;
      SumXneg[j][type] = 0.0;
      for( ii = 0; ii < n; ii++ ){
	i = observed_target[ii];
	short x = (short)X_work[0][i*p+j] + (short)X_work[1][i*p+j];
	double xx = getXij( type, mean[j][type], sd[j][type],
			    x, missing[i*p+j] );
	if( logistic ){
	  if( double(y[i])*xx > 0 )
	    SumXpos[j][type] += abs(xx);
	  else
	    SumXneg[j][type] += abs(xx);
	}
      }
      //      cout <<  SumXpos[j][type] << " " << SumXneg[j][type] << endl;
    }
  }
  
  vector<double> SumXpos_d(p_d), SumXneg_d(p_d);
  for( jj = 0; jj < p_d; jj++ ){
    j=ptr_d[jj];
    for( ii = 0; ii < n; ii++ ){
      i = observed_target[ii];
      if( logistic && !missing_d[i*p_d+j] ){
	if( y[i] * X_d[i*p_d+j] > 0 )
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
      for( unsigned int ii = 0; ii < n; ii++ ){
	unsigned int i = observed_target[ii];
	error[i] = 0.0;
      }
      MoveDoubleCovariates( beta_d, X_d, missing_d, p_d,
			    SumXpos_d, SumXneg_d, ptr_d,
			    change_in_status );
      MoveShortCovariates( beta, beta_type, X_work, missing,
			   p_work, p, mean, sd,
			   SumXpos, SumXneg, ptr_j, change_in_status );
      MoveIntercept( beta0 );
      for( unsigned int ii = 0; ii < n; ii++ ){
	unsigned int i = observed_target[ii];
	sumeta += abs(eta[i]);
	sumdeltaeta += abs(error[i]);
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
      //      cout << sumdeltaeta << " " << sumeta <<  " " << change_in_status << " " << sumdeltaeta/(1.0+sumeta) << endl;
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
   x=abs(b)/gamma;
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
   		   -(lambda+0.5)*(lambda+0.5)*D2*D2/(D1*D1) ) / (gamma*gamma);
}
   
double CLG::getXij( short type, double mean, double sd, short X_work, bool missing )
{  
  double xx=0;
  if( !missing ){
    if( type == 0 )
      xx = (double)X_work;
    else if( type == 1 ){
      if( X_work == 2 )
	xx = 2;
      else
	xx = 0;
    }else if( type == 2 ){
      if( X_work == 0 )
	xx = 0;
      else
	xx = 2;
    }else{
      if( X_work == 1 )
	xx = 2;
      else
	xx = 0;
    }
    xx = (xx-mean)/sd;
  }else
    xx = 0;
  return xx;
}

void CLG::getEta( double beta0, const double *beta,
		  const short *beta_type, const double *beta_d,
		  const vector< vector<bool> >& X_work,
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
   for( unsigned long ii = 0; ii < n; ii++ ){
     unsigned long i = observed_target[ii];
     eta[i] = beta0*y[i];
     for( unsigned int j = 0; j < p; j++ ){
       if( !missing[i*p+j] ){
	 short x = (short)X_work[0][i*p+j] + (short)X_work[1][i*p+j];
	 X = getXij( beta_type[j], mean[j][beta_type[j]], sd[j][beta_type[j]], x, missing[i*p+j] );
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

double CLG::getLogLikelihood()
{
   double l=0.0,theta;
   for( unsigned int ii = 0; ii < n; ii++ ){
     unsigned int i = observed_target[ii];
     theta = 1.0 / ( 1.0 + exp(-eta[i]) );
     l += log(theta);
   }
   return l;
}

void CLG::getEta_lin(double beta0, const double *beta,
		     const short *beta_type, const double *beta_d,
		     const vector< vector<bool> >& X_work,
		     const vector<bool>& missing,
		     unsigned long p,
		     const vector<vector<double> >& mean,
		     const vector<vector<double> >& sd,
		     const vector<double>& X_d,
		     const vector<bool>& missing_d,
		     unsigned long p_d )
{
   double X;
   for( unsigned int ii = 0; ii < n; ii++ ){
     unsigned long i = observed_target[ii];
      eta[i] = beta0;
      for( unsigned int j = 0; j < p; j++ ){
	if( !missing[i*p+j] ){
	  short x = (short)X_work[0][i*p+j] + (short)X_work[1][i*p+j];
	  X = getXij( beta_type[j], mean[j][beta_type[j]], sd[j][beta_type[j]], x, missing[i*p+j] );
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
   double l=0.0;
   for( unsigned long ii = 0; ii < n; ii++ ){
     unsigned long i = observed_target[ii];
     l += -0.5 * ( yy[i] - eta[i] ) * ( yy[i] - eta[i] );
     //     cout << yy[i] << " " << eta[i] << " " << l << endl;
   }
   return l;
}

double GetF( double d, double eta )
{
  double F;
  if( eta > d )
    F =  1.0 / ( 2.0 + exp(abs(eta)-d) + exp(d-abs(eta)) );
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
      for( unsigned long ii = 0; ii < n; ii++ ){
	unsigned long i = observed_target[ii];
	if( isnan(eta[i]) )
	  cout << "here\n";
	if( !missing[i*p_d+j] ){
	  if( logistic ){
	    F = GetF( abs(delta_d[j]*X_d[i*p_d+j]), eta[i] );
	    df += X_d[i*p_d+j]*(double)y[i]/(1.0+exp(eta[i]));
	    ddf += X_d[i*p_d+j]*X_d[i*p_d+j]*F;
	  }else{
	    df += ( yy[i]-eta[i] ) * X_d[i*p_d+j] * tau;
	    ddf += X_d[i*p_d+j]*X_d[i*p_d+j] * tau;
	  }
	}
      }
      double d = abs(df);
      if( d > DerivativeAtOrigin_d || beta_d[j] != 0 || model_d == 2 ){
	Move( movement, delta_beta, beta_d[j], delta_d[j],
	      DerivativeAtOrigin_d, DDerivativeAtOrigin_d,
	      df, ddf, lambda_d, gamma_d, model_d );
	beta_d[j] += delta_beta;
	if( delta_beta != 0 && beta_d[j] == 0 )
	  change_in_status = true;
      }
    }
    if( movement ){
      for( unsigned long ii = 0; ii < n; ii++ ){
	unsigned long i = observed_target[ii];
	if( !missing[i*p_d+j] ){
	  if( logistic ){
	    change = delta_beta*X_d[i*p_d+j]*y[i];
	  }else{
	    change = delta_beta*X_d[i*p_d+j];
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
    if( model_d < 2 ){
      delta_d[j] *= 0.5;
      double d = abs(2.0*delta_beta);
      if( d > delta_d[j] )
	delta_d[j] = abs(2.0*delta_beta);
    }
  }
}

void CLG::MoveIntercept( double &beta0 )
{
  double df=0.0, ddf=0.0, change, delta_beta;
  for( unsigned long ii = 0; ii < n; ii++ ){
    unsigned long i = observed_target[ii];
    if( isnan(eta[i]) )
      cout << "here\n";
    if( logistic ){
      ddf += exp(eta[i])/((1.0+exp(eta[i]))*(1.0+exp(eta[i])));
      df += (double)y[i]/(1.0+exp(eta[i]));
    }else{
      df -= eta[i];
    } 
  }

  if( logistic )
    delta_beta = (df)/(ddf);
  else
    delta_beta = df/n - beta0;

  beta0 += delta_beta;
  for( unsigned long ii = 0; ii < n; ii++ ){
     unsigned long i = observed_target[ii];
     if( logistic ){
       change = delta_beta*y[i];
     }else{
       change = delta_beta;
     }
     if( isnan(change) )
       cout << "here\n";
     eta[i] += change;
     //     error[i] += change;
     if( eta[i] > EtaMax )
       EtaMax = eta[i];
     if( eta[i] < EtaMin )
       EtaMin = eta[i];
  }
}

void CLG::MoveShortCovariates(double* beta,  short* beta_type,
			      const vector< vector<bool> >& X_work,
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
	double d1 = abs(SumXpos[j][type]/(1+exp(EtaMin))-SumXneg[j][type]/(1+exp(EtaMax)));
	double d2 = abs(SumXpos[j][type]/(1+exp(EtaMax))-SumXneg[j][type]/(1+exp(EtaMin)));
	// cout << SumXpos[j][type] << " " << SumXneg[j][type] << endl;
	// cout << EtaMin << " " << EtaMax << endl;
	// cout << SumXpos[j][type]/(1+exp(EtaMin)) << " " << SumXneg[j][type]/(1+exp(EtaMax)) << endl;
	// cout<<  SumXpos[j][type]/(1+exp(EtaMax)) << " " << SumXneg[j][type]/(1+exp(EtaMin)) << endl;
	// cout << d1 << " " << d2 << " " << DerivativeAtOrigin[j][type] << endl;
	if( d1  > DerivativeAtOrigin[j][type] || d2 > DerivativeAtOrigin[j][type] ){
	  test = true;
	  //	  cout << j << " " << type << endl;
	}
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
	  beta_current=0.0;
	  // vector<double> probs(4);
	  // for( int i=0; i<4; i++ )
	  //   probs[i] = 0.25;
	  // type = SampleFromDiscrete3( probs );
	  type++;
	}
	df[type]=0.0;
	ddf[type]=0.0;
	for( unsigned long ii = 0; ii < n; ii++ ){
	  unsigned long i = observed_target[ii];
	  if( isnan(eta[i]) )
	    cout << "here\n";
	  if( !missing[i*p+j] ){
	    short x = (short)X_work[0][i*p+j] + (short)X_work[1][i*p+j];
	    X = getXij( type, mean[j][type], sd[j][type],
			    x, missing[i*p+j] );
	    if( logistic ){
	      F = GetF( abs(delta[j]*X), eta[i] );
	      df[type] += X*(double)y[i]/(1.0+exp(eta[i]));
	      ddf[type] += X*X*F;
	    }else{
	      df[type] += ( yy[i]-eta[i] ) * X *tau;
	      ddf[type] += X*X*tau;
	    }
	  }
	}
      }while( beta[j] == 0 && type < disease_models-1 );
      if( beta[j] != 0 )
      	type = beta_type[j];
      else if( disease_models == 1 )
	type = 0;
      else if( disease_models > 1 && BF_implement ){
	double w = .21*.21;
	vector<double> BF(4,0);
      	for( short i=0; i < disease_models; i++ ){
	  double d = abs(df[i]);
	  if( !isnan(DerivativeAtOrigin[j][i]) & !isnan(df[i]) &
	      ( d > DerivativeAtOrigin[j][i]) ){
	    BF[i] = exp( 0.5 * w * Info[j][i] * d * d / (w + Info[j][i]) )
	      / sqrt( 1 + w/Info[j][i] );
	  }
      	}
	type = SampleFromDiscrete3( BF );
      }else{
	type = 0;
	//	cout << df[type] << " " << ddf[type] << endl;
	for( short i=1; i < disease_models; i++ ){
	  //	  cout << df[i] << " " << ddf[i] << endl;
	  if( abs(df[type]) < abs(df[i]) ){
	    type = i;
	  }
	}
      }
      //      cout << type << endl;      

      double d = abs(df[type]);
      if( model==2 || d > DerivativeAtOrigin[j][type] || beta_current != 0 ){
	Move( movement, delta_beta[type], beta_current, delta[j],
	      DerivativeAtOrigin[j][type], DDerivativeAtOrigin[j][type],
	      df[type], ddf[type], lambda[j][type], gamma[j][type], model );
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
      for( unsigned long ii = 0; ii < n; ii++ ){
	unsigned long i = observed_target[ii];
	if( !missing[i*p+j] ){
	  if( delta_beta[type] != 0 ){
	    short x = (short)X_work[0][i*p+j] + (short)X_work[1][i*p+j];
	    X = getXij( beta_type[j], mean[j][beta_type[j]],
			sd[j][beta_type[j]], x, missing[i*p+j] );
	    if( logistic ){
	      change = delta_beta[type]*X*y[i];
	    }else{
	      change = delta_beta[type]*X;
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
    if( model < 2 ){
      delta[j] *= 0.5;
      for( short type1=0; type1<disease_models; type1++ ){
	double d = abs(2.0*delta_beta[type1]);
	if( d > delta[j] )
	  delta[j] = abs(2.0*delta_beta[type1]);
      }
    }
  }
}
