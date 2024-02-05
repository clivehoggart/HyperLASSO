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

#include "HLasso.h"

extern "C" {
   extern void pbdv_(double*,double*,double*,double*);
}

double LogDeterminant(  unsigned int, gsl_matrix* );

#define PR(x) cerr << #x << " = " << x << endl;
using namespace std;
float chisq( vector<short>& y, vector<short>& g, vector<bool> missing,
	     float, float V );
double getNEGLogPrior( double*, unsigned long, double, double );
double getNEGLogPrior( double*, unsigned long, vector<double>, vector<double> );
double getLaplaceLogPrior( double*, unsigned long, double );
double getLaplaceLogPrior( double*, unsigned long, vector<double> );
double getGausLogPrior( double*, unsigned long, double );
double getGausLogPrior( double*, unsigned long, vector<double> );

HLasso::HLasso( bool ilogistic,
                char* infilename, char* targetfilename, char* covfilename,
                unsigned short iprior, double ipenalty, double shape,
		double scale, unsigned short iprior_d, double ipenalty_d,
		double shape_d, double scale_d, bool standardize,
		bool standardize_d, short idisease_models, double alpha,
		double alpha_d, double iBF, char* GW_BFfile, char* BFfile )
{
  logistic = ilogistic;
  NumberOfOutputModels=0;
  cv_sets=100;
  cv_size=300;
  disease_models = idisease_models;
  p2 = 0;
  vector<double> penalty;
  Matrix_s target;

   if( *targetfilename ){
     readFile(targetfilename, target);
     cout << "Dimensions of target file: "
	  << target.size()  << " " << target[0].size() << endl;
     n = target.size();;
     y.resize(n);
     yy.resize(n);
     for( unsigned int i = 0; i < n; i++ ){
       if( target[0].size() == 1 ){
	 if( logistic ){
	   if( target[i][0] == "0" || target[i][0] == "1" ){
	     y[i] = 2*StringConvertor::toInt(target[i][0])-1;
	     missing_target.push_back(false);
	   }else if( target[i][0] == "NA" ){
	     y[i] = 0;
	     missing_target.push_back(true);
	   }else{
	     cout << "Controls coded as 0, cases coded as 1, do not recognize "
		  << target[i][0] << endl;
	     exit(0);
	   }
	 }else{
	   if( target[i][0] != "NA" ){
		yy[i] = StringConvertor::toFloat(target[i][0]);
		y[i] = 1;
		missing_target.push_back(false);
	   }else{
	     yy[i] = 0;
	     y[i] = 0;
	     missing_target.push_back(true);
	   }
	 }
       }else{
	 if( logistic ){
	   if( target[0][i] == "0" || target[0][i] == "1" ){
	     y[i] = 2*StringConvertor::toInt(target[0][i])-1;
	     missing_target.push_back(false);
	   }else if( target[0][i] == "NA" ){
	     y[i] = 0;
	     missing_target.push_back(true);
	   }else{
	     cout << "Controls coded as 0, cases coded as 1, do not recognize "
		  << target[0][i] << endl;
	     exit(0);
	   }
	 }else{
	   if( target[0][i] != "NA" ){
	     yy[i] = StringConvertor::toFloat(target[0][i]);
	     y[i] = 1;
	     missing_target.push_back(false);
	   }else{
	     yy[i] = 0;
	     y[i] = 0;
	     missing_target.push_back(true);
	   }
	 }
       }
     }
   }

   readFile2(infilename, covfilename, !*targetfilename );
   if( p != 0 ){
     if( X.size() / p != n ){
       cout << "Dimesnions of target and genotypes file do not match.\n";
       cout << "Input file has " << X.size()/p << " rows." << endl;
       cout << "Target file has " << n << " rows." << endl;
       exit(0);
     }
   }
   if( p2 != 0 ){
     if( X_d.size() / p2 != n ){
       cout << "Dimesnions of target and phenotypes file do not match.\n";
	 cout << "Input file has " << X_d.size()/p2 << " rows." << endl;
	 cout << "Target file has " << n << " rows." << endl;
	 exit(0);
       }
   }
   cout << n << " "<< p <<  " " << p2 << endl;
   

   mean.resize(p);
   sd.resize(p);
   double s2;
   unsigned short xx;
   unsigned long obs=0;
   for( unsigned int j = 0; j < p; j++ ){
      mean[j].resize(disease_models);
      sd[j].resize(disease_models);
      for( short type = 0; type < disease_models; type++ ){
         s2 = 0.0;
         obs = 0;
         for( unsigned int i = 0; i < n; i++ ){
            if( !missing[i*p+j] ){
               obs++;
               if( type == 0 )
                  xx = X[i*p+j];
	       else if( type == 1 ){
		 if( X[i*p+j] == 2 )
		   xx = 1;
		 else
		   xx = 0;
	       }else if( type == 2 ){
		 if( X[i*p+j] == 0 )
		   xx = 0;
		 else
		   xx = 1;
	       }else{
		 if( X[i*p+j] == 1 )
		   xx = 1;
		 else
		   xx = 0;
	       }
               mean[j][type] += xx;
               s2 += xx*xx;
            }
         }
         mean[j][type] /= obs;
	 s2 = sqrt(s2/obs - mean[j][type]*mean[j][type]);
	 if( s2 == 0 ){
	   cout << "Warning, variable " << Labels[j] << " at position " << j  << " is monomorphic for effect type " << type << endl;
	   sd[j][type] = 1;
	 }else{
	   sd[j][type] = s2;
	 }
      }
      if( !standardize ){
	s2 = sd[j][0];
	for( short type = 0; type < disease_models; type++ ){
	  sd[j][type] /= s2;
	}
      }
   }
   mean_d.resize(p2);
   sd_d.resize(p2);
   for( unsigned int j = 0; j < p2; j++ ){
     mean_d[j] = 0.0;
     sd_d[j] = 0.0;
     for( unsigned int i = 0; i < n; i++ ){
       mean_d[j] += X_d[i*p2+j];
       sd_d[j] += X_d[i*p2+j]*X_d[i*p2+j];
     }
     mean_d[j] /= n;
     sd_d[j] = sqrt(sd_d[j]/n - mean_d[j]*mean_d[j]);
     if( sd_d[j] == 0 ){
       cout << "Warning, variable " << Labels_d[j] << " at position " << j  << " is monomorphic, therefore cannot standardize." << endl;
       sd_d[j] = 1;
     }
     if( !standardize_d )
       sd_d[j] = 1.0;
     for( unsigned int i = 0; i < n; i++ ){
       X_d[i*p2+j] = (X_d[i*p2+j]-mean_d[j])/sd_d[j];
     }
   }

   beta_shrinkage = new double[p];
   beta_type = new short[p];
   double sum=0,sum2=0;
   obs = 0;
   if( logistic ){
      data_mean = 0;
      data_sd = 1;
      for( unsigned int i = 0; i < n; i++ ){
         sum += y[i]+1;
         if( y[i]==1 )
            case_ptr.push_back(i);
         else
            control_ptr.push_back(i);
      }
      case_ptr.push_back(9999);
      control_ptr.push_back(9999);
      beta0 = log((double)sum/(2*n-sum));
      if( alpha != 0 ){
 	penalty.assign( p, sqrt( sum*(2*n-sum)/(4*n) ) * gsl_cdf_ugaussian_Pinv( 1-alpha/2) );
	standardize = true;
      }
   }else{
     for( unsigned int i = 0; i < n; i++ ){
       if( !missing_target[i] ){
	 obs++;
	 sum += yy[i];
	 sum2 += yy[i]*yy[i];
       }
     }
     data_mean = sum/obs;
     data_sd = sqrt(sum2/obs - data_mean*data_mean);
     cout << data_mean << " " << obs << endl;
     cout << "sd: " << data_sd << endl;
     for( unsigned int i = 0; i < n; i++ ){
       yy[i] = (yy[i]-data_mean)/data_sd;
     }
     beta0 = 0;
     if( alpha != 0 )
       penalty.assign( p, sqrt(n) * gsl_cdf_ugaussian_Pinv( 1-alpha/2) );
   }
   if( iBF != 0 ){
     if( *GW_BFfile )
       SetGW_BF( iBF, GW_BFfile );
     else
       BF.assign( p, iBF );
     if( *BFfile )
       SetBF( BFfile );

     // w -- parameter used to fix Bayes Factor, notation as in Wakefield
     double w = .21*.21;
     for( unsigned int j = 0; j < p; j++ ){
       double v, n0, n1;
       if( logistic ){
	 n0 = (double)case_ptr.size()-1.0;
	 n1 = (double)control_ptr.size()-1.0;	  
	 v = (n0+n1)/( n0 * n1 * mean[j][0] * (1.0-mean[j][0]/2) );
       }else
	 v = 1 / (n* .5*mean[j][0] * (1.0-mean[j][0]/2));
       penalty.push_back( sqrt(2.0 * log(BF[j] * sqrt(1.0+w/v)) * (1.0 + v/w) / v) );
       cout << Labels[j] << " " << BF[j] << " " << 1.0-mean[j][0]/2
	    << " " << penalty[j] << endl;
     }
   }

   cout << "Intercept: " << beta0 << endl;
   for( unsigned int j = 0; j < p; j++ ){
      beta_type[j] = 0;
      beta_shrinkage[j] = 0;
   }
   beta_d = new double[p2];
   for( unsigned int j = 0; j < p2; j++ ){
      beta_d[j] = 0;
   }

   prior = iprior;
   prior_d = iprior_d;
   if( penalty.size() == 0 )
     penalty.assign( p, ipenalty );
   //   cout << penalty[0] << " " << penalty[10] << endl;
   penalty_d = ipenalty_d;

   if( prior==0 ){
     lambda.assign( p , shape );
     double x=0, v, temp, D1, D2;
     if( scale != 0 ){
       v = -2.0*lambda[0]-1.0;
       pbdv_(&v,&x,&D1,&temp);
       v--;
       pbdv_(&v,&x,&D2,&temp);
       gamma.assign( p, scale );
       penalty[0] = (2.0*lambda[0]+1)*D2 / (gamma[0]*D1);
     }
     else{
       for( unsigned int j = 0; j < p; j++ ){
	 v = -2.0*lambda[j]-1.0;
	 pbdv_(&v,&x,&D1,&temp);
	 v--;
	 pbdv_(&v,&x,&D2,&temp);
	 gamma.push_back( (2.0*lambda[j]+1)*D2 / (penalty[j]*D1) );
       }
     }
   }else{
     lambda.assign( p , ipenalty );
   }
   
  if( *covfilename ){
     if( prior_d==0 ){
       lambda_d = shape_d;
       double x=0, v, temp, D1, D2;
       v = -2.0*lambda_d-1.0;
       pbdv_(&v,&x,&D1,&temp);
       v--;
       pbdv_(&v,&x,&D2,&temp);
       if( alpha_d != 0 ){
	 if( !logistic )
	   penalty_d = sqrt(n) * gsl_cdf_ugaussian_Pinv( 1-alpha_d/2);
	 else
	   penalty_d = sqrt( sum*(2*n-sum)/(4*n) ) * gsl_cdf_ugaussian_Pinv( 1-alpha_d/2);
       }
       if( scale_d != 0 ){
         gamma_d = scale_d;
	 penalty_d = (2.0*lambda_d+1)*D2 / (gamma_d*D1);
       }
       else{
	 gamma_d = (2.0*lambda_d+1)*D2/(penalty_d*D1);
      }
     }else{
       lambda_d = penalty_d;
     }
     if( lambda_d <= 0 || penalty_d <= 0 ){
       cout << "Must specify parameters for shrinkage of covariates." << endl;
       exit(0);
     }
   }
   
   if( p != 0 ){
     cout << "Prior for SNP data: ";
     if( prior==0 ){
       cout << "NEG with parameters:" << endl;
       cout << "shape = " << lambda[0] << " scale = " << gamma[0]
	    <<  " penalty = " << penalty[0] << endl;
     }
     if( prior==1 )
       cout << "DE with parameter " << lambda[0] << endl;
     if( prior==2 )
       cout << "Gausaian with parameter " << lambda[0] << endl;
     if( standardize )
       cout << "With standardized data." << endl;
   }else
     cout << "No SNP data." << endl;
   
   cout << "Prior for data: ";
   if( prior_d==0 )
     cout << "NEG" << endl;
   if( prior_d==1 )
     cout << "DE" << endl;
   if( prior_d==2 )
     cout << "Gausaian" << endl;
   cout << "lambda = " << lambda_d << " gamma = " << gamma_d
	<<  " penalty = " << penalty_d << endl;
   if( standardize_d )
     cout << "With standardized data." << endl;
}

HLasso::~HLasso()
{
   outfile.close();
   delete [] beta_shrinkage;
   delete [] beta_type;
   delete [] beta_d;
}

void HLasso::permute_outcome()
{
   for( unsigned int j = 0; j < n; j++ ){
      short random = int(myrand()*(float)(n));
      if( logistic ){
         unsigned long temp = y[random];
         y[random] = y[j];
         y[j] = temp;
      }else{
         double temp = yy[random];
         yy[random] = yy[j];
         yy[j] = temp;
      }
   }
}

void HLasso::setOutputFiles( char* outfilename )
{
   cout << "Opening output file " << outfilename << endl;
   outfile.open(outfilename);
}

void HLasso::runCLG(unsigned short iSamples,
		    bool mle_start, bool permute )
{
  unsigned short iter, Samples = iSamples;
   vector<unsigned long> indicator, indicator_d, ind_comb;
   double LL, LP;
   vector <unsigned long> ptr_j, ptr_d;;
   
   shrinkage_fit.SetModel( prior, disease_models, beta0 );
   mle_fit.SetModel( 2, 1, beta0 );
   if( logistic ){
     mle_fit.SetResponse( n, y, missing_target );
     shrinkage_fit.SetResponse( n, y, missing_target );
   }else{
      mle_fit.SetResponse( n, yy, missing_target );
      shrinkage_fit.SetResponse( n, yy, missing_target );
   }

   shrinkage_fit.SetPriorParams( lambda, gamma );
   shrinkage_fit.SetPriorParams_d( lambda_d, gamma_d, prior_d );
//    mle_fit.SetPriorParams( .001, gamma );
//    mle_fit.SetPriorParams_d( .001, 0, 2 );

   for( unsigned int j = 0; j < p; j++ ){
     ptr_j.push_back(j);
   }
   for( unsigned int j = 0; j < p2; j++ ){
     ptr_d.push_back(j);
   }

   if( permute ){
     for( unsigned int j = 0; j < ptr_j.size(); j++ ){
       long random = int(myrand()*(float)(ptr_j.size()));
       unsigned long temp = ptr_j[random];
       ptr_j[random] = ptr_j[j];
       ptr_j[j] = temp;
     }
     for( unsigned int j = 0; j < ptr_d.size(); j++ ){
       long random = int(myrand()*(float)(ptr_d.size()));
       unsigned long temp = ptr_d[random];
       ptr_d[random] = ptr_d[j];
       ptr_d[j] = temp;
     }
   }

   for(  unsigned short kk = 0; kk < Samples; kk++ ){
      cout << "Iteration: " << kk << endl;
      for( unsigned int j = 0; j < p; j++ ){
	beta_shrinkage[j] = 0;
      }
      for( unsigned int j = 0; j < p2; j++ ){
	beta_d[j] = 0;
      }
      if( mle_start ){
	iter = mle_fit.Maximize( beta0, beta_shrinkage, beta_type,
				 beta_d, ptr_j, ptr_d, mean, sd,
				 X, missing, X_d, missing_d, LL );
	cout << "Found mle.\n";
      }
      iter = shrinkage_fit.Maximize( beta0, beta_shrinkage, beta_type,
				     beta_d, ptr_j, ptr_d, mean, sd,
				     X, missing, X_d, missing_d, LL );

      if( prior == 0 )
	LP = getNEGLogPrior( beta_shrinkage, p, lambda, gamma );
      else if( prior == 1 )
	LP = getLaplaceLogPrior( beta_shrinkage, p, lambda );
      else
	LP = getGausLogPrior( beta_shrinkage, p, lambda );
      
      if( prior_d == 0 )
	LP += getNEGLogPrior( beta_d, p2, lambda_d, gamma_d );
      else if( prior_d == 1 )
	LP += getLaplaceLogPrior( beta_d, p2, lambda_d );
      else
	LP += getGausLogPrior( beta_d, p2, lambda_d );
      
      for( unsigned int j = 0; j < p; j++ )
	if( beta_shrinkage[j] != 0.0 ){
	  indicator.push_back(j);
	}
      
      for( unsigned int j = 0; j < p2; j++ )
         if( beta_d[j] != 0.0 ){
            indicator_d.push_back(j);
         }
            
      pair< model_list::const_iterator, bool > new_element;
      ind_comb=indicator;
      for( unsigned int j = 0; j < indicator_d.size(); j++ ){
	ind_comb.push_back( indicator_d[j] + p );
      }
      outfile.setf(ios::fixed);
      new_element = ModelList.insert( model_list::value_type(ind_comb,LP+LL) );
      if( new_element.second ){
	if( ModelList.size() != 1 )
	  outfile << ",\n";
         outfile << "structure(c(";
         outfile << LL << "," << LL+LP << "," << -2.0*LL + double(indicator_d.size()+indicator.size())*log(double(n)) << ",0,";
	 float intercept = beta0*data_sd + data_mean;
         for( unsigned int j = 0; j < indicator_d.size(); j++ ){
	   unsigned int jj = indicator_d[j];
	   outfile << jj << ",\"" << Labels_d[jj] << "\","
		   <<scientific<< beta_d[jj]*data_sd/sd_d[jj] << "," << 0 << ",";
	   intercept -= beta_d[jj]*mean_d[jj]/sd_d[jj];
         }
         for( unsigned int j = 0; j < indicator.size(); j++ ){
	   unsigned int jj = indicator[j];
	   outfile << jj  << ",\"" << Labels[jj] << "\","
		   <<scientific<< data_sd*beta_shrinkage[jj]/sd[jj][beta_type[jj]]
		   << "," << beta_type[jj] << ",";
	   intercept -= beta_shrinkage[jj]*mean[jj][beta_type[jj]]/sd[jj][beta_type[jj]];
         }
	 intercept *= data_sd;
	 outfile << 0 << ",\"Intercept\"," << intercept << "," << 0;
         outfile << "), .Dim = c(4,"
		 << indicator.size()+indicator_d.size()+2 << "))" << endl;
      }
      indicator.clear();
      indicator_d.clear();
      if( kk < Samples-1 ){
         for( unsigned int j = 0; j < ptr_j.size(); j++ ){
            long random = int(myrand()*(float)(ptr_j.size()));
            unsigned long temp = ptr_j[random];
            ptr_j[random] = ptr_j[j];
            ptr_j[j] = temp;
         }
         for( unsigned int j = 0; j < ptr_d.size(); j++ ){
            long random = int(myrand()*(float)(ptr_d.size()));
            unsigned long temp = ptr_d[random];
            ptr_d[random] = ptr_d[j];
            ptr_d[j] = temp;
         }
      }
   }

   NumberOfOutputModels = ModelList.size();
   ModelList.clear();
}

void HLasso::OpenOutFiles()
{
   outfile << "structure(list(length =\n";
//    if(VSP){
//       outfile2 << "structure(list(length =\n";
//    }
}


void HLasso::CloseOutFiles()
{
   outfile << "),.Names = c(";
   for(unsigned int kk = 0; kk < NumberOfOutputModels-1; kk++ ){
      outfile << "\"\",";
   }
   outfile << "\"\"))\n";
}

double getNEGLogPrior( double *beta, unsigned long p, double lambda, double gamma )
{
   double v,x,f,df,l=0;
   v = -2.0*(lambda+0.5);
   for( unsigned int j = 0; j < p; j++ ){
      x = abs(beta[j])/gamma;
      pbdv_(&v,&x,&f,&df);
      l += log(f);
   }
   return l;
}

double getNEGLogPrior( double *beta, unsigned long p,
		       vector<double> lambda, vector<double> gamma )
{
   double v,x,f,df,l=0;
   for( unsigned int j = 0; j < p; j++ ){
     v = -2.0*(lambda[j]+0.5);
     x = abs(beta[j])/gamma[j];
     pbdv_(&v,&x,&f,&df);
     l += log(f);
   }
   return l;
}

double getLaplaceLogPrior( double *beta, unsigned long p, double lambda )
{
   double l=0;
   for( unsigned int j = 0; j < p; j++ ){
     l -= lambda*abs(beta[j]);
   }
   return l;
}

double getLaplaceLogPrior( double *beta, unsigned long p, vector<double> lambda )
{
   double l=0;
   for( unsigned int j = 0; j < p; j++ ){
     l -= lambda[j]*abs(beta[j]);
   }
   return l;
}

double getGausLogPrior( double *beta, unsigned long p, double lambda )
{
   double l=0;
   for( unsigned int j = 0; j < p; j++ ){
      l -= beta[j]*beta[j];
   }
   l *= 0.5*lambda;
   return l;
}

double getGausLogPrior( double *beta, unsigned long p, vector<double> lambda )
{
   double l=0;
   for( unsigned int j = 0; j < p; j++ ){
      l -= 0.5*lambda[j]*beta[j]*beta[j];
   }
   return l;
}

void HLasso::SetGW_BF( double iBF, const char *fname )
{
  vector<double> updateBF;
  vector<string> BFlabel;
  ifstream in(fname);
  double minProb = 0.5;
  if( in.is_open() )
    cout << "Loading " << fname << endl;

  if( *fname ){
    try{
      StringSplitter splitter;
      string line;
      while( getline(in, line) ){
	if( !isWhiteLine(line.c_str()) ){
	  Vector_s vline(splitter.split(line.c_str()));
	  BFlabel.push_back(vline[0]);
	  updateBF.push_back(StringConvertor::toFloat(vline[1]));
	}
	if( updateBF.back() < minProb )
	  minProb = updateBF.back();
      }
    }catch (...) {
      in.close();
      throw;
    }
  }
  BF.assign( p, iBF * (1.0-minProb)/minProb);
  
  for( unsigned int i = 0; i < Labels.size(); i++ ){
    unsigned int j=0;
    while( j < BFlabel.size() && Labels[i] != BFlabel[j] )
      j++;
    if( j < BFlabel.size() )
      BF[i] = iBF * (1.0 - updateBF[i]) / updateBF[i];
    else
       cout << "No weight for SNP " << Labels[i] << " in set Bayes Factor." << endl;
  }
}

void HLasso::SetBF( const char *fname )
{
  vector<double> updateBF;
  vector<string> BFlabel;
  ifstream in(fname);
  if( in.is_open() )
    cout << "Loading " << fname << endl;

  try{
    StringSplitter splitter;
    string line;
    while( getline(in, line) ){
      if( !isWhiteLine(line.c_str()) ){
	Vector_s vline(splitter.split(line.c_str()));
	BFlabel.push_back(vline[0]);
	updateBF.push_back(StringConvertor::toInt(vline[1]));
      }
    }
  }catch (...) {
    in.close();
    throw;
  }
  
  for( unsigned int i = 0; i < updateBF.size(); i++ ){
    unsigned int j=0;
    while( j < p && Labels[j] != BFlabel[i] )
      j++;
    if( j < p )
      BF[j] = updateBF[i];
    else
      cout << "Unknown SNP " << BFlabel[i] << " in set Bayes Factor." << endl;
  }
}

void HLasso::readFile2( const char *fname, const char *fname2, bool target )
{
  p = 0; p2 = 0;
  ifstream in(fname);
  if( in.is_open() )
    cout << "Loading " << fname << endl;
  ifstream in2(fname2);
  if( in2.is_open() )
    cout << "Loading " << fname2 << endl;
  if( !in.is_open() && !in2.is_open() ){
    cout << "Error: Must specify a genotypes or covariates file.\n";
    exit(0);
  }
  
  try {
    StringSplitter splitter;
    string line, line2;
    int j=0;
    
    while( getline(in, line) ){
      if( !isWhiteLine(line.c_str()) ){
	Vector_s vline(splitter.split(line.c_str()));
	if( j==0 ){
	  if( target )
	    vline.erase( vline.begin() );
	  Labels = vline;
	  p = vline.size();
	}else{
	  vector<float> vec;
	  if( vline.size()-target != p ){
	    cout << "Error: input file has only " << vline.size()
		 << " elements on line " << j << endl;
	    exit(0);
	  }
	  for( unsigned int i=0; i<vline.size(); i++ ){
	    if( target && i == 0 ){
	      if( logistic ){
		if( vline[i] == "0" || vline[i] == "1" ){
		  y.push_back(2*StringConvertor::toInt(vline[i])-1);
		  missing_target.push_back(false);
		}else if( vline[i] == "NA" ){
		  y[i] = 0;
		  missing_target.push_back(true);
		}else{
		  cout << "Controls coded as 0, cases coded as 1, do not recognize "
		       << vline[i] << endl;
		  exit(0);
		}
	      }
	      else
		if( vline[i] != "NA" ){
		  yy.push_back(StringConvertor::toFloat(vline[i]));
		  missing_target.push_back(false);
		}else{
		  yy[i] = 0;
		  missing_target.push_back(true);
		}
	    }
	    else if( vline[i] == "0" || vline[i] == "1" || vline[i] == "2" ){
	      X.push_back(StringConvertor::toInt(vline[i]));
	      missing.push_back(false);
	    }else if( vline[i] == "NA" || vline[i] == "9" || vline[i] == "-1" ){
	      X.push_back(0);
	      missing.push_back(true);
                }else{
	      cout << "Genotypes coded as 0, 1, 2, missing values coded as NA, 9, -1, do not recognize "
		   << vline[i] << endl;
	      exit(0);
	    }
	  }
	}
	j++;
      }
    }
    if( target )
      n = j-1;
    j = 0;
    while( getline(in2, line2) ){
      if( !isWhiteLine(line2.c_str()) ){
	Vector_s vline2(splitter.split(line2.c_str()));
	if( j==0 ){
	  p2 = vline2.size();
	  Labels_d = vline2;
	}else{
	  vector<float> vec;
	  if( vline2.size() != p2 ){
	    cout << "Error: covariates file has only " << vline2.size()
		 << " elements on line " << j << endl;
	    exit(0);
	  }
	  for( unsigned int i=0; i<vline2.size(); i++ ){
	    if( vline2[i] != "NA" ){
	      X_d.push_back(StringConvertor::toFloat(vline2[i]));
	      missing_d.push_back(false);
	    }else{
	      X_d.push_back(0);
	      missing_d.push_back(true);
	    }
	  }
	}
	j++;
      }
    }
  } catch (...) {
    in.close();
    throw;
 }
}
