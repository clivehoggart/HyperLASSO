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
double getNEGLogPrior( double, double, double );
double getLaplaceLogPrior( double, double );
double getGausLogPrior( double, double );
void StringToChar( char *c, string s );

HLasso::HLasso( bool ilogistic,
                char* infilename, char* targetfilename, char* covfilename,
		char* plinkfilename )
{
  X.resize(2);
  logistic = ilogistic;
  NumberOfOutputModels=0;
  cv_sets=100;
  cv_size=300;
  p2 = 0;
  Matrix_s target;

   if( *targetfilename ){
     readFile(targetfilename, target);
     cout << "Dimensions of target file: "
	  << target.size()  << " " << target[0].size() << endl;
     n = target.size();
     y.resize(n);
     yy.resize(n);
     string t;
     long observation = 0;
     for( unsigned int i = 0; i < n; i++ ){
       if( target[0].size() == 1 )
	 t = target[i][0];
       else
	 t = target[0][i];
       if( logistic ){
	 if( t == "0" || t == "1" ){
	   y[i] = 2*StringConvertor::toInt(t)-1;
	   missing_target.push_back(false);
	   observed_target.push_back(observation);
	 }else if( target[i][0] == "NA" ){
	   y[i] = 0;
	   missing_target.push_back(true);
	 }else{
	   cout << "Controls coded as 0, cases coded as 1, do not recognize "
		<< t << endl;
	   exit(0);
	 }
       }else{
	 if( t != "NA" ){
	   yy[i] = StringConvertor::toFloat(t);
	   y[i] = 1;
	   missing_target.push_back(false);
	   observed_target.push_back(observation);
	 }else{
	   yy[i] = 0;
	   y[i] = 0;
	   missing_target.push_back(true);
	 }
       }
       observation++;
     }
   }
   
   if( *plinkfilename )
     infilename = plinkfilename;
   readFile2( infilename, covfilename, !*targetfilename, *plinkfilename );
   SelectedRows = observed_target;
   if( p != 0 ){
     if( X[0].size() / p != n ){
       cout << "Dimensions of target and genotypes file do not match.\n";
       cout << "Genotypes file has " << X[0].size()/p << " rows." << endl;
       cout << "Target file has " << n << " rows." << endl;
       exit(0);
     }
   }
   if( p2 != 0 ){
     if( X_d.size() / p2 != n ){
       cout << "Dimensions of target and covariates file do not match.\n";
       cout << "Covariates file has " << X_d.size()/p2 << " rows." << endl;
       cout << "Target file has " << n << " rows." << endl;
       cout << "Genotypes file has " << X[0].size()/p << " rows." << endl;
       exit(0);
     }
   }
   cout << n << " "<< p <<  " " << p2 << endl;
   SumSelected.resize(p);
   SumSelected_d.resize(p2);
}

void HLasso::SetParams( unsigned short iprior, double ipenalty, double shape,
			double scale, unsigned short iprior_d, double ipenalty_d,
			double shape_d, double scale_d, bool standardize,
			bool standardize_d, short idisease_models, double alpha,
			double alpha_d, double iBF, char* GW_BFfile, char* BFfile )
{
  vector< vector<double> > penalty;
  if( alpha != 0 )
    standardize = true;
  disease_models = idisease_models;
  mean.resize(p);
  info.resize(p);
  penalty.resize(p);
  lambda.resize(p);
  gamma.resize(p);
   sd.resize(p);
   double s2;
   unsigned short xx, XX;
   unsigned long obs=0;
   for( unsigned int j = 0; j < p; j++ ){
     cout << "Reading SNP " << j << endl;
      mean[j].resize(disease_models);
      info[j].resize(disease_models);
      sd[j].resize(disease_models);
      penalty[j].resize(disease_models,ipenalty);
      lambda[j].resize(disease_models,shape);
      gamma[j].resize(disease_models,scale);
      for( short type = 0; type < disease_models; type++ ){
	s2 = 0.0;
         obs = 0;
	 mean[j][type] = 0.0;
	 info[j][type] = 0.0;
	 for( unsigned long ii = 0; ii < SelectedRows.size(); ii++ ){
	   unsigned long i = SelectedRows[ii];
            if( !missing[i*p+j] ){
               obs++;
	       XX = (short)X[0][i*p+j] + (short)X[1][i*p+j];
               if( type == 0 )
                  xx = XX;
	       else if( type == 1 ){
		 if( XX == 2 )
		   xx = 2;
		 else
		   xx = 0;
	       }else if( type == 2 ){
		 if( XX == 0 )
		   xx = 0;
		 else
		   xx = 2;
	       }else{
		 if( XX == 1 )
		   xx = 2;
		 else
		   xx = 0;
	       }
               mean[j][type] += xx;
               s2 += xx*xx;
            }
         }
         mean[j][type] /= obs;
	 s2 = sqrt(s2/obs - mean[j][type]*mean[j][type]);
	 cout << s2 << endl;
	 if( s2 == 0 ){
	   cout << "Warning, variable " << Labels[j] << " at position " << j  << " is monomorphic for effect type " << type << endl;
	   sd[j][type] = 1;
	 }else if( isnan(s2) ){
	   cout << "Warning, variable " << Labels[j] << " at position " << j  << " has all missing values." << endl;
	   mean[j][type] = 0;
	   sd[j][type] = 1;
	 }else
	   sd[j][type] = s2;
      }
   }
   mean_d.resize(p2);
   sd_d.resize(p2);
   for( unsigned int j = 0; j < p2; j++ ){
     mean_d[j] = 0.0;
     sd_d[j] = 0.0;
     obs=0;
     for( unsigned long ii = 0; ii < SelectedRows.size(); ii++ ){
       unsigned long i = SelectedRows[ii];
       if( !missing_d[i*p2+j] ){
	 obs++;
	 mean_d[j] += X_d[i*p2+j];
	 sd_d[j] += X_d[i*p2+j]*X_d[i*p2+j];
       }
     }
     mean_d[j] /= obs;
     sd_d[j] = sqrt(sd_d[j]/obs - mean_d[j]*mean_d[j]);
     if( sd_d[j] == 0 ){
       cout << "Warning, variable " << Labels_d[j] << " at position " << j  << " is monomorphic, therefore cannot standardize." << endl;
       sd_d[j] = 1;
     }
     if( !standardize_d )
       sd_d[j] = 1.0;
     for( unsigned long ii = 0; ii < SelectedRows.size(); ii++ ){
       unsigned long i = SelectedRows[ii];
       X_d[i*p2+j] = (X_d[i*p2+j]-mean_d[j])/sd_d[j];
     }
   }

   beta_shrinkage = new double[p];
   beta_type = new short[p];
   long n0=0,n1=0;
   double sum=0,sum2=0;
   obs = 0;
   if( logistic ){
      data_mean = 0;
      data_sd = 1;
      for( unsigned long ii = 0; ii < SelectedRows.size(); ii++ ){
	unsigned long i = SelectedRows[ii];
         if( y[i]==1 )
	   n1++;
	 else if( y[i]==-1 )
	   n0++;
      }
      beta0 = log((double)n1/(double)n0);
      if( alpha != 0 ){
	for( unsigned int j = 0; j < p; j++ )
	  penalty[j].assign( disease_models, sqrt( (double)(n0*n1)/(n0+n1) ) * gsl_cdf_ugaussian_Pinv( 1-alpha/2) );
      }
   }else{
     for( unsigned long ii = 0; ii < SelectedRows.size(); ii++ ){
       unsigned long i = SelectedRows[ii];
       sum += yy[i];
       sum2 += yy[i]*yy[i];
     }
     data_mean = sum/SelectedRows.size();
     data_sd = sqrt(sum2/SelectedRows.size() - data_mean*data_mean);
     cout << data_mean << " " << SelectedRows.size() << endl;
     cout << "sd: " << data_sd << endl;
     for( unsigned long ii = 0; ii < SelectedRows.size(); ii++ ){
       unsigned long i = SelectedRows[ii];
       yy[i] = (yy[i]-data_mean)/data_sd;
     }
     beta0 = 0;
     if( alpha != 0 ){
       for( unsigned int j = 0; j < p; j++ )
	 penalty[j].assign( disease_models, sqrt(n) * gsl_cdf_ugaussian_Pinv( 1-alpha/2) );
     }
   }

   for( unsigned int j = 0; j < p; j++ ){
     for( short type = 0; type < disease_models; type++ ){
       if( logistic ){
	 info[j][type] = (double)(n0+n1)/( (double)n0 * (double)n1 * sd[j][type] * sd[j][type] );
       }else{
	 info[j][type] = 1 / (SelectedRows.size() * sd[j][type] * sd[j][type] );
       }
     }
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
	for( short type = 0; type < disease_models; type++ ){
	  double k=1;
	  if( mean[j][type] > 2 ){
	    k *= 100;
	  }
	  penalty[j][type] = ( k * sqrt(2.0 * log(BF[j] * sqrt(1.0+w/info[j][type])) * (1.0 + info[j][type]/w) / info[j][type]) );
	}
     }
   }

   cout << "Intercept: " << beta0 << endl;
   for( unsigned int j = 0; j < p; j++ ){
      beta_type[j] = 0;
      beta_shrinkage[j] = 0;
      if( !standardize )
	for( short type = 0; type < disease_models; type++ )
	  sd[j][type] = 1.0;
   }
   beta_d = new double[p2];
   for( unsigned int j = 0; j < p2; j++ ){
      beta_d[j] = 0;
   }

   prior = iprior;
   prior_d = iprior_d;
   if( penalty.size() == 0 )
     for( unsigned int j = 0; j < p; j++ )
       penalty[j].assign( disease_models, ipenalty );
   penalty_d = ipenalty_d;

   if( prior==0 ){
     for( unsigned int j = 0; j < p; j++ )
       lambda[j].assign( disease_models, shape );
     double x=0, v, temp, D1, D2;
     if( scale != 0 ){
       v = -2.0*lambda[0][0]-1.0;
       pbdv_(&v,&x,&D1,&temp);
       v--;
       pbdv_(&v,&x,&D2,&temp);
       for( unsigned int j = 0; j < p; j++ )
	 gamma[j].assign( disease_models, scale );
     }
     else{
       for( unsigned int j = 0; j < p; j++ ){
	 for( short type = 0; type < disease_models; type++ ){
	   v = -2.0*lambda[j][type]-1.0;
	   pbdv_(&v,&x,&D1,&temp);
	   v--;
	   pbdv_(&v,&x,&D2,&temp);
	   gamma[j][type] = ( (2.0*lambda[j][type]+1)*D2 / (penalty[j][type]*D1) );
	   if( isnan(gamma[j][type]) ){
	     cout << "bollox" << endl;
	   }
	 }
       }
     }
   }else{
     lambda = penalty;
   }
   
   if( p2 != 0 ){
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
	   penalty_d = sqrt( (double)(n0*n1)/(n0+n1) ) * gsl_cdf_ugaussian_Pinv( 1-alpha_d/2);
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
       cout << "shape = " << lambda[0][0] << " scale = " << gamma[0][0]
	    <<  " penalty = " << penalty[0][0] << endl;
     }
     if( prior==1 )
       cout << "DE with parameter " << lambda[0][0] << endl;
     if( prior==2 )
       cout << "Gausaian with parameter " << lambda[0][0] << endl;
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
   selectionprobfile.close();
   delete [] beta_shrinkage;
   delete [] beta_type;
   delete [] beta_d;
}

void HLasso::permute_outcome()
{
  bool tmp, maintain_pheno_cov_order=true;
  for( unsigned int i = 0; i < n; i++ ){
    short random = int(myrand()*(float)(n));
    if( logistic ){
      unsigned long temp = y[random];
      y[random] = y[i];
      y[i] = temp;
    }else{
      double temp = yy[random];
      yy[random] = yy[i];
      yy[i] = temp;
    }
    if( maintain_pheno_cov_order ){
      for( unsigned int j = 0; j < p2; j++ ){
	double temp = X_d[random*p2+j];
	X_d[random*p2+j] = X_d[i*p2+j];
	X_d[i*p2+j] = temp;
	tmp = missing_d[random];
	missing_d[random] = missing_d[i];
	missing_d[i] = tmp;
      }
    }
  }
}

void HLasso::runCLG(unsigned short iSamples,
		    bool mle_start, bool permute, bool ifirst_model )
{
  unsigned short iter, Samples = iSamples;
  vector<unsigned long> indicator, indicator_d, ind_comb;
  double LL, LP;
  vector <unsigned long> ptr_j, ptr_d;;

  bool first_model;
  if( ifirst_model )
    first_model=true;
  else
    first_model=false;

   shrinkage_fit.SetModel( prior, disease_models, beta0 );
   //   mle_fit.SetModel( 2, 1, beta0 );
   if( logistic ){
     //     mle_fit.SetResponse( y, SelectedRows );
     shrinkage_fit.SetResponse( y, SelectedRows );
   }else{
     //      mle_fit.SetResponse( yy, SelectedRows );
      shrinkage_fit.SetResponse( yy, SelectedRows );
   }

   shrinkage_fit.SetPriorParams( lambda, gamma, info, BF.size()!=0 );
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
     LP = 0.0;
      cout << "Iteration: " << kk << endl;
      for( unsigned int j = 0; j < p; j++ ){
	beta_shrinkage[j] = 0;
      }
      for( unsigned int j = 0; j < p2; j++ ){
	beta_d[j] = 0;
      }
      // if( mle_start ){
      // 	iter = mle_fit.Maximise( beta0, beta_shrinkage, beta_type,
      // 				 beta_d, ptr_j, ptr_d, mean, sd,
      // 				 X, missing, X_d, missing_d, LL );
      // 	cout << "Found mle.\n";
      // }
      iter = shrinkage_fit.Maximise( beta0, beta_shrinkage, beta_type,
				     beta_d, ptr_j, ptr_d, mean, sd,
				     X, missing, X_d, missing_d, LL );

      if( prior == 0 )
	for( unsigned int j = 0; j < p; j++ ){
	  LP += getNEGLogPrior( beta_shrinkage[j], lambda[j][beta_type[j]], gamma[j][beta_type[j]] );
	}
      if( prior == 1 )
	for( unsigned int j = 0; j < p; j++ )
	  LP += getLaplaceLogPrior( beta_shrinkage[j], lambda[j][beta_type[j]] );
      if( prior == 2 )
	for( unsigned int j = 0; j < p; j++ )
	  LP += getGausLogPrior( beta_shrinkage[j], lambda[j][beta_type[j]] );
      
      // if( prior_d == 0 )
      // 	LP += getNEGLogPrior( beta_d, p2, lambda_d, gamma_d );
      // else if( prior_d == 1 )
      // 	LP += getLaplaceLogPrior( beta_d, p2, lambda_d );
      // else
      // 	LP += getGausLogPrior( beta_d, p2, lambda_d );
      
      for( unsigned int j = 0; j < p; j++ )
	if( beta_shrinkage[j] != 0.0 ){
	  indicator.push_back(j);
	  SumSelected[j]++;
	}
      
      for( unsigned int j = 0; j < p2; j++ )
	if( beta_d[j] != 0.0 ){
	  indicator_d.push_back(j);
	  SumSelected_d[j]++;
	}
            
      pair< model_list::const_iterator, bool > new_element;
      ind_comb=indicator;
      for( unsigned int j = 0; j < indicator_d.size(); j++ ){
	ind_comb.push_back( indicator_d[j] + p );
      }
      //      outfile.setf(ios::fixed);
      new_element = ModelList.insert( model_list::value_type(ind_comb,LP+LL) );
      if( new_element.second ){
	if( !first_model ){
	  outfile << ",\n";
	}
	outfile << "structure(c(";
	outfile << LL << "," << LL+LP << "," << iter << ",0,";
	float intercept = beta0*data_sd + data_mean;
	for( unsigned int j = 0; j < indicator_d.size(); j++ ){
	  unsigned int jj = indicator_d[j];
	  outfile << jj << ",\"" << Labels_d[jj] << "\","
		  <<scientific<< beta_d[jj]*data_sd/sd_d[jj] << "," << 0 << ",";
	  intercept -= beta_d[jj]*mean_d[jj]*data_sd/sd_d[jj];
	}
	for( unsigned int j = 0; j < indicator.size(); j++ ){
	  unsigned int jj = indicator[j];
	  outfile << jj  << ",\"" << Labels[jj] << "\","
		  <<scientific<< data_sd*beta_shrinkage[jj]/sd[jj][beta_type[jj]]
		  << "," << beta_type[jj] << ",";
	  intercept -= beta_shrinkage[jj] * mean[jj][beta_type[jj]] * data_sd/sd[jj][beta_type[jj]];
	}
	outfile << 0 << ",\"Intercept\"," << intercept << "," << 0;
	outfile << "), .Dim = c(4,"
		<< indicator.size()+indicator_d.size()+2 << "))" << endl;
	first_model=false;
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

void HLasso::WriteSelectionProbs(unsigned short perms)
{
  for( unsigned long i=0; i<p; i++ ){
    selectionprobfile << Labels[i] << " " << (double)SumSelected[i]/perms << endl;
  }
}

void HLasso::OpenSelectionProbFile(char* outfilename)
{
  char filename[256];
  string str_outfile =  outfilename;
  
  str_outfile = str_outfile + "_SelectionProbs.dat";
  StringToChar( filename, str_outfile );
  cout << "Opening selection prob file " << filename << endl;
  selectionprobfile.open(filename);
}

void HLasso::OpenOutFiles(char* outfilename)
{
  char filename[256];
  string str_outfile =  outfilename;
  
  str_outfile = str_outfile + "_models.R";
  StringToChar( filename, str_outfile );
  cout << "Opening output file " << filename << endl;
  outfile.open(filename);
  outfile << "structure(list(length =\n";
}


void HLasso::CloseOutFiles(unsigned int samples)
{
   outfile << "),.Names = c(";
   for(unsigned int kk = 0; kk < (samples-1); kk++ ){
      outfile << "\"\",";
   }
   outfile << "\"\"))\n";
}

void HLasso::SampleRows()
{
  unsigned long n = observed_target.size(), size = (observed_target.size()/2);
  vector<unsigned long> vec = sample( n, size );
  SelectedRows.resize(size);
  for( unsigned long i=0; i<size; i++ ){
    SelectedRows[i] = observed_target[vec[i]];
  }
}

double getNEGLogPrior( double beta, double lambda, double gamma )
{
  double v,x,ff,f,df;
   v = -2.0*lambda - 1.0;
   x = abs(beta)/gamma;
   pbdv_(&v,&x,&ff,&df);
   f = log(ff);
   return f;
}

double getLaplaceLogPrior( double beta, double lambda )
{
   double l = -lambda*abs(beta);
   return l;
}

double getGausLogPrior( double beta, double lambda )
{
   double l = -0.5*lambda*beta*beta;
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

void HLasso::readFile2( const char *fname, const char *fname2, bool target, bool plink )
{
  p = 0; p2 = 0;
  ifstream in(fname);
  if( in.is_open() )
    cout << "Loading " << fname << endl;
  else if(*fname){
    cout << "Cannot open " << fname << endl;
    exit(0);
  }
  ifstream in2(fname2);
  if( in2.is_open() )
    cout << "Loading " << fname2 << endl;
  else if(*fname2){
    cout << "Cannot open " << fname2 << endl;
    exit(0);
  }
  if( !in.is_open() && !in2.is_open() ){
    cout << "Error: Must specify a genotypes or covariates file.\n";
    exit(0);
  }
  // cout << "max_size: " << X.max_size() << "\n";
  // cout << "max_size: " << missing.max_size() << "\n";
  
  try{
    StringSplitter splitter;
    string line, line2;    
    unsigned long observation=0;
    while( getline(in, line) ){
      if( !isWhiteLine(line.c_str()) ){
	Vector_s vline(splitter.split(line.c_str()));
	//	cout << "max_size: " << vline.max_size() << "\n";
	if( observation==0 ){
	  if( plink )
	    vline.erase( vline.begin(), vline.begin()+6 );
	  if( target && !plink )
	    vline.erase( vline.begin() );
	  Labels = vline;
	  p = vline.size();
	}else{
	  if( plink )
	    vline.erase( vline.begin(), vline.begin()+5 );
	  if( !target && plink )
	    vline.erase( vline.begin() );
	  if( (vline.size() - target) != p ){
	    cout << "Error: input file has only " << vline.size()
		 << " elements on line " << observation << endl;
	    exit(0);
	  }
	  vector<float> vec;
	  for( unsigned long i=0; i<vline.size(); i++ ){
	    if( target && i == 0 ){
	      short pheno;
	      if( logistic ){
		if( vline[i] == "0" || vline[i] == "1" || ( (vline[i] == "2") && plink) ){
		  pheno = 2*StringConvertor::toInt(vline[i])-1;
		  if( plink )
		    pheno -= 2;
		  y.push_back(pheno);
		  missing_target.push_back(false);
		  observed_target.push_back(observation-1);
		}else if( vline[i] == "NA" ){
		  y.push_back(0);
		  missing_target.push_back(true);
		}
		if( pheno !=-1 && pheno != 1 ){
		  cout << "Controls coded as 0, cases coded as 1, do not recognize "
		        << vline[i] << endl;
		  exit(0);
		}
	      }else{
		if( vline[i] != "NA" ){
		  yy.push_back(StringConvertor::toFloat(vline[i]));
		  missing_target.push_back(false);
		  observed_target.push_back(observation-1);
		}else{
		  yy.push_back(0);
		  missing_target.push_back(true);
		}
	      }
	    }else if( vline[i] == "NA" || vline[i] == "9" || vline[i] == "-1" ){
	      X[0].push_back(false);
	      X[1].push_back(false);
	      missing.push_back(true);
	    }
	    else{
	      missing.push_back(false);
	      if( vline[i] == "0" ){
		X[0].push_back(false);
		X[1].push_back(false);
	      }else if( vline[i] == "1" ){
		X[0].push_back(true);
		X[1].push_back(false);
	      }else if( vline[i] == "2" ){
		X[0].push_back(true);
		X[1].push_back(true);
	      }else{
		cout << "Genotypes must be coded 0,1 or 2, or NA, 9 or -1 for missing. Do not recognize "
		     << vline[i] << endl;
		exit(0);
	      }
	    }
	  }
	}
	observation++;
      }
    }
    int observation2=0;
    while( getline(in2, line2) ){
      if( !isWhiteLine(line2.c_str()) ){
    	Vector_s vline2(splitter.split(line2.c_str()));
	if( plink )
	  vline2.erase( vline2.begin(), vline2.begin()+2 );
    	if( observation2==0 ){
    	  p2 = vline2.size();
    	  Labels_d = vline2;
    	}else{
    	  vector<float> vec;
    	  if( vline2.size() != p2 ){
    	    cout << "Error: covariates file has only " << vline2.size()
    		 << " elements on line " << observation2 << endl;
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
    	observation2++;
      }
    }
    if( target )
      n = observation-1;
  } catch (...) {
    cout << X.size() << endl;
    in.close();
    throw;
  }
}

void StringToChar( char *c, string s )
{
   int len;
   len = s.length();
   s.copy( c, len, 0 );
   c[ len ] = 0;
}
