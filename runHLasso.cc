/*
 *  This file is part of the HyperLasso program.
 *  Copyright (c) Clive J Hoggart   <c.hoggart@imperial.ac.uk>
 *
 *  HyperLasso is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License, version 2,
 *  as published by the Free Software Foundation.
 *
 *  HyperLasso is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 *  02111-1307 USA
 *
 */ 

#include "HLasso.h"

using namespace std;

int main( int argc, char* argv[] )
{
  bool logistic=true, std=false, std_d=false, mle_start=false, permute=false;
  unsigned short model=0, model_d=0, null=0, disease_models=1, stability=0;
   unsigned int iter=1;
   double penalty=0.0, scale=0.0, shape=0.0, penalty_d=0.0, scale_d=0.0, shape_d=0.0, step = 0.0, alpha=0.0, alpha_d=0.0, BF = 0.0;
   char infile[256], targetfile[256], outfile[256], infile_d[256], modefile[256], BFfile[256], GW_BFfile[256], plinkfile[256];
   *targetfile=0;
   *infile=0;
   *BFfile=0;
   *GW_BFfile=0;
   *infile_d=0;
   *plinkfile=0;
   float threshold = 0;
   
   int na=1;
   while(na < argc)
   {
      if ( 0 == strcmp(argv[na],"-penalty") ){
         penalty = atof(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-penalty_d") ){
         penalty_d = atof(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-threshold") ){
         threshold = atof(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-shape") ){
         shape = atof(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-scale") ){
         scale = atof(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-shape_d") ){
         shape_d = atof(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-scale_d") ){
         scale_d = atof(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-step") ){
         step = atof(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-iter") ){
         iter = atoi(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-plink") ){
         strcpy(plinkfile,argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-genotypes") ){
         strcpy(infile,argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-target") ){
         strcpy(targetfile,argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-covariates") ){
         strcpy(infile_d,argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-o") ){
         strcpy(outfile,argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-mode") ){
         strcpy(modefile,argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-linear") ){
         logistic = false;
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-std") ){
         std = true;
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-std_d") ){
         std_d = true;
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-dom") ){
         disease_models=3;
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-dom+het") ){
         disease_models=4;
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-null") ){
         null = atoi(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-stability") ){
	stability = atoi(argv[++na]);
	if (na+1==argc) break;
	na++;
      }
      else if ( 0 == strcmp(argv[na],"-mlestart") ){
         mle_start = true;
         if (na+1==argc) break;
         na++;
      }
// fit = 0 -> Hyper-Lasso
// fit = 1 -> Lasso
// fit = 2 -> Ridge
      else if ( 0 == strcmp(argv[na],"-model") ){
         model = atoi(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-model_d") ){
         model_d = atoi(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-seed")){
         smyrand( (unsigned long)(atoi(argv[++na])) );
         if (na+1==argc) break;
         na++;
	 permute = true;
      }
      else if ( 0 == strcmp(argv[na],"-typeIerr")){
         alpha = atof(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-BF")){
         BF = atof(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-BFfile")){
	strcpy(BFfile,argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-GW_BFfile")){
	strcpy(GW_BFfile,argv[++na]);
         if (na+1==argc) break;
         na++;
      }
      else if ( 0 == strcmp(argv[na],"-typeIerr_d")){
         alpha_d = atof(argv[++na]);
         if (na+1==argc) break;
         na++;
      }
//       else if ( 0 == strcmp(argv[na],"-startvals")){
//          Matrix_s start;
//          readFile( argv[++na], start );
//          for( unsigned long i = 0; i < start.size(); i++ ){
//             ptr_j.push_back( StringConvertor::toInt(start[i][0]) );
//          }
//          if (na+1==argc) break;
//          na++;
//       }
      else{
         cout << "Unknown option: " << argv[na] << endl;
         exit(0);
      }
   }

   if( shape > 25 )
     model = 1;
   if( shape_d > 25 )
     model_d = 1;
   
   if( strlen(outfile) != 0 ){
     HLasso HL( logistic, infile, targetfile, infile_d, plinkfile );
     HL.SetParams( model, penalty, shape, scale,
		   model_d, penalty_d, shape_d, scale_d,
		   std, std_d, disease_models, alpha, alpha_d,
		   BF, GW_BFfile, BFfile );
     short samples=0;
     if( null != 0 )
       samples = null;
     else if( stability != 0 ){
       samples = stability;
       HL.OpenSelectionProbFile(outfile);
     }
     else
       samples = 1;
     HL.OpenOutFiles(outfile);
     bool first_model=true;
     for( int i = 0; i < samples; i++ ){
       cout << i << endl;
       if( null != 0 )
	 HL.permute_outcome();
       if( stability != 0 ){
	 HL.SampleRows();
	 HL.SetParams( model, penalty, shape, scale,
		       model_d, penalty_d, shape_d, scale_d,
		       std, std_d, disease_models, alpha, alpha_d,
		       BF, GW_BFfile, BFfile );
       }
       HL.runCLG( iter, mle_start, permute, first_model );
       first_model=false;
       permute=true;
     }
     if( stability != 0 )
       HL.WriteSelectionProbs(stability);
     int ii = samples*iter;
     HL.CloseOutFiles(ii);
   }
}
