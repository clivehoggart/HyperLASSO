# HyperLASSO
Simultaneous analysis of many SNPs and covariates using penalised regression
----------------------------------------------------------------------------

Contents of this file:

1. Program description
2. Download, Installation and Running
3. Examples

1. Program description
-----------------------

The runHLasso program implements the HyperLasso, or HLasso for short,
an algorithm that finds approximate modes of a penalised likelihood
function for linear or logistic regression, when the penalty function
is the Normal-Exponential-Gamma (NEG) probability density.
Equivalently, the HLasso finds posterior modes in Bayesian regression
with the NEG prior.  Because only the posterior mode is found, the
HLasso does not implement full Bayesian inference, but the Bayesian
language of "prior" and "posterior" can provide intuitive alternatives
to "penalty" and "penalised likelihood".

The NEG has two parameters, shape and scale.  In the limit as both
increase such that sqrt(2*shape)/scale is a constant, say lambda, the
NEG converges to the Laplace (double exponential) distribution with
rate parameter lambda.  In this special case the HLasso reduces to the
Lasso procedure of: Tibshirani R (1996) Regression shrinkage and
selection via the lasso. Journal of the Royal Statistical Society,
series B 58: 267-88.

HLasso is useful when there are many more predictors than
observations, because the penalty function can overcome the problem of
over-fitting that would undermine any usefulness of standard logistic
regression.  When the shape parameter is small (say < 1; we find 0.05
to be the smallest practical value, below this computational problems
arise), the NEG penalty function has the properties: (a) sharp peak at
zero, reflecting a strong "penalty" on any non-zero regression
coefficient; and (b) a very high gradient near zero, but heavy
(i.e. flattish) tails away from zero, corresponding to little
shrinkage effect on any non-zero value.  In Bayesian language, this
corresponds to a strong prior on effect sizes very close to zero, but
little prior information about the magnitude of non-zero effects.  In
contrast, the Laplace penalty implemented in the Lasso has lighter
tails, with more curvature away from the origin, so that non-zero
parameter values tend to be strongly "shrunk" towards zero, reflecting
a prior belief that non-zero effects are small.  This has the
undesirable consequence that if there are two strongly-correlated
predictors, only one of which is causal, the mode may correspond to
including both predictors in the model so that the causal effect is
spread over the two, and its size is thus greatly underestimated at
the true causal predictor.  The HLasso is more inclined to choose one
of the two predictors, and unless the correlation between the
predictors is perfect this will usually be the true causal predictor.

Although it can also solve other regression problems with many
predictors, HLasso has been designed to analyse genome-wide genetic
association studies (GWAS) in which the predictors are SNP genotypes
coded as 0, 1, and 2.  The case-control phenotype is the default, but
HLasso can also accommodate a continuous phenotype that is modelled
using the usual linear regression assumption of Gaussian (i.e. normal)
residuals.  The situation described above, with highly-correlated
predictors, often arises with tightly-linked SNPs in a GWAS.  If
several such SNPs are associated with phenotype, it could be that
there are multiple tightly-linked causal SNPs, or that there is a
single causal SNP and the others display association only through LD
with it.  HLasso can typically distinguish these scenarios
automatically, and include in the fitted model either multiple SNPs,
or just one, as appropriate.

A note on data input matrices:

The covariates for the regression analyses are supplied to HLasso in
two files, intended to be the genetic (i.e. SNP) and non-genetic
covariates respectively; see below for format.  However, for many
purposes HLasso treats all covariates the same way, and in effect
appends the non-genetic covariates after the genetic covariates.  The
distinction between the two files is that the genetic covariates are
stored as short integers, and specific genetic models are applied
assuming three distinct values; non-genetic covariates are stored as
double-precision reals and are only analysed as linear covariates.  A
feature of HLasso is that it can cope with up to 1 million SNPs in a
few hours of computation time per posterior mode estimate; we
typically evaluate 100 mode estimates.

Credits and further information:

runHLasso was written by Clive Hoggart with input from the other
authors of the paper below.  Funding came from the UK Medical Research
Council.  For further details and the results of simulation and real
data studies illustrating its power, see:

Hoggart CJ, Whittaker JC, De Iorio M, Balding DJ (2008) Simultaneous
Analysis of all SNPs in Genome-Wide and Re-sequencing Association
Studies, PLoS Genetics 4(7): e1000130;
doi:10.1371/journal.pgen.1000130

If this document and the PLoS Genetics paper do not answer your
queries, please contact Clive c.hoggart@imperial.ac.uk or David
Balding d.balding@imperial.ac.uk


2. Download, Installation and Running
-------------------------------------

The HLasso code will compile under the Linux OS.  It requires that
gsl, gcc and g77 are already installed.  First download HLasso.tar.gz
from the EBI website. Then to install type:

	tar -zxf HLasso.tar.gz
	cd HyperLasso
	make

HLasso is run from the command line.  We first list all possible
commandline options, then give commands to run three example analyses
of the toy data set (gen.dat, cov.dat, phen.dat, cont_phen.dat)
provided with the code in directory data.

-genotypes Path name for genotypes file. Rows index individuals and
	columns index SNPs. Genotypes are coded as 0, 1 and 2. First
	row contains SNP names, but there is no column for individual
	identifiers. Columns are separated by any white space.
	Missing values can be coded as NA, 9 or -1.

-covariates Path name for file containing non-genetic covariates. Rows
	index individuals and columns index covariates. All covariates
	are analysed as linear (factor covariates with more than 2
	levels must be coded as multiple binary covariates).  Missing
	values are coded as above.

-target Path name for file containing the (single) phenotype variable,
	can be either row or column.  Unless -linear is set, only 0 and
	1 values are allowed.  If -linear is set, HLasso automatically
	standardises the phenotype variable.

-linear Switch to fit linear regression model (continuous phenotype);
	default is logistic regression (case-control phenotype). If
	selected, the model has an additional parameter: the precision
	tau (=1/variance) of the target variable, conditional on the
	covariates. Initially tau=1, but is updated to its maximum
	likelihood estimate after each cycle through the covariates,
	and the value is output at each iteraction.

-iter   Number of mode searches. Default is 1. The search procedure is
	deterministic for a given order of covariates (both genetic
	and non-genetic), but the mode identified can vary according
	to this search order. The first iteration uses the search
	order supplied by the user (genetic covariates then
	non-genetic covariates, both in the order of input).
	Subsequent iterations use a random permutation of all
	covariates.  The search order remains fixed through each
	iteration, i.e. as the algorithm cycles through the covariates
	until a convergence criterion is met; we set this criterion in
	the same way as: Genkin A, Lewis DD, Madigan D (2007)
	Large-scale Bayesian logistic regression for text
	categorization. Technometrics 49(3): 291-304.

-o	Path for output file.  Output is an R list which can be read
	into R using the dget command. Each element of the list is a
	matrix describing a unique mode found by the search. First
	column records the number of cycles required to find the mode,
	the log-posterior and the log-likelihood of the mode (the
	final entry of this column is always 0). Subsequent columns
	describe the covariates selected in the model, first the
	non-genetic and then the genetic covariates:

	row 1 - position (rank) of the covariate in eah of the input
	    files, starting at zero, for nongenetic and then genetic
	    covariates. 
	row 2 - the covariate name supplied in the input files.
	row 3 - the value of the regression coefficient (recall that
	    continuous phenotypes are standardised by HLasso; to
	    recover coefficient estimates for the original scale,
	    multiply by the phenotype standard deviation).
	row 4 - for non-genetic covariates, always 0.  For genetic
	    covariates, indicator for effect type (i.e. the genetic
	    model):

			  0 - additive (genotypes coded 0,1,2),
			  1 - recessive (genotypes coded 0,0,1),
			  2 - dominant (genotypes coded 0,1,1),
			  3 - heterozygous (genotypes coded 0,1,0).

-dom	Switch to include tests for dominant, recessive and
	heterozygous effects for the genetic covariates (additive
	effects are always tested). If set, for a SNP that is not
	currently included in the model, the four effect types shown
	above are all considered, and only the one with the largest
	gradient of the log-likelihood at the origin is considered for
	entry into the model.  If included in the model, only the
	parameter value is updated at subsequent visits of the HLasso
	to this SNP, and not the effect type; the latter can only
	change within an iteration by the parameter value reverting to
	zero, and then a different model being selected, in subsequent
	cycles.

-std 	Switch to standardise the genetic covariates. Default is not to
   	standardize.

-std_d  Switch to standardise the non-genetic covariates. Default is
   	not to standardize (the d denotes that these can be double
   	precision).

-model 	Indicates which shrinkage prior is being used for the genetic
	covariates: 0 NEG (default) 1 DE 2 Gaussian.

-model_d As above but for the non-genetic covariates.

-shape	Shape parameter of NEG prior for genotypes.

-shape_d As above but for the non-genetic covariates.

-penalty Either the rate parameter of a DE prior, or precision
	of a normal prior (also called tau, =1/variance). If used
	with NEG prior, the scale parameter of NEG will be set such
	that the derivative of the density at the origin will be the
	same as for a DE distribution with this rate parameter.

-penalty_d As above but for the non-genetic covariates.

-scale Scale parameter of NEG prior for genotype data, not required if
       lambda is set.

-scale_d As above but for the non-genetic covariates.


3. Examples
------------

These examples use the following data files that are supplied with the
code in directory data:

gen.dat Genotype SNP data for 1000 individuals at 100 independent SNPs.

cov.dat	Data for three covariates for the 1000 individuals.

phen.dat Binary outcome for the 1000 individuals generated from the following model:

P(Affected) = logit( alpha + SNP1 + SNP2/2 + SNP3/4 + SNP4/8 + cov1/2 + cov2/4 +cov3/8 )
(alpha chosen to give approximately equal numbers of cases and controls)

cont_phen.dat Continuous outcome for the 1000 individuals generated
	      from the following model:

Trait ~ N( SNP1 + SNP2/2 + SNP3/4 + SNP4/8 + cov1/2 + cov2/4 +cov3/8, 1 )

These datasets are not in any sense realistic, and are only used to
show how to run the program and interpret its output, not to demonstrate
the power of HLasso.  Note that # SNPs << # individuals, and so
standard logistic regression should work well - we use this for
comparison below.

Example 1:

./runHLasso -genotypes data/gen.dat -target data/phen.dat -shape 0.1 -lambda 50 -o test.R -std

This fits a logistic regression model to the binary trait.  Only the
standardized SNP covariates are fitted, with an NEG prior assigned to
the regression coefficients.

The above command generates the following chatter output to the terminal:

Shrinkage parameters for SNP data:
lambda = 0.1 gamma = 0.018069 penalty = 50
Shrinkage parameters for other data:
lambda = 0 gamma = inf penalty = 0
Loading phen.dat
Dimensions of target file: 1000 1
Loading gen.dat
1000 100 0
Intercept: -0.619039
Opening output file test.R
Iteration: 0

The final line above indicates that only 1 iteration of algorithm is
used (default), so only one mode will be reported. The results viewed in R are:

R version 2.7.1 (2008-06-23)
Copyright (C) 2008 The R Foundation for Statistical Computing
ISBN 3-900051-07-0

.....

> dget("test.R")
[[1]]
     [,1]       [,2]       [,3]       [,4]       
[1,] "12"       "0"        "2"        "0"        
[2,] "-613.98"  "SNP1"       "SNP3"   "Intercept"
[3,] "-627.781" "0.471651" "0.330511" "-1.25968" 
[4,] "0"        "0"        "0"        "0"      

So HLasso has identified SNPs 1 and 3 as being associated with the
phenotype, and not SNPs 2 and 4; this is similar to the results of a
logistic regression which finds that SNPs 1 and 3 are both significant
at 10^3, whereas neither 2 nor 4 are significant at 10^-2.  Repeating
the HLasso command with -iter 100 does not identify any new modes. The
(correct) additive model was associated with the two SNPs identified,
and the parameter estimates were approximately 0.47 and 0.33.  The
negative value for intercept indicates fewer cases than controls (here
350 and 650, respectively).

Example 2:

./runHLasso -target data/phen.dat -covariates data/cov.dat -shape 0.1 -lambda 50 -o test.R -model_d 2 -lambda_d 0.1 -std -std_d -genotypes data\gen.dat

This fits a logistic regression model, but this time with the
non-genetic covariates as well as the SNP data, and all covariates are
standardised. The SNP regression coefficients are assigned
a NEG prior and the covariates are assigned a normal prior.

> dget("test.R")
[[1]]
     [,1]       [,2]       [,3]      [,4]       [,5]       [,6]     [,7]       
[1,] "11"       "0"        "1"       "2"        "0"        "2"      "0"        
[2,] "-585.01"  "cov1"     "cov2"    "cov3"     "SNP1"     "SNP3"   "Intercept"
[3,] "-589.028" "0.512079" "0.34897" "0.097446" "0.482773" "0.3614" "-1.29689" 
[4,] "0"        "0"        "0"       "0"        "0"        "0"      "0"      

So all three covariates are correctly found to be associated with the
phenotype, with decreasing effect sizes.  SNPs 1 and 3 are again found
to be associated with the phenotype and, perhaps because of the
residual variation explained by the covariates, their effect size
estimates are now slightly higher.

Example 3:

./runHLasso -target data/cont_phen.dat -covariates data/cov.dat -shape .1 -lambda 100 -o test.R -lambda_d 50 -shape 0.1 -genotypes data\gen.dat -linear

This fits a linear regression model to the continuous trait using
non-standardized SNP and covariate data. Both types of data are
assigned an NEG prior, but with different parameters.
