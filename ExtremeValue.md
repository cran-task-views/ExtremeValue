---
name: ExtremeValue
topic: Extreme Value Analysis
maintainer: Christophe Dutang
email: dutangc@gmail.com
version: 2023-03-30
source: https://github.com/cran-task-views/ExtremeValue/
---


Extreme values modelling and estimation are an important challenge in
various domains of application, such as environment, hydrology, finance,
actuarial science, just to name a few. The restriction to the analysis
of extreme values may be justified since the extreme part of a sample
can be of a great importance. That is, it may exhibit a larger risk
potential such as high concentration of air pollutants, flood, extreme
claim sizes, price shocks in the four previous topics respectively. The
statistical analysis of extreme may be spread out in many packages
depending on the topic of application. In this task view, we present the
packages from a methodological side.

Applications of extreme value theory can be found in other task views:
for financial and actuarial analysis in the
`r view("Finance")` task view, for environmental analysis in
the `r view("Environmetrics")` task view. General
implementation of probability distributions is studied in the
`r view("Distributions")` task view.

The maintainers gratefully acknowledge E. Gilleland, M. Ribatet and A.
Stephenson for their review for extreme value analysis packages (2013)
Kevin Jaunatre for his helpful advice
and Achim Zeileis for his useful comments. If you think information is
not accurate or if we have omitted a package or important information
that should be mentioned here, please send an e-mail or submit an issue
or pull request in the GitHub repository linked above.


### Table of contents
- [Univariate Extreme Value Theory](#UnivEVT)
  - [Bayesian approach](#UniBayesian)
  - [Block Maxima approach](#UniBlockMaxima)
  - [Extremal index estimation approach](#UniExtremalIndex)
  - [Mixture distribution or composite distribution approach](#UniMixture)
  - [Peak-Over-Threshold by GPD approach](#UniPOT)
  - [Record models](#UniRecord)
  - [Regression models](#UniRegression)
  - [Threshold selection](#UniThreshold)
  
- [Bivariate Extreme Value Theory](#BivEVT) 
  - [Copula approach](#BiCopula)
  - [Maxima approach](#BiMaxima)
  - [Peak-Over-Threshold by GPD approach](#BiPOT)
  - [Tail dependence coefficient approach](#BiTailDependence)
  
- [Multivariate Extreme Value Theory](#MultiEVT) 
  - [Bayesian approach](#MultiBayesian)
  - [Copula approach](#MultiCopula)
  - [Multivariate Maxima](#MultiMaxima)
  - [Peak-Over-Threshold by GPD approach](#BiPOT)
  - [Tail dependence coefficient approach](#MultiTailDependence)        
  - [Statistical tests](#MultiTests)

- [Classical graphics](#Graphics)
  
  

[Univariate Extreme Value Theory:]{#UnivEVT}
--------------------------------------------

Several packages export the probability functions (quantile, density,
distribution and random generation) for the Generalized Pareto and the 
Generalized Extreme Value distributions, often sticking to the 
classical prefixing rule (with prefixes `"q"`, `"d"`, `"p"`, `"r"`) 
and allowing the use of the formals such as `log` and `lower tail`, 
see the [Distributions](https://cran.r-project.org/view=Distributions) 
task view for details. 
Several strategies can be used for the numeric evaluation of these 
functions in the small shape (near exponential) case.  
Also, some implementations allow the use of parameters in 
vectorized form and some can provide the derivatives w.r.t. the parameters.
Nevertheless, the `r pkg("nieve")` package provides symbolic
differentiation for two EVT probability distribution (GPD and GEV) 
in order to compute the log-likelihood. 

-   ### [Bayesian approach:]{#UniBayesian}

    -   The package `r pkg("extRemes")` also provides bayesian estimation.    
    -   The package `r pkg("MCMC4Extremes")` proposes some functions 
        to perform posterior estimation for some distribution, with an
        emphasis to extreme value distributions.
    -   The package `r pkg("revdbayes")` provides the
        Bayesian analysis of univariate extreme value models using
        direct random sampling from the posterior distribution, that is,
        without using MCMC methods.    
    -   The package `r pkg("texmex")` fit GPD models by using maximum 
        (optionally penalised-)likelihood, or Bayesian estimation, and 
        both classes of models may be fitted with covariates in any/all model parameters.         
    
|package        | function      | models[^1]  | covariates  | sampling[^2]  | prior choice  | generic functions |
|:--------------|:--------------|:------------|:------------|:--------------|:--------------|:------------------|
|`extRemes`     | `fevd`        | 1--4,*      | all         | RWMH          | custom        | plot, summary       |
|`MCMC4Extremes`|`ggev`,`gpdp`  | 1--2,*      | no          | RWMH          | fixed         | plot, summary       |
|`revdbayes`    | `rpost`       | 1--4        | no          | RU            | custom        | plot, summary       |
|`texmex`       | `evm`         | 1--2,*      | all         | IMH           | gaussian      | plot, summary, density,correlogram|
        
[^1] model family: generalized extreme value distribution (1), generalized Pareto distribution (2), inhomogeneous Poisson process (3), order statistics/r-largest (4) or custom/other (*).    

[^2] sampling: random walk Metropolis--Hastings (RWMH), exact sampling ratio-of-uniform (RU), independent Metropolis--Hastings (IMH) 

-   ### [Block Maxima approach:]{#UniBlockMaxima}

    -   The package `r pkg("climextRemes")` provides functions for 
        fitting GEV via point process fitting for extremes in climate data,
        providing return values, return probabilities, and return 
        periods for stationary and nonstationary models. 
    -   The package `r pkg("evd", priority = "core")`
        provides functions for a wide range of univariate distributions.
        Modelling function allow estimation of parameters for standard
        univariate extreme value methods.
    -   The package `r pkg("evir", priority = "core")`
        performs modelling of univariate GEV distributions by maximum
        likelihood fitting.
    -   The package `r pkg("extRemes")` provides EVDs
        univariate estimation for block maxima model approache by MLE.
        It also incorporates a non-stationarity through the parameters
        of the EVDs and L-moments estimation for the stationary case for
        the GEV distributions. Finally, it has also Bayes estimation
        capabilities. A separate package
        `r pkg("in2extRemes")` provides some GUI interfaces
        to `r pkg("extRemes")`.
    -   The package `r pkg("extremeStat")` includes
        functions to fit multiple GEV distributions types available in
        the package `r pkg("lmomco")` using linear moments
        to estimate the parameters.
    -   The package `r pkg("fExtremes")` provides univariate
        data processing and modelling. It includes clustering, block
        maxima identification and exploratory analysis. The estimation
        of stationary models for the GEV is provided by maximum
        likelihood and probability weighted moments.
    -   The package `r pkg("ismev")` provides a collection
        of three functions to fit the GEV (diagnostic plot, MLE,
        likelihood profile) and follows the book of Coles (2001).    
    -   The package `r pkg("lmom")` has functions to fit
        probability distributions from GEV distributions to data using
        the low-order L-moments.
    -   The package `r pkg("lmomRFA")` extends package
        `r pkg("lmom")` and implements all the major
        components for regional frequency analysis using L-moments.
    -   The package `r pkg("QRM")` provides a function to fit GEV    
        in Quantitative Risk Management perspective.
    -   The package `r pkg("Renext")` provides various
        functions to fit the GEV distribution using an aggregated marked
        POT process.


Summary of GEV density functions and GEV fitting functions

|  package	      |density function | location	  |scale	    |shape	    |fit function	      |argdata    |outputS4	  |outputS3       |outputS3par   |
|-----------------|:--------|:------------|:----------|:----------|:----------|:---------|:-----------|:--------------|:-------------|
| climextRemes	  |NA       |`location`	  |`scale`	  |`shape`	  |`fit_gev`	|`y`	     |NA	        |`mle`	        |NA             |
| evd	            |`dgev`	    |`loc`	      |`scale`	  |`shape`	  |`fgev`	    |`x`	     |NA	        |`estimate`	    |NA             |
| evir	          |`dgev`	    |`mu`	        |`sigma`	  |`xi`	      |`gev`	    |`data`	   |NA	        |`par.ests`	    |NA             |
| extraDistr	    |`dgev`	    |`mu`	        |`sigma`	  |`xi`	      |NA	        |NA	       |NA	        |NA	            |NA             |  
| extRemes	      |`devd`	    |`loc`	      |`scale`	  |`shape`	  |`fevd`	    |`x`	     |NA	        |`results`	    |`par`          |
| fExtremes	      |`dgev`	    |`mu`	        |`beta`	    |`xi`       |`gevFit`	  |`x` 	     |`fit` 	    |`par.ests`	    |NA             |
| ismev	          |NA	      |NA	          |NA	        |NA	        |`gev.fit`	|`xdat`	   |NA	        |`mle`	        |NA             |
| lmomco	        |`pdfgev`	    |`xi`	        |`alpha`	  |`kappa`	  |NA	        |NA	        |NA	        |NA	            |NA             |
| QRM	            |`dGEV`	    |`mu`	        |`sigma`	  |`xi`	      |`fit.GEV`	|`maxima`	  |NA	        |`par.ests`	    |NA             |
| revdbayes	      |`dgev`	    |`loc`	      |`scale`	  |`shape`	  |NA	        |NA	        |NA	        |NA	            |NA             |
| SpatialExtremes	|`dgev`	    |`loc`	      |`scale`	  |`shape`	  |NA	        |NA	        |NA	        |NA	            |NA             |
| texmex	        |`dgev`	    |`mu`	        |`sigma`	  |`xi`  	    |`evm`	    |`y`	     |NA	        |`coefficients`	|NA             |
| TLMoments	      |`dgev`	    |`loc`	      |`scale`	  |`shape`	  |NA	        |NA	        |NA	        |NA	            |NA             |


-   ### [Extremal index estimation approach:]{#UniExtremalIndex}

    -   The package `r pkg("evd")` implements univariate
        estimation for extremal index estimation approach.
    -   the package `r pkg("evir")` includes extremal index
        estimation.
    -   The package `r pkg("extRemes")` also provides EVDs
        univariate estimation for the block maxima and poisson point
        process approache by MLE. It also incorporates a
        non-stationarity through the parameters.
    -   The package `r pkg("extremefit")` provides
        modelization of exceedances over a threshold in the Pareto type
        tail. It computes an adaptive choice of the threshold.  
    -   The package `r pkg("ExtremeRisks")` provides risk measures
        such as Expectile, Value-at-Risk, for univariate independent 
        observations and temporal dependent observations. 
        The statistical inference is performed through parametric 
        and non-parametric estimators. 
        Inferential procedures such as confidence intervals, confidence 
        regions and hypothesis testing are obtained by exploiting the asymptotic theory.     
    -   The package `r pkg("fExtremes")` provides univariate
        data processing and modelling. It includes extremal index
        estimation.
    -   The package `r pkg("mev")` provides extremal index
        estimators based on interexceedance time (MLE and iteratively
        reweigthed least square estimators of Suveges (2007)). It
        provides the information matrix test statistic proposed by
        Suveges and Davison (2010) and MLE for the extremal index.
    -   The package `r pkg("ReIns")` provides functions for
        extremal index and splicing approaches in a reinsurance
        perspective.
    -   The package `r pkg("ptsuite")` implements various
        estimation methods for the shape parameter of Pareto distributed
        data.
    -   The package `r pkg("evgam")` implements a moment-based estimator
        of extremal index based on Ferro and Segers (2003).


-   ### [Mixture distribution or composite distribution approach:]{#UniMixture}

    -   The package `r pkg("evmix")` provides kernel density
        estimation and extreme value modelling. It also implements
        mixture extreme value models and includes help on the choice of the
        threshold within those models using MLE: either parametric / GPD,
        semi-parametric / GPD or non-parametric / GPD.
        

-   ### [Peak-Over-Threshold by GPD approach:]{#UniPOT}

    -   The package `r pkg("ercv")` provides a methodology to fit a 
        generalized Pareto distribution, together with an automatic 
        threshold selection algorithm.
    -   The package `r pkg("eva")` provides Goodness-of-fit tests 
        for selection of r in the r-largest order statistics
        and threshold selection.
    -   The package `r pkg("evd")` includes univariate
        estimation for GPD approach by MLE.
    -   The package `r pkg("evir")` performs modelling of
        univariate GPD by maximum likelihood fitting.
    -   The package `r pkg("extRemes")` provides EVDs
        univariate estimation for GPD approach by MLE. A
        non-stationarity through the parameters of the EVDs and
        L-moments estimation for the stationnary case for the GPD
        distributions is also included.
    -   The package `r pkg("extremeStat")` includes
        functions to fit multiple GPD distributions types available in
        the package `r pkg("lmomco")` using linear moments
        to estimate the parameters.
    -   The package `r pkg("fExtremes")` includes the
        estimation of stationary models for the GPD by maximum
        likelihood and probability weighted moments.
    -   The package `r pkg("ismev")` provides a collection
        of three functions to fit the GPD (diagnostic plot, MLE over a
        range of thresholds, likelihood profile) and follows the book of
        Coles (2OO1).    
    -   The package `r pkg("lmom")` includes functions to
        fit probability distributions from GPD to data using the
        low-order L-moments.
    -   The package `r pkg("lmomRFA")` extends package
        `r pkg("lmom")` and implements all the major
        components for regional frequency analysis using L-moments.
    -   The package `r pkg("mev")` provides functions to
        simulate data from GPD and multiple method to estimate the
        parameters (optimization, MLE, Bayesian methods and the method
        used in the `r pkg("ismev")` package).        
    -   The package `r pkg("POT")` provides multiple
        estimators of the GPD parameters (MLE, L-Moments, method of
        median, minimum density power divergence). L-moments diagrams
        and from the properties of a non-homogeneous Poisson process
        techniques are provided for the selection of the threshold.
    -   The package `r pkg("QRM")` provides functions to fit
        and graphically assess the fit of the GPD.
    -   The package `r pkg("ReIns")` provides a
        function to fit the GPD distribution as well as the extended
        Pareto distribution.    
    -   The package `r pkg("Renext")` provides various
        functions to fit and assess the GPD distribution using an
        aggregated marked POT process.
    -   The package `r pkg("SpatialExtremes")` provides a
        function to fit the GPD distribution.       
    -   The package `r pkg("SpatialExtremes")` provides different 
        approaches for fitting/selecting the threshold in generalized 
        Pareto distributions. Most of them are based on minimizing the 
        AMSE-criterion or at least by reducing the bias of the assumed GPD-model.
    -   The package `r pkg("texmex")` fit GPD models by using maximum 
        (optionally penalised-)likelihood, or Bayesian estimation, and 
        both classes of models may be fitted with covariates in any/all model parameters.     
    -   The package `r pkg("NHPoisson")` provides a function to fit
        non-homogeneous Poisson processes for peak over threshold analysis. 
    
    
Summary of GPD density functions and GPD fitting functions

|package	        |density function | location  |scale|shape	  |fit function |argdata    | argthres    |outputS4	  |outputS3       |outputS3par   |
|-----------------|:--------|:---------|:----------|:-----------|:------------|:----------|:------------|:----------|:--------------|:--------------|
| ercv	          |NA	      |NA	        |NA	        |NA	        |`fitpot`	    |`data`	    |`threshold`	|NA	        |`coeff`	      |NA           |
| eva	            |`dgpd`	  |`loc`	    |`scale`	  |`shape`	  |`gpdFit`	    |`data`	    |`threshold`	|NA	        |`par.ests`	    |NA           |
| evd	           | `dgpd`	  |`loc`	    |`scale`	  |`shape`	  |`fpot`	      |`x`	      |`threshold`	|NA	        |`estimate`	    |NA           |
| evir	          |`dgpd`	  |`mu`	      |`beta`	    |`xi`	      |`gpd`	      |`data`	    |`threshold`	|NA	        |`par.ests`	    |NA           |
| extraDistr	    |`dgpd`	  |`mu`	      |`sigma`	  |`xi`	      |NA	          |NA	        |NA	          |NA	        |NA	            |NA           |
| extRemes	      |`devd` 	|`loc`	    |`scale`	  |`shape`	  |`fevd`	      |`x`	      |`threshold`	|NA	        |`results`	    |`par`           |
| fExtremes	      |`dgpd`	  |`mu`  	    |`beta`	    |`xi`	      |`gpdFit`	    |`x`	      |`u`	        |`fit`	    |`fit`	        |`par`           |
| ismev	          |NA	      |NA	        |NA	        |NA	        |`gpd.fit`	  |`xdat`	    |`threshold`	|NA	        |`mle`	         |NA           |
| lmomco	        |`pdfgpa`	|`xi`	      |`alpha`	  |`kappa`	  |NA	          |NA	        |NA	          |NA	        |NA	            |NA           |
| mev	            |NA	      |NA	        |`scale`	  |`shape`	  |`fit.gpd`	  |`xdat`	    |`threshold`	|NA	        |`estimate`	    |NA           |
| POT	            |`dgpd`	  |`loc`	    |`scale`	  |`shape`	  |`fitgpd`	    |`data`	    |`threshold`	|NA	        |`fitted.values`|NA           |
| QRM	            |`dGPD`	  |NA	        |`beta`	    |`xi`	      |`fit.GPD`	  |`data`	    |`threshold`	|NA	        |`par.ests`	    |NA           |
| ReIns	          |`dgpd`	  |`mu`  	    |`sigma`	  |`gamma`	  |`GPDfit`	    |`data`	    |NA	          |NA	        |NA	            |NA           |
| Renext	        |`dGPD`	  |`loc`	    |`scale`	  |`shape`	  |`fGPD`	      |`x`	      |NA	          |NA	        |`estimate`	    |NA           |
| revdbayes	      |`dgp`	  |`loc`	    |`scale`	  |`shape`	   | NA	        |NA	        |NA	          |NA	        |NA	            |NA           |
| SpatialExtremes	|`dgpd`	  |`loc`	    |`scale`	  |`shape`	  |`gpdmle`	    |`x`	      |`threshold`	|NA	        |NA	            |NA           |
| tea	            |`dgpd`	  |`loc`	    |`scale`	  |`shape`	  |`gpdFit`	    |`data`	    |`threshold`	|NA	        |`par.ests`	    |NA           |
| texmex	        |`dgpd`	  |`u`	      |`sigma`	  |`xi`	      |`evm`	      |`y`	       |`th`	      |NA	        |`coefficients`  |NA           |
| TLMoments	      |`dgpd`	  |`loc`	    |`scale`	  |`shape`	   | NA	        |NA	        |NA	          |NA	        |NA	            |NA           |
 

-   ### [Record models:]{#UniRecord}
    -   `r pkg("RecordTest")` studies the analysis of record-breaking events 
    and provides non-parametric modeling/testing of a non-stationary behaviour 
    in (extreme) records. 
    -   `r pkg("evir")` provides only a function `records()` for extracting records.


-   ### [Regression models:]{#UniRegression}

    -   The package `r pkg("VGAM")` offers additive
        modelling for extreme value analysis. The estimation for vector
        generalised additive models (GAM) is performed using a backfitting
        algorithm and employs a penalized likelihood for the smoothing
        splines. It is the only package known to the authors that
        performs additive modelling for a range of extreme value
        analysis. It includes both GEV and GP distributions.
    -   The package `r pkg("ismev")` provides a collection
        of functions to fit a point process with explanatory variables
        (diagnostic plot, MLE) and follows the book of Coles (2001).
    -   The package `r pkg("texmex")` fit GPD models by using maximum 
        (optionally penalised-)likelihood, or Bayesian estimation, and 
        both classes of models may be fitted with covariates in any/all model parameters.    
    -   The package `r pkg("evgam")` provides methods for fitting various 
        extreme value distributions with parameters of generalised additive 
        model (GAM) form. 

    
-   ### [Threshold selection:]{#UniThreshold}

    -   The package `r pkg("threshr")` deals with the
        selection of thresholds using a Bayesian leave-one-out
        cross-validation approach in order to compare the predictive
        performance resulting from a set of thresholds.
    -   The package `r pkg("ercv")` provides a methodology to fit a 
        generalized Pareto distribution, together with an automatic 
        threshold selection algorithm.
     -   The package `r pkg("POT")` provides multiple
        estimators of the GPD parameters (MLE, L-Moments, method of
        median, minimum density power divergence). L-moments diagrams
        and from the properties of a non-homogeneous Poisson process
        techniques are provided for the selection of the threshold.
      
        
[Bivariate Extreme Value Theory:]{#BivEVT}
------------------------------------------
        
-   ### [Copula approach:]{#BiCopula}

    -   The package `r pkg("copula")` provides utilities for
        exploring and modelling a wide range of commonly used copulas,
        see also the `r view("Distributions")` task view
        (copula section).        
    -   The pacage `r pkg("fCopulae")` provides utilities to fit
        bivariate extreme copulas.

-   ### [Maxima approach:]{#BiMaxima}

    -   The package `r pkg("evd")` provides functions for
        multivariate distributions. Modelling function allow estimation
        of parameters for class of bivariate extreme value
        distributions. Both parametric and non-parametric estimation of
        bivariate EVD can be performed.
    -   Nonparametric estimation of the spectral measure 
        using a sample of pseudo-angles is available in the package 
        `r pkg("extremis")` in the bivariate setting.    

-   ### [Peak-Over-Threshold by GPD approach:]{#BiPOT}

    -   The package `r pkg("evd")` implements bivariate
        threshold modelling using censored likelihood methodology.
    -   The single multivariate implementation in the package
        `r pkg("evir")` is a bivariate threshold method.
    -   The package `r pkg("extremefit")` provides
        modelization of exceedances over a threshold in the Pareto type
        tail depending on a time covariate. It provides an adaptive
        choice of the threshold depending of the covariate.
    -   The package `r pkg("POT")` provides estimators of
        the GPD parameters in the bivariate case.
    
-   ### [Tail dependence coefficient approach:]{#BiTailDependence}

    -   The package `r pkg("RTDE")` implements bivariate
        estimation for the tail dependence coefficient.

[Multivariate Extreme Value Theory:]{#MultiEVT}
-----------------------------------------------

-   ### [Bayesian approach:]{#MultiBayesian}

    -   The package `r pkg("SpatialExtremes")` provides 
        tools for the statistical modelling of spatial extremes 
        using Bayesian hierarchical models (fitting, checking,
        selection).
    -   The package `r pkg("ExtremalDep")` also provides function to fit 
        a multivariate extreme value using Bayesian inference.
        
-   ### [Copula approach:]{#MultiCopula}

    -   The package `r pkg("SpatialExtremes")` provides functions
        to estimate a copula-based model to spatial extremes as well
        as model checking and selection.
    -   The package `r pkg("copula")` provides utilities for
        exploring and modelling a wide range of commonly used copulas.
        Extreme value copulas and non-parametric estimates of extreme
        value copulas are implemented. See also the
        `r view("Distributions")` task view (copula
        section).
    -   The package `r pkg("SimCop")` has functionalities 
        for simulation of some bivariate extreme value 
        distributions and the multivariate logistic model, 
        or Gumbel copula.

-   ### [Multivariate Maxima:]{#MultiMaxima}

    -   The package `r pkg("lmomco")` is similar to the
        `r pkg("lmom")` but also implements recent advances
        in L-moments estimation, including L-moments for censored data,
        trimmed L-moments and L-moment for multivariate analysis for GEV
        distributions.
    -   The package `r pkg("SpatialExtremes")` provides functions
        to fit max-stable processes to data using pairwise likelihood
        or spatial GEV models possibly with covariates.
    -   A set of procedures for modelling parametrically 
        and non-parametrically the dependence structure of 
        multivariate extreme-values is provided in `r pkg("ExtremalDep")`.
    -   The `r pkg("BMAmevt")` package implements a Bayesian
        nonparametric model that uses a trans-dimensional Metropolis algorithm 
        for fitting a Dirichlet mixture to the spectral measure 
        based on pseudo-angles.

-   ### [Peak-Over-Threshold by GPD approach:]{#MultiPOT}

    -   The package `r pkg("lmomco")` also implements
        L-moments multivariate analysis for GPD distributions.
    -   The package `r pkg("graphicalExtremes")` develops a     
        statistical methodology for sparse multivariate extreme 
        value models. Methods are provided for exact simulation 
        and statistical inference for multivariate Pareto 
        distributions on graphical structures.
    
        
-   ### [Tail dependence coefficient approach:]{#MultiTailDependence}  

    -   The package `r pkg("SpatialExtremes")` provides functions
        to estimate non parametrically the extremal coefficient function
         as well as model checking and selection.
    -   The package `r pkg("ExtremeRisks")` provides risk measures
        such as Expectile, Value-at-Risk, for multivariate independent 
        marginals.
    -   The package `r pkg("tailDepFun")` provides functions 
        implementing minimal distance estimation methods for 
        parametric tail dependence models.

    
-   ### [Statistical tests:]{#MultiTests}

    -   The `r pkg("copula")` package includes three tests 
        of max-stability assumption.
        

[Classical graphics:]{#Graphics}
--------------------------------


Graphics for univariate extreme value analysis
  
|Graphic name               |Packages       |Function names                       |
|---------------------------|:--------------|:------------------------------------|
|Dispersion index plot       |POT           |`diplot`                       |
|Distribution fitting plot   |extremeStat   |`distLplot`                       |
|Hill plot                   |evir          |`hill`                         |
|Hill plot                   |evmix         |`hillplot`                       |
|Hill plot                   |extremefit    |`hill`                         |
|Hill plot                   |QRM           |`hillPlot`                       |
|Hill plot                   |ReIns         |`Hill`                         |
|Hill plot                   |ExtremeRisks  |`HTailIndex`                         |
|L-moment plot               |POT           |`lmomplot`                       |
|Mean residual life plot     |POT           |`mrlplot`                       |
|Mean residual life plot     |evd           |`mrlplot`                       |
|Mean residual life plot     |evir          |`meplot`                       |
|Mean residual life plot     |evmix         |`mrlplot`                       |
|Mean residual life plot     |ismev         |`mrl.plot`                       |
|Mean residual life plot     |QRM           |`MEplot`                       |
|Mean residual life plot     |ReIns         |`MeanExcess`                     |
|Pickand's plot              |evmix         |`pickandsplot`                   |
|QQ Pareto plot              |POT           |`qplot`                        |
|QQ Pareto plot              |RTDE          |`qqparetoplot`                   |
|QQ Pareto plot              |QRM           |`plotFittedGPDvsEmpiricalExcesses`   |
|QQ Pareto plot              |ReIns         |`ParetoQQ`                       |
|QQ Exponential plot         |QRM           |`QQplot`                       |
|QQ Exponential plot         |ReIns         |`ExpQQ`                        |
|QQ Exponential plot         |Renext        |`expplot`                       |
|QQ Lognormal plot           |ReIns         |`LognormalQQ`                   |
|QQ Weibull plot             |ReIns         |`WeibullQQ`                     |
|QQ Weibull plot             |Renext        |`weibplot`                     |
|Risk measure plot           |QRM           |`RMplot`                       |
|Threshold choice plot       |evd           |`tcplot`                       |
|Threshold choice plot       |evmix         |`tcplot`                       |
|Threshold choice plot       |POT           |`tcplot`                       |
|Threshold choice plot       |QRM           |`xiplot`                       |
|Return level plot           |POT           |`retlev`                       |
|Return level plot           |POT           |`Return`                       |
|Return level plot           |Renext        |`plot,lines`                       |


  



Graphics for multivariate extreme value analysis

| Graphic                                       | Package | Function     |
|:----------------------------------------------|:--------|:-------------|
| Angular densities plot                        | `ExtremalDep`   | `AngDensPlot`   |
| Bivariate threshold choice plot               | `evd`   | `bvtcplot`   |
| Dependence measure (chi) plot                 | `POT`   | `chimeas`    |
| Dependence measure (chi) plot                 | `evd`   | `chiplot`    |
| Dependence diagnostic plot within time series | `POT`   | `tsdep.plot` |
| Extremal index plot                           | `POT`   | `exiplot`    |
| Extremal index plot                           | `evd`   | `exiplot`    |
| (2D)map for a max-stable process              | `SpatialExtremes`|`map`|
| madogram for a max-stable process             | `SpatialExtremes`|`madogram`|
| madogram for a max-stable process             | `ExtremalDep`|`madogram`|
| F-madogram for a max-stable process           | `SpatialExtremes`|`fmadogram`|
| lambda-madogram for a max-stable process      | `SpatialExtremes`|`lmadogram`|
| Multidimensional Hill plot                    | `ExtremeRisks`|`MultiHTailIndex`|
| Pickands' dependence function plot            | `POT`   | `pickdep`    |
| Pickands' dependence function plot            | `ExtremalDep`   | `bbeed`    |
| QQ-plot for the extremal coefficient          | `SpatialExtremes`|`qqextcoeff`|
| Spectral density plot                         | `POT`   | `specdens`   |







[Classical books and review papers:]{#References}
-------------------------------------------------

-   E. Gilleland, M. Ribatet, A. Stephenson (2013). A Software Review
    for Extreme Value Analysis, *Extremes* , **16** , 103-119.
-   R.-D. Reiss, M. Thomas (2007). *Statistical Analysis of Extreme
    Values with Applications to Insurance, Finance, Hydrology and Other
    Fields* , Springer-Verlag.
-   L. de Haan, A. Ferreira (2006). *Extreme Value Theory: An
    Introduction* , Springer-Verlag.
-   Stephenson AG, Gilleland E (2006). Software for the analysis of extreme events: The current state and future directions. *Extremes*, **8**, 87–109.
-   J. Beirlant, Y. Goegebeur, J. Teugels, J. Segers (2004). *Statistics
    of Extremes: Theory and Applications* , John Wiley & Sons.
-   B. Finkenstaedt, H. Rootzen (2004). *Extreme Values in Finance,
    Telecommunications, and the Environment* , Chapman & Hall/CRC.
-   S. Coles (2001). *An Introduction to Statistical Modeling of Extreme
    Values* , Springer-Verlag.
-   P. Embrechts, C. Klueppelberg, T. Mikosch (1997). *Modelling
    Extremal Events for Insurance and Finance* , Springer-Verlag.
-   S.I. Resnick (1987). *Extreme Values, Regular Variation and Point
    Processes* , Springer-Verlag.
-   Smith, R.L. (1987). Approximations in extreme value theory.
    Technical report 205, Center for Stochastic Process, University of
    North Carolina, 1--34.
-   Suveges (2007) Likelihood estimation of the extremal index.
    Extremes, 10(1), 41-55.
-   Suveges and Davison (2010), Model misspecification in peaks over
    threshold analysis. Annals of Applied Statistics, 4(1), 203-221.



### Links to review papers

-   [Gilleland, Eric, Mathieu Ribatet, and Alec G. Stephenson, A software review for extreme value analysis Extremes 16(1) (2013): 103-119](https://doi.org/10.1007/s10687-012-0155-0)
-   [Alec G. Stephenson and Eric Gilleland, Software for the analysis of extreme events: The current state and future directions. Extremes 8:87–109 (2006)](https://doi.org/10.1007/s10687-006-7962-0)
