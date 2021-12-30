---
name: ExtremeValue
topic: Extreme Value Analysis
maintainer: Christophe Dutang, Kevin Jaunatre
email: Christophe.Dutang@ensimag.fr
version: 2021-12-29
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
and Achim Zeileis for his useful comments. If you think information is
not accurate or if we have omitted a package or important information
that should be mentioned here, please send an e-mail or submit an issue
or pull request in the GitHub repository linked above.

[Univariate Extreme Value Theory:]{#UnivEVT}
--------------------------------------------

-   ### Block Maxima approach:

    -   The package `r pkg("evd", priority = "core")`
        provides functions for a wide range of univariate distributions.
        Modelling function allow estimation of parameters for standard
        univariate extreme value methods.
    -   The package `r pkg("evdbayes")` provides the
        Bayesian analysis of univariate extreme value models using MCMC
        methods. It uses likelihood to estimate the parameters of the
        GEV distributions.
    -   The package `r pkg("revdbayes")` provides the
        Bayesian analysis of univariate extreme value models using
        direct random sampling from the posterior distribution, that is,
        without using MCMC methods.
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
    -   The package `r pkg("lmom")` has functions to fit
        probability distributions from GEV distributions to data using
        the low-order L-moments.
    -   The package `r pkg("lmomRFA")` extends package
        `r pkg("lmom")` and implements all the major
        components for regional frequency analysis using L-moments.
    -   The package `r pkg("ismev")` provides a collection
        of three functions to fit the GEV (diagnostic plot, MLE,
        likelihood profile) and follows the book of Coles (2001).
    -   The package `r pkg("mev")` has a function using the
        Smith (1987) penultimate approximation for block maxima
        approach.
    -   The package `r pkg("Renext")` provides various
        functions to fit the GEV distribution using an aggregated marked
        POT process.

-   ### Peak-Over-Threshold by GPD approach:

    -   The package `r pkg("evd")` includes univariate
        estimation for GPD approach by MLE.
    -   The Bayesian analysis of univariate extreme value models using
        MCMC methods in the package `r pkg("evdbayes")`
        includes the likelihood to estimate GP distributions.
    -   The package `r pkg("revdbayes")` provides the
        Bayesian analysis of univariate extreme value models using
        direct random sampling from the posterior distribution, that is,
        without using MCMC methods.
    -   The package `r pkg("evir")` performs modelling of
        univariate GPD by maximum likelihood fitting.
    -   The package `r pkg("evmix")` provides kernel density
        estimation and extreme value modelling. It also implements
        extreme value models and includes help on the choice of the
        threshold within those models using MLE.
    -   The package `r pkg("extremefit")` provides
        modelization of exceedances over a threshold in the Pareto type
        tail. It computes an adaptive choice of the threshold.
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
    -   The package `r pkg("lmom")` includes functions to
        fit probability distributions from GPD to data using the
        low-order L-moments.
    -   The package `r pkg("lmomRFA")` extends package
        `r pkg("lmom")` and implements all the major
        components for regional frequency analysis using L-moments.
    -   The package `r pkg("POT")` provides multiple
        estimators of the GPD parameters (MLE, L-Moments, method of
        median, minimum density power divergence). L-moments diagrams
        and from the properties of a non-homogeneous Poisson process
        techniques are provided for the selection of the threshold.
    -   The package `r pkg("ismev")` provides a collection
        of three functions to fit the GPD (diagnostic plot, MLE over a
        range of thresholds, likelihood profile) and follows the book of
        Coles (2OO1).
    -   The package `r pkg("mev")` provides functions to
        simulate data from GPD and multiple method to estimate the
        parameters (optimization, MLE, Bayesian methods and the method
        used in the `r pkg("ismev")` package).
    -   The package `r pkg("QRM")` provides functions to fit
        and graphically assess the fit of the GPD.
    -   The package `r pkg("Renext")` provides various
        functions to fit and assess the GPD distribution using an
        aggregated marked POT process.
    -   The package `r pkg("threshr")` deals with the
        selection of thresholds using a Bayesian leave-one-out
        cross-validation approach in order to compare the predictive
        performance resulting from a set of thresholds.

-   ### Extremal index estimation approach:

    -   The package `r pkg("evd")` implements univariate
        estimation for extremal index estimation approach.
    -   The package `r pkg("evdbayes")` includes point
        process characterisation
    -   the package `r pkg("evir")` includes extremal index
        estimation.
    -   The package `r pkg("extRemes")` also provides EVDs
        univariate estimation for the block maxima and poisson point
        process approache by MLE. It also incorporates a
        non-stationarity through the parameters.
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

-   ### Regression models

    -   The package `r pkg("VGAM")` offers additive
        modelling for extreme value analysis. The estimation for vector
        generalised additive models is performed using a backfitting
        algorithm and employs a penalized likelihood for the smoothing
        splines. It is the only package known to the authors that
        performs additive modelling for a range of extreme value
        analysis. It includes both GEV and GP distributions.
    -   The package `r pkg("ismev")` provides a collection
        of functions to fit a point process with explanatory variables
        (diagnostic plot, MLE) and follows the book of Coles (2001).

-   ### Copula approach:

    -   The package `r pkg("copula")` provides utilities for
        exploring and modelling a wide range of commonly used copulas,
        see also the `r view("Distributions")` task view
        (copula section).

-   ### Mixture distribution or composite distribution approach:

    -   The package `r pkg("evmix")` provides several
        functions to fit mixture distributions: either parametric / GPD,
        semi-parametric / GPD or non-parametric / GPD.

[Bivariate Extreme Value Theory:]{#BivEVT}
------------------------------------------

-   ### Block Maxima approach:

    -   The package `r pkg("evd")` provides functions for
        multivariate distributions. Modelling function allow estimation
        of parameters for class of bivariate extreme value
        distributions. Both parametric and non-parametric estimation of
        bivariate EVD can be performed.

-   ### Peak-Over-Threshold by GPD approach:

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

-   ### Tail dependence coefficient approach:

    -   The package `r pkg("RTDE")` implements bivariate
        estimation for the tail dependence coefficient.

[Multivariate Extreme Value Theory:]{#MultiEVT}
-----------------------------------------------

-   ### Block Maxima approach:

    -   The package `r pkg("lmomco")` is similar to the
        `r pkg("lmom")` but also implements recent advances
        in L-moments estimation, including L-moments for censored data,
        trimmed L-moments and L-moment for multivariate analysis for GEV
        distributions.

-   ### Peak-Over-Threshold by GPD approach:

    -   The package `r pkg("lmomco")` also implements
        L-moments multivariate analysis for GPD distributions.

-   ### Copula approach:

    -   The package `r pkg("copula")` provides utilities for
        exploring and modelling a wide range of commonly used copulas.
        Extreme value copulas and non-parametric estimates of extreme
        value copulas are implemented. See also the
        `r view("Distributions")` task view (copula
        section).

[Classical graphics:]{#Graphics}
--------------------------------

\

  --------------------------- ------------- ------------------------------------
  *Graphic name*              *Packages*    *Function names*
  Dispersion index plot       POT           `diplot`
  Distribution fitting plot   extremeStat   `distLplot`
  Hill plot                   evir          `hill`
  Hill plot                   evmix         `hillplot`
  Hill plot                   extremefit    `hill`
  Hill plot                   QRM           `hillPlot`
  Hill plot                   ReIns         `Hill`
  L-moment plot               POT           `lmomplot`
  Mean residual life plot     POT           `mrlplot`
  Mean residual life plot     evd           `mrlplot`
  Mean residual life plot     evir          `meplot`
  Mean residual life plot     evmix         `mrlplot`
  Mean residual life plot     ismev         `mrl.plot`
  Mean residual life plot     QRM           `MEplot`
  Mean residual life plot     ReIns         `MeanExcess`
  Pickand's plot              evmix         `pickandsplot`
  QQ Pareto plot              POT           `qplot`
  QQ Pareto plot              RTDE          `qqparetoplot`
  QQ Pareto plot              QRM           `plotFittedGPDvsEmpiricalExcesses`
  QQ Pareto plot              ReIns         `ParetoQQ`
  QQ Exponential plot         QRM           `QQplot`
  QQ Exponential plot         ReIns         `ExpQQ`
  QQ Exponential plot         Renext        `expplot`
  QQ Lognormal plot           ReIns         `LognormalQQ`
  QQ Weibull plot             ReIns         `WeibullQQ`
  QQ Weibull plot             Renext        `weibplot`
  Risk measure plot           QRM           `RMplot`
  Threshold choice plot       evd           `tcplot`
  Threshold choice plot       evmix         `tcplot`
  Threshold choice plot       POT           `tcplot`
  Threshold choice plot       QRM           `xiplot`
  Return level plot           POT           `retlev`
  Return level plot           POT           `Return`
  Return level plot           Renext        `plot,lines`
  --------------------------- ------------- ------------------------------------

  :  Graphics for univariate extreme value analysis

\
\

  ----------------------------------------------- ----- --------------
  Bivariate threshold choice plot                 evd   `bvtcplot`
  Dependence measure (chi) plot                   POT   `chimeas`
  Dependence measure (chi) plot                   evd   `chiplot`
  Dependence diagnostic plot within time series   POT   `tsdep.plot`
  Extremal index plot                             POT   `exiplot`
  Extremal index plot                             evd   `exiplot`
  Pickands' dependence function plot             POT   `pickdep`
  Spectral density plot                           POT   `specdens`
  ----------------------------------------------- ----- --------------

  :  Graphics for multivariate extreme value analysis

\

[Classical books and review papers:]{#References}
-------------------------------------------------

-   E. Gilleland, M. Ribatet, A. Stephenson (2013). A Software Review
    for Extreme Value Analysis, *Extremes* , **16** , 103-119.
-   R.-D. Reiss, M. Thomas (2007). *Statistical Analysis of Extreme
    Values with Applications to Insurance, Finance, Hydrology and Other
    Fields* , Springer-Verlag.
-   L. de Haan, A. Ferreira (2006). *Extreme Value Theory: An
    Introduction* , Springer-Verlag.
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



### Links
-   [Download statistics per view](http://rpackages.io/views)
-   [Gilleland, Eric, Mathieu Ribatet, and Alec G. Stephenson, A software review for extreme value analysis Extremes 16(1) (2013): 103-119](http://link.springer.com/article/10.1007/s10687-012-0155-0)