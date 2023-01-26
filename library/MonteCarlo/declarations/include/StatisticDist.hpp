/**
  * @file StatisticDist.hpp 
  *
  * @brief declarations of available distributions to be used 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef STATISTIC_DIST_DECLARATIONS 
#define STATISTIC_DIST_DECLARATIONS 


// Standard Normal Distribution 

namespace MonteCarlo {


    template < typename R >
    /**
      * Compute CDF of a given standard normal random variable 
      *
      * @tparam R a type of floating point e.g. double 
      *
      * @param x a standard normal random variables in R 
      * @return P( X < x ) for X ~ N( 0, 1 ) 
      */
    R StdNormCDF ( const R x );


    template < typename Z, typename R >
    /**
      * Compute inverse of standard normal CDF for a given quantile 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * 
      * @return z such that CDF(z) = quantile 
      */
    R InvStdNormCDF ( const R quantile );


} // MonteCarlo : Standard Normal Distribution 


// Log Normal Distribution 

namespace MonteCarlo {


    template < typename R >
    /**
      * Compute inverse of log normal CDF for a given quantile 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      *
      * @param mu mean of distribution 
      * @param sig standard deviation of distribution 
      * 
      * @return z such that CDF(z) = quantile 
      */
    R InvLogNormCDF ( const R quantile, const R mu, const R sig );


} // MonteCarlo : Log Normal Distribution 


#ifndef STATISTIC_DIST_IMPLEMENTATIONS 
    #include "StatisticDist_imp.hpp" 
#endif 

#endif // STATISTIC_DIST_DECLARATIONS 

