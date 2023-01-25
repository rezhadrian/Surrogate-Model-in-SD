/**
  * @file SupplementaryMaths_MC.hpp
  *
  * @brief declarations of math functions not implemented in cmath 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef SUPPLEMENTARY_MATHS_MC_DECLARATIONS 
#define SUPPLEMENTARY_MATHS_MC_DECLARATIONS 


namespace MonteCarlo {


    template < typename Z, typename R >
    /**
      * Raise number to given power recursively 
      *
      * @tparam Z a type of non-negative integers e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      */
    R Power ( const R number, const Z power );


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


    template < typename R >
    /**
      * Compute inverse of complement error function 
      *
      * @tparam R a type of floating point e.g. double 
      * @return x such that erfc(z) = x  
      */
    R InvErfc ( const R p );


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


} // MonteCarlo : SupplementaryMaths 


#endif // SUPPLEMENTARY_MATHS_MC_DECLARATIONS 

