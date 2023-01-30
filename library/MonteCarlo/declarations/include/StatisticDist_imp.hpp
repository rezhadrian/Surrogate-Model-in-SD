/**
  * @file StatisticDist_imp.hpp 
  *
  * @brief implementations  of available distributions to be used 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef STATISTIC_DIST_IMPLEMENTATIONS 
#define STATISTIC_DIST_IMPLEMENTATIONS 

#ifndef STATISTIC_DIST_DECLARATIONS 
    #include "StatisticDist.hpp" 
#endif 

#include "LibrariesLoader_MC.hpp" 


namespace MonteCarlo {


    template < typename Z, typename R >
    /**
      * Raise number to given power recursively 
      *
      * @tparam Z a type of non-negative integers e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      */
    R Power ( const R number, const Z power );


    template < typename R >
    /**
      * Compute inverse of complement error function 
      *
      * @tparam R a type of floating point e.g. double 
      * @return x such that erfc(z) = x  
      */
    R InvErfc ( const R p );


} // MonteCarlo : Additional math functions not in cmath 


namespace MonteCarlo {

    template < typename R >
    R StdNormCDF ( const R x ) {

        return 0.5 * std::erfc ( -x * std::sqrt(0.5) );

    }

} // MonteCarlo : InvStdNormCDF 


namespace MonteCarlo {

    template < typename Z, typename R >
    /**
      * The algorithm used is the Beasley-Springer-Moro algorithm 
      * 
      * Reference:
      *
      * Glasserman, P. (2004). Monte Carlo methods in financial engineering. 
      * New York, NY : Springer.
      */
    R InvStdNormCDF ( const R quantile ) {

        if ( quantile < 0.0 ) {

            throw std::runtime_error ( 
                "InverseSingleCDF: quantile must be non-negative" 
            );

        }

        if ( quantile >= 1.0 ) {

            throw std::runtime_error ( 
                "InverseSingleCDF: quantile must be below unity" 
            );

        }

        if ( quantile == 1.0 ) {

            return 1.0;

        }

        const R a[4] = {

              2.50662823884,
            -18.61500062529,
             41.39119773534,
            -25.44106049637

        };

        const R b[4] = {

             -8.47351093090,
             23.08336743743,
            -21.06224101826,
              3.13082909833

        };

        const R c[9] = {

            0.3374754822726147,
            0.9761690190917186,
            0.1607979714918209,
            0.0276438810333863,
            0.0038405729373609,
            0.0003951896511919,
            0.0000321767881768,
            0.0000002888167364,
            0.0000003960315187

        };

        if ( quantile < 0.5 ) {
            return -1.0 * InvStdNormCDF<Z,R> ( 1 - quantile );
        }

        if ( quantile >= 0.5 && quantile <= 0.92 ) {

            R   numerator = 0.0;
            R denominator = 1.0;

            for ( auto i = 0; i < 4; i++ ) {
                  numerator += a[i] * Power<Z,R> ( quantile - 0.5, 2*i + 1 ); 
                denominator += b[i] * Power<Z,R> ( quantile - 0.5, 2*i + 2 ); 
            }

            return numerator / denominator;

        }

        R result = 0.0;

        for ( auto i = 0; i < 9; i++ ) {
            result += c[i] * Power<Z,R> (
                std::log ( - std::log (1-quantile) ),
                i
            );
        }

        return result;

    }

} // MonteCarlo : InvStdNormCDF 


namespace MonteCarlo {

    template < typename R >
    R InvLogNormCDF ( const R quantile, const R mu, const R sig ) {

        return std::exp (-1.41421356237309505*sig*InvErfc(2.0*quantile)+mu);

    }

} // MonteCarlo : InvLogNormCDF 


namespace MonteCarlo {

    template < typename Z, typename R >
    R Power ( const R number, const Z power ) {

        R result = 1.0;

        for ( auto i = 0; i < power; i++ ) {
            result *= number;
        }

        return result;

    }

} // MonteCarlo : Power 


namespace MonteCarlo {

    template < typename R >
    R InvErfc ( const R p ) {

        R x, err, t, pp;

        if ( p >= 2.0 ) return -100;
        if ( p <= 0.0 ) return  100;

        pp = ( p < 1.0 )? p : 2.0 - p;

        t = std::sqrt ( -2.0 * std::log ( pp / 2.0 ) );

        x = -0.70711*((2.30753+t*0.27061)/(1.0+t*(0.99229+t*0.04481))-t );

        for ( auto i = 0; i < 2; i++ ) {

            err = std::erfc(x) - pp;
            x  += err/(1.12837916709551257*std::exp(-x*x) - x*err);

        }

        return ( p < 1.0 )? x : -x;

    }

} // MonteCarlo : InvErfc 


#endif // STATISTIC_DIST_IMPLEMENTATIONS 
