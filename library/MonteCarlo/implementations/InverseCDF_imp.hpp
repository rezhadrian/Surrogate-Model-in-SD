/**
  * This file is part of a Surrogate Model library for structural dynamics. 
  *
  * The library is free: you can distribute it and/or modify it.
  * The library is distributed in the hope that it can be useful and helpful,
  * particularly for students and self - learning individuals.
  * 
  * The library comes WITHOUT ANY WARRANTY: without even the implied warranty 
  * of MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE.
  *
  * The library is licensed under the GNU General Public License v3.0
  *
  *
  * The library is based on works by:
  * - Felix Schneider
  * - Iason Papaioannou
  * 
  * Chair of Structural Mechanics 
  * Engineering Risk Analysis Group 
  *
  * Technische Universitaet Muenchen
  *
  * www.cee.ed.tum.de/bm
  * www.cee.ed.tum.de/era 
  */ 

/**
  * The algorithm used is the Beasley-Springer-Moro algorithm 
  * 
  * Reference:
  *
  * Glasserman, P. (2004). Monte Carlo methods in financial engineering. 
  * New York, NY : Springer.
  */

/**
  * @file InverseCDF_imp.hpp
  *
  * @brief implement functions to compute inverse CDF of standard normal dist
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef INVERSE_CDF_IMPLEMENTATIONS 
#define INVERSE_CDF_IMPLEMENTATIONS 

#ifndef INVERSE_CDF_DECLARATIONS 
    #include "InverseCDF.hpp" 
#endif 


namespace MonteCarlo {

    template < typename Z, typename R >
    R InverseSingleCDF ( const R quantile ) {

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
            return -1.0 * InverseSingleCDF<Z,R> ( 1 - quantile );
        }

        if ( quantile >= 0.5 && quantile <= 0.92 ) {

            R   numerator = 0.0;
            R denominator = 1.0;

            for ( auto i = 0; i < 4; i++ ) {
                  numerator += a[i] * Power<Z,R> ( 2*i + 1, quantile - 0.5 ); 
                denominator += b[i] * Power<Z,R> ( 2*i + 2, quantile - 0.5 ); 
            }

            return numerator / denominator;

        }

        R result = 0.0;

        for ( auto i = 0; i < 9; i++ ) {
            result += c[i] * Power<Z,R> (
                i,
                std::log ( - std::log (1-quantile) ) 
            );
        }

        return result;

    }

} // MonteCarlo : InverseSingleCDF 


#endif // INVERSE_CDF_IMPLEMENTATIONS 

