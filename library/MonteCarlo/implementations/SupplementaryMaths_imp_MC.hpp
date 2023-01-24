/**
  * @file SupplementaryMaths_imp_MC.hpp
  *
  * @brief implement additional maths functions for Monte Carlo
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef SUPPLEMENTARY_MATHS_MC_IMPLEMENTATIONS 
#define SUPPLEMENTARY_MATHS_MC_IMPLEMENTATIONS 

#ifndef MONTE_CARLO_DECLARATIONS 
    #include "MonteCarlo.hpp" 
#endif 


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


#endif // SUPPLEMENTARY_MATHS_MC_IMPLEMENTATIONS 

