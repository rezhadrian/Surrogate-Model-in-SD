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

#ifndef MONTE_CARLO_DECLARATIONS 
    #include "MonteCarlo.hpp" 
#endif 


namespace MonteCarlo {

    template < typename Z, typename R, typename Function >
    void InverseCDFs ( Vector<R>& RandomVariables, const Function ICDF ) {

        std::transform (

            RandomVariables.begin(),
            RandomVariables.end(),
            RandomVariables.begin(),

            [ICDF]( const auto m ) {
                return ICDF ( m );
            }

        );

    }

} // MonteCarlo : InverseCDFs 


#endif // INVERSE_CDF_IMPLEMENTATIONS 

