/**
  * @file VariableGeneration_imp.hpp
  *
  * @brief implement function to generate correlated random variables 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef VARIABLE_GENERATION_IMPLEMENTATIONS 
#define VARIABLE_GENERATION_IMPLEMENTATIONS 

#ifndef MONTE_CARLO_DECLARATIONS 
    #include "MonteCarlo.hpp" 
#endif 


namespace MonteCarlo {


    template < class LTriangularMatrix, typename Z, typename R >
    Vector<R> CombineRVs ( 

        const LTriangularMatrix& L, 
        Vector<R>& RVs, 
        const Z dim 

    ) {

        Vector<R> result ( RVs.size(), 0.0 );

        Z N = RVs.size() / dim;

        for ( auto k = 0; k < N; k++ ) {
        for ( auto i = 0; i < dim; i++ ) {
        for ( auto j = 0; j < dim; j++ ) {

            result[i+k*dim] += (
                L(i,j) * RVs[k+j*N]
            );


        }
        }
        }

        return result;

    }


} // MonteCarlo : VariableGeneration 


#endif // VARIABLE_GENERATION_IMPLEMENTATIONS  

