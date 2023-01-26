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


    template < typename Z, typename R, typename Function >
    /**
      * Apply given Inverse CDF function to a vector of RVs
      * Results overwrites inputs 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * @tparam Function takes one floating point and return another 
      */
    void InverseCDFs ( Vector<R>& RandomVariables, const Function ICDF );


    template < typename Z, typename R, typename Function >
    /**
      * Apply given CDF function to a vector of RVs
      * Results overwrites inputs 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * @tparam Function takes one floating point and return another 
      */
    void ComputeCDFs ( Vector<R>& RandomVariables, const Function CDF );


} // MonteCarlo : InverseCDF 


namespace MonteCarlo {


    template < class LTriangularMatrix, typename Z, typename R >
    /**
      * Combine uncorrelated RVs to generate correlated RVs 
      *
      * @tparam LTriangularMatrix lower triangular matrix class 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * 
      * @param L triangular mat. from cholesky decomp. of correlation mat.  
      */
    Vector<R> CombineRVs ( 

        const LTriangularMatrix& L, 
        Vector<R>& RVs,
        const Z dim

    );


} // MonteCarlo : VariableGeneration 


namespace MonteCarlo {

    template < typename Z, typename R, typename Function >
    void InverseCDFs ( Vector<R>& RandomVariables, const Function ICDF ) {

        std::transform (

            RandomVariables.begin(), RandomVariables.end(),
            RandomVariables.begin(),

            [ICDF]( const auto m ) {
                return ICDF ( m );
            }

        );

    }

} // MonteCarlo : InverseCDFs 


namespace MonteCarlo {


    template < typename Z, typename R >
    void ConvertLHS ( Vector<R>& LHSResult ) {

        auto ICDF = [](const auto m){
            return MonteCarlo::InvStdNormCDF<Z,R> ( m ) ;
        };

        std::transform (

            LHSResult.begin(), LHSResult.end(),
            LHSResult.begin(),

            [ICDF]( const auto m ) {
                return ICDF ( m );
            }

        );

    }


} // MonteCarlo : ConvertLHS 


namespace MonteCarlo {


    template < typename Z, typename R, 

        class SymMatrix, 
        class LTriangularMatrix,
        class Function 

    >
    Vector<R> ConvertStdNorm ( 

        const Vector<R>& StdNormRVs, 
        const SymMatrix& Correlation, 
        const Function & ICDF, 
        const Z dim 

    ) {

        typedef SymMatrix SM;
        typedef LTriangularMatrix LT;

        auto L = CholeskyDecompose<LT,SM> ( Correlation );
        auto result = CombineRVs<LT,Z,R> ( L, StdNormRVs, dim );

        auto CDF = [](const auto m){
            return MonteCarlo::StdNormCDF<Z,R> ( m ) ;
        };

        std::transform (

            result.begin(), result.end(),
            result.begin(),

            [CDF,ICDF]( const auto m ) {
                return ICDF ( CDF ( m ) );
            }

        );

        return result; 

    }


} // MonteCarlo : ConvertStdNorm 


namespace MonteCarlo {

    template < typename Z, typename R, typename Function >
    void ComputeCDFs ( Vector<R>& RandomVariables, const Function CDF ) {

        std::transform (

            RandomVariables.begin(), RandomVariables.end(),
            RandomVariables.begin(),

            [CDF]( const auto m ) {
                return CDF ( m );
            }

        );

    }

} // MonteCarlo : ComputeCDFs 


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

