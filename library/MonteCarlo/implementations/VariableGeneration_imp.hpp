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


} // MonteCarlo : Declarations of SubFunctions 


namespace MonteCarlo {

    template < typename Z, typename R >
    void ConvertLHStoStdNorm ( Vector<R>& LHSResult ) {

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

} // MonteCarlo : ConvertLHStoStdNorm 


namespace MonteCarlo {

    template < typename R, class LTMatrix >
    #ifdef MC_COMPLEX 
        Vector< std::complex<R> > GenerateRVs (
    #else 
        Vector<R> GenerateRVs ( 
    #endif

        Vector<R>& StdNormRVs, 
        const LTMatrix& L, 
        const Vector< std::function<R(R)> >& ICDFs,
        const size_t dim 

    ) {

        size_t nSamples = StdNormRVs.size() / dim;

        typedef LTMatrix LT;

        auto CorrelatedStdNorm = CombineRVs<R,LT> ( L, StdNormRVs, dim );
        
        #ifdef MC_COMPLEX 
            Vector< std::complex<R> > result ( CorrelatedStdNorm.size(), 0.0 );
        #else 
            Vector<R> result ( CorrelatedStdNorm.size(), 0.0 );
        #endif 

        auto CDF = [](const auto m){
            return MonteCarlo::StdNormCDF<R> ( m ) ;
        };

        for ( auto i = 0; i < nSamples; i++ ) {

            std::transform (

                CorrelatedStdNorm.begin() + i * dim,
                CorrelatedStdNorm.begin() + i * dim + dim, 
                ICDFs.begin(),

                result.begin() + i * dim,

                [CDF]( const auto m, const auto& ICDF ) {

                #ifdef MC_COMPLEX 
                    return std::complex <R> ( ICDF ( CDF ( m.real() ) ) );
                #else
                    return ICDF ( CDF ( m ) );
                #endif 

                }

            );

        }

        return result; 

    }

} // MonteCarlo : GenerateRVs 


namespace MonteCarlo {

    template < typename R, class LTriangularMatrix >
    Vector<std::complex<R>> CombineRVs ( 

        const LTriangularMatrix& L, 
        Vector<R>& RVs, 
        const size_t dim 

    ) {

        Vector<std::complex<R>> result ( RVs.size(), 0.0 );

        size_t N = RVs.size() / dim;

        for ( auto k = 0; k < N; k++ ) {
        for ( auto i = 0; i < dim; i++ ) {
        for ( auto j = 0; j < dim; j++ ) {

            result[i+k*dim] +=

                std::complex<R>(L(i,j) * RVs[j + k * dim]);

        }
        }
        }

        return result;

    }

} // MonteCarlo : CombineRVs 


#endif // VARIABLE_GENERATION_IMPLEMENTATIONS  

