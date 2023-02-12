/**
  * @file VariableGeneration_imp.hpp
  *
  * @brief 
  * Implementations of functions to generate correlated random variables 
  * 
  * @anchor _VariableGeneration_imp_hpp_ 
  *
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef VARIABLE_GENERATION_IMPLEMENTATIONS 
#define VARIABLE_GENERATION_IMPLEMENTATIONS 

#ifndef MONTE_CARLO_DECLARATIONS 
    #include "MonteCarlo.hpp" 
#endif 



namespace MonteCarlo {

    template < typename Z, typename R >
    void ConvertLHStoStdNorm ( Vector<R>& LHSResult ) {

        boost::math::normal dist ( R(0.0), R(1.0) );

        std::transform (

            LHSResult.begin(), LHSResult.end(),
            LHSResult.begin(),

            [&dist]( const auto m ) {

                return quantile ( dist, m );

            }

        );

    }

} // MonteCarlo : ConvertLHStoStdNorm 


namespace MonteCarlo {

    template < typename R, typename C >
    Vector<C> CombineRVs ( 

        const MatrixXT<R>& Correl, 
        Vector<R>& RVs, 
        const size_t dim 

    ) {

        // obtain lower triangular from cholesky decomp.
        MatrixXT<R> L = Correl.llt().matrixL();

        // wrap RVs vector to perform multiplication
        Eigen::Map<MatrixXT<R>> RVMap ( RVs.data(), dim, RVs.size()/dim );

        // create result vector and wrap to receive multiplication result 
        Vector<C> result ( RVs.size() );
        Eigen::Map<MatrixXT<C>> MResult ( result.data(), dim, RVs.size()/dim);

        MResult = L * RVMap;

        return result;

    }

} // MonteCarlo : CombineRVs 


namespace MonteCarlo {

    template < typename Z, typename R >
    Vector<R> GenerateRVs (

        Vector<R>& StdNormRVs, 
        const MatrixXT<R>& Correl, 
        const Vector< std::function<R(R)> >& ICDFs,
        const Z dim 

    ) {

        // obtain lower triangular from cholesky decomp.
        MatrixXT<R> L = Correl.llt().matrixL();

        // wrap RVs vector to perform multiplication
        Eigen::Map<MatrixXT<R>> RVMap ( 
            StdNormRVs.data(), dim, StdNormRVs.size()/dim 
        );

        auto CorrelatedStdNorm = L * RVMap;


        Vector<R> result ( StdNormRVs.size(), 0.0 );


        boost::math::normal dist ( R(0.0), R(1.0) );

        auto CDF = [dist](const auto m){
            return boost::math::cdf ( dist, m );
        };


        for ( auto i = 0; i < CorrelatedStdNorm.rows(); i++ ) {
        for ( auto j = 0; j < CorrelatedStdNorm.cols(); j++ ) {

            result[i+j*dim] = ICDFs[i] ( CDF ( CorrelatedStdNorm(i,j) ) );

        }
        }

        return result; 

    }

} // MonteCarlo : GenerateRVs 


#endif // VARIABLE_GENERATION_IMPLEMENTATIONS  

