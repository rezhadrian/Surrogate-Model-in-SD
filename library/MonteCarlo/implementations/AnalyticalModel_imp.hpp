/**
  * @file AnalyticalModel.hpp
  *
  * @brief implement analytical Frequency Response Function (FRF) models 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef ANALYTICAL_MODEL_IMPLEMENTATIONS 
#define ANALYTICAL_MODEL_IMPLEMENTATIONS 

#ifndef MONTE_CARLO_DECLARATIONS 
    #include "MonteCarlo.hpp" 
#endif 


namespace MonteCarlo {


    template < typename R, typename C >
    #ifdef MC_COMPLEX 
        C EvaluateFRF ( const R omega, const Vector<C>& params ) {
    #else
        C EvaluateFRF ( const R omega, const Vector<R>& params ) {
    #endif

        R I_c  = 3.976078202199582e-04;
        R A_g  = 0.15;
        R H    = 4;
        R L    = 10;
        R Zeta = 0.02;

        #ifdef MC_COMPLEX 
            auto E  = params[0].real();
            auto Mu = params[1].real();
        #else 
            auto E  = params[0];
            auto Mu = params[1];
        #endif

        if ( E <= 0.0 ) {
            throw std::runtime_error (
                "EvaluateFRF: Elasticity modulus must be positive"
            );
        }

        if ( Mu <= 0.0 ) {
            throw std::runtime_error (
                "EvaluateFRF: Mu must be positive"
            );
        }

        auto k = 24 * E * I_c / ( H * H * H );
        auto m = Mu * A_g * L;

        R omega_n = std::sqrt ( k / m );
        
        C FRF_numerator;
        C FRF_denominator;

        FRF_numerator.real ( omega * omega );

        FRF_denominator.real ( omega_n * omega_n - omega * omega );
        FRF_denominator.imag ( 2 * Zeta * omega_n * omega );

        return FRF_numerator / FRF_denominator;

    }


} // MonteCarlo : EvaluateFRF 


namespace MonteCarlo {


    template < typename Z, typename R, typename C, class Function >
    Vector<C> EvaluateModel ( 

        const R omega, 

        #ifdef MC_COMPLEX
            const Vector<C>& RVs, 
        #else 
            const Vector<R>& RVs, 
        #endif

        const Function& Model, 
        const Z dim 

    ) {

        if ( dim < 1 ) {
            throw std::runtime_error (
                "EvaluateModel: dimension must be positive"
            );
        }

        Z nSample = RVs.size() / dim;

        Vector<C> result ( nSample );

        for ( auto i = 0; i < nSample; i++ ) {

            #ifdef MC_COMPLEX
                Vector<C> params ( 
            #else 
                Vector<R> params ( 
            #endif

                RVs.begin() + i * dim, 
                RVs.begin() + i * dim + dim 

            );

            result[i] = Model ( omega, params );

        }

        return result;

    }


} // MonteCarlo : EvaluateModel 


#endif // ANALYTICAL_MODEL_IMPLEMENTATIONS 

