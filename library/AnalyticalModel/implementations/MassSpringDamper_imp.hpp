/**
  * @file MassSpringDamper_imp.hpp 
  *
  * @brief 
  * Implementations of Mass-Spring-Damper model 
  * 
  * @anchor _MassSpringDamper_imp_hpp_ 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef MASS_SPRING_DAMPER_IMPLEMENTATIONS 
#define MASS_SPRING_DAMPER_IMPLEMENTATIONS 

#ifndef ANALYTICAL_MODEL_DECLARATIONS 
    #include "AnalyticalModel.hpp" 
#endif 


namespace Analytical {

    template < typename Z, typename R, typename C >
    MassSpringDamper<Z,R,C>::MassSpringDamper (

        const Vector<R> Masses, 
        const Vector<R> Dampers, 
        const Vector<R> Springs

    ) : 

        Masses_  ( Vector<R> ( Masses.size()  ) ), 
        Dampers_ ( Vector<R> ( Dampers.size() ) ), 
        Springs_ ( Vector<R> ( Springs.size() ) ), 

        Dim_ ( Masses.size() )

    {

        std::copy ( 

            Masses.begin(), Masses.end(), 
            Masses_.begin()

        );

        std::copy (

            Dampers.begin(), Dampers.end(), 
            Dampers_.begin()

        );

        std::copy ( 

            Springs.begin(), Springs.end(), 
            Springs_.begin()

        );

    }


    template < typename Z, typename R, typename C > 
    Vector<R> MassSpringDamper<Z,R,C>::ComputeMassMatrix () const {

        Vector<R> result ( Dim_ * Dim_, 0.0 );

        for ( auto i = 0; i < Dim_; i++ ) {

            result[i+i*Dim_] = Masses_[i];

        }

        return result;

    }


    template < typename Z, typename R, typename C >
    Vector<R> MassSpringDamper<Z,R,C>::ComputeDampingMatrix () const {

        if ( Dim_ == 1 ) {

            return Dampers_;

        }

        Vector<R> result ( Dim_ * Dim_, 0.0 );


        result[0] = Dampers_[0] + Dampers_[1];
        result[1] =             - Dampers_[1];


        for ( auto i = 1; i < Dim_ - 1; i++ ) {

            result[i+i*Dim_-1] = -Dampers_[i]; 
            result[i+i*Dim_  ] =  Dampers_[i] + Dampers_[i+1];
            result[i+i*Dim_+1] =              - Dampers_[i+1];

        }

        auto i = Dim_ - 1;

        result[i+i*Dim_-1] = -Dampers_[i]; 
        result[i+i*Dim_  ] =  Dampers_[i];

        return result;

    }


    template < typename Z, typename R, typename C > 
    Vector<R> MassSpringDamper<Z,R,C>::ComputeStiffnessMatrix () const {

        if ( Dim_ == 1 ) {

            return Springs_;

        }

        Vector<R> result ( Dim_ * Dim_, 0.0 );


        result[0] = Springs_[0] + Springs_[1];
        result[1] =             - Springs_[1];


        for ( auto i = 1; i < Dim_ - 1; i++ ) {

            result[i+i*Dim_-1] = -Springs_[i]; 
            result[i+i*Dim_  ] =  Springs_[i] + Springs_[i+1];
            result[i+i*Dim_+1] =              - Springs_[i+1];

        }

        auto i = Dim_ - 1;

        result[i+i*Dim_-1] = -Springs_[i]; 
        result[i+i*Dim_  ] =  Springs_[i];

        return result;

    }


    template < typename Z, typename R, typename C >
    Vector<C> MassSpringDamper<Z,R,C>::ComputeDynamicStiffness ( 

        const R omega 

    ) const {

        Vector<R> result ( Dim_ * Dim_, 0.0 );

        auto MMatrix = ComputeMassMatrix ();
        auto CMatrix = ComputeDampingMatrix ();
        auto KMatrix = ComputeStiffnessMatrix ();

        std::transform (

            MMatrix.begin(), MMatrix.end(), 
            CMatrix.begin(), 
            KMatrix.begin(), 

            result.begin(), 

            [omega](const auto m, const auto c, const auto k) {

                return C ( k - omega * omega * m, omega * c );

            }

        );

        return result;

    }
    

}


#endif // MASS_SPRING_DAMPER_IMPLEMENTATIONS 

