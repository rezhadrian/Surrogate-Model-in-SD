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

    template < typename C >
    using MatrixXC = Eigen::Matrix< C, Eigen::Dynamic, Eigen::Dynamic >;

    template < typename C >
    using VectorXC = Eigen::Vector< C, Eigen::Dynamic>;


    template < typename Z, typename R, typename C >
    MassSpringDamper<Z,R,C>::MassSpringDamper (

        const Vector<R> Masses, 
        const Vector<R> Dampers, 
        const Vector<R> Springs

    ) {

        if ( Masses.size() != Dampers.size() ) {
            throw std::runtime_error (
                "MassSpringDamper: Number of masses and dampers must be equal"
            );
        }

        if ( Masses.size() != Springs.size() ) {
            throw std::runtime_error (
                "MassSpringDamper: Number of masses and springs must be equal"
            );
        }

        Dim_ = Masses.size(); 

        if ( Dim_ == 1 ) {

            MassMatrix_.push_back      (  Masses[0] );
            DampingMatrix_.push_back   ( Dampers[0] );
            StiffnessMatrix_.push_back ( Springs[0] );

            return; 

        }

        // Use private function for conciseness 
        MassMatrix_      = MassMatrix      ( Masses  );
        DampingMatrix_   = DampingMatrix   ( Dampers );
        StiffnessMatrix_ = StiffnessMatrix ( Springs );

    } // Constructor 


    template < typename Z, typename R, typename C > 
    Vector<R> MassSpringDamper<Z,R,C>::MassMatrix (

        const Vector<R>& Masses 

    ) const {

        auto Dim = Masses.size(); 

        Vector<R> MassMatrix ( Dim * Dim, 0.0 );

        for ( auto i = 0; i < Dim; i++ ) {

            MassMatrix[i+i*Dim] = Masses[i];

        }

        return MassMatrix;

    } // MassMatrix 


    template < typename Z, typename R, typename C >
    Vector<R> MassSpringDamper<Z,R,C>::DampingMatrix ( 

        const Vector<R>& Dampers 

    ) const {

        auto Dim = Dampers.size();

        Vector<R> DampingMatrix ( Dim * Dim, 0.0 );

        DampingMatrix[0] = Dampers[0] + Dampers[1];
        DampingMatrix[1] =            - Dampers[1];

        for ( auto i = 1; i < Dim - 1; i++ ) {

            DampingMatrix[i+i*Dim-1] = -Dampers[i]; 
            DampingMatrix[i+i*Dim  ] =  Dampers[i] + Dampers[i+1];
            DampingMatrix[i+i*Dim+1] =             - Dampers[i+1];

        }

        auto i = Dim - 1;

        DampingMatrix[i+i*Dim-1] = -Dampers[i]; 
        DampingMatrix[i+i*Dim  ] =  Dampers[i];

        return DampingMatrix;

    } // DampingMatrix 


    template < typename Z, typename R, typename C > 
    Vector<R> MassSpringDamper<Z,R,C>::StiffnessMatrix ( 

        const Vector<R>& Springs
            
    ) const {

        auto Dim = Springs.size();

        Vector<R> StiffnessMatrix ( Dim * Dim, 0.0 );

        StiffnessMatrix[0] = Springs[0] + Springs[1];
        StiffnessMatrix[1] =            - Springs[1];


        for ( auto i = 1; i < Dim - 1; i++ ) {

            StiffnessMatrix[i+i*Dim-1] = -Springs[i]; 
            StiffnessMatrix[i+i*Dim  ] =  Springs[i] + Springs[i+1];
            StiffnessMatrix[i+i*Dim+1] =             - Springs[i+1];

        }

        auto j = Dim - 1;

        StiffnessMatrix[j+j*Dim-1] = -Springs[j]; 
        StiffnessMatrix[j+j*Dim  ] =  Springs[j];

        return StiffnessMatrix; 

    } // StiffnessMatrix 


    template < typename Z, typename R, typename C >
    Vector<C> MassSpringDamper<Z,R,C>::DynamicStiffness (

        const R omega, const Vector<R>& AdditionalSprings 

    ) const {

        auto AdditionalStiffness = StiffnessMatrix ( AdditionalSprings );

        Vector<C> result ( Dim_ * Dim_, 0.0 );

        for ( auto i = 0; i < result.size(); i++ ) {

            result[i] = C (

                StiffnessMatrix_[i] + AdditionalStiffness[i] -
                omega * omega * MassMatrix_[i]

                , 

                omega * DampingMatrix_[i]

            );

        }

        return result;

    } // DynamicStiffness 

    template < typename Z, typename R, typename C >
    Vector<C> MassSpringDamper<Z,R,C>::DynamicStiffness (

        const R omega

    ) const {

        Vector<C> result ( Dim_ * Dim_, 0.0 );

        for ( auto i = 0; i < result.size(); i++ ) {

            result[i] = C (

                StiffnessMatrix_[i] -
                omega * omega * MassMatrix_[i]

                , 

                omega * DampingMatrix_[i]

            );

        }

        return result;

    } // DynamicStiffness 


    template < typename Z, typename R, typename C >
    Vector<C> MassSpringDamper<Z,R,C>::ComputeResponse ( 

        const Vector<C>& Force, 
        const R omega, 
        const typename Vector<R>::const_iterator FirstSpring, 
        const typename Vector<R>::const_iterator LastSpring 

    ) const {

        Vector<C> result ( Dim_, 0.0 );

        auto KDynamic = DynamicStiffness ( 
            omega, Vector<R> ( FirstSpring, LastSpring ) 
        );

        Eigen::Map<MatrixXC<C>, Eigen::RowMajor> KD ( 

            KDynamic.data(), Dim_, Dim_

        );

        Eigen::Map<const VectorXC<C>> Load ( Force.data(), Dim_ );
        Eigen::Map<VectorXC<C>> Disp ( result.data(), Dim_ );

        Disp = KD.ldlt().solve ( Load );

        return result;

    } // ComputeResponse 

} // Analytical : MassSpringDamper 

#endif // MASS_SPRING_DAMPER_IMPLEMENTATIONS 

