/**
  * @file AnalyticalModel.hpp 
  *
  * @brief 
  * Declarations of collection of analytical models 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef ANALYTICAL_MODEL_DECLARATIONS 
#define ANALYTICAL_MODEL_DECLARATIONS 

#include "LibrariesLoader_AM.hpp" 

/**
  * @namespace Analytical 
  * 
  * @brief 
  * Contains a collection of analytical structural dynamic models. 
  * 
  * @anchor _Analytical_  
  */
namespace Analytical {

    template < typename T >
    using Vector = std::vector<T>;


    template < typename Z, typename R, typename C > 
    /**
      * @class Model 
      *
      * @brief 
      * Abstract class of analytical structural dynamic model. 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating number e.g. double 
      * @tparam C a type of floating complex number e.g. std::complex<float> 
      * 
      * @anchor _Model_ 
      */
    class Model {

        public: 

        /**
          * @brief 
          * Compute structure response to harmonic load. 
          * 
          * @tparam R a type of floating number e.g. double 
          * @tparam C a type of floating complex number e.g. std::complex<float> 
          * 
          * @param Force harmonic load vector 
          * @param omega angular velocity 
          * @param AdditionalSprings offset to default springs 
          * 
          * @return displacement vector 
          */
        virtual Vector<C> ComputeResponse ( 

            const Vector<C>& Force,
            const R omega, 
            const typename Vector<R>::const_iterator FirstSpring, 
            const typename Vector<R>::const_iterator LastSpring

        ) const = 0;

        virtual Z Dim () const = 0;
 
        virtual Vector<R> StiffnessMatrix ( const Vector<R>& Springs ) const = 0;

        virtual Vector<C> DynamicStiffness ( const R omega ) const = 0;

    }; // Model 


    template < typename Z, typename R, typename C > 
    /**
      * @class MassSpringDamper 
      * 
      * @brief 
      * SD model consisting of point masses, springs, and dampers. @n 
      * Implemented in @ref _MassSpringDamper_imp_hpp_ 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating number e.g. double 
      * @tparam C a type of floating complex number e.g. std::complex<float> 
      */
    class MassSpringDamper final : public Model<Z,R,C> {

        Vector<R> MassMatrix_;
        Vector<R> DampingMatrix_;
        Vector<R> StiffnessMatrix_;

        Z Dim_;

        public:

        /**
          * @brief 
          * Generate a mass-spring-damper SD model. 
          * 
          * @tparam R a type of floating number e.g. double 
          * 
          * @param Masses  a vector of point masses 
          * @param Dampers a vector of damping coefficient 
          * @param Springs a vector of spring stiffness
          */
        MassSpringDamper ( 

            const Vector<R> Masses, 
            const Vector<R> Dampers, 
            const Vector<R> Springs 

        );

        Z Dim () const override { return Dim_; }

        /**
          * @brief 
          * Compute structure response to harmonic load. 
          * 
          * @tparam R a type of floating number e.g. double 
          * @tparam C a type of floating complex number e.g. std::complex<float> 
          * 
          * @param Force harmonic load vector 
          * @param omega angular velocity 
          * @param FirstSpring iterator of a vector of springs 
          * @param LastSpring  iterator of a vector of springs 
          * 
          * @return displacement vector 
          */
        Vector<C> ComputeResponse ( 

            const Vector<C>& Force, 
            const R omega,
            const typename Vector<R>::const_iterator FirstSpring, 
            const typename Vector<R>::const_iterator LastSpring 

        ) const override;

        Vector<R> StiffnessMatrix ( const Vector<R>& Springs ) const override;

        Vector<C> DynamicStiffness ( const R omega ) const override;

        private: 

        Vector<R> MassMatrix      ( const Vector<R>& Masses  ) const;
        Vector<R> DampingMatrix   ( const Vector<R>& Dampers ) const;

        Vector<C> DynamicStiffness ( 
            const R omega, const Vector<R>& AdditionalSprings 
        ) const;

    };

} // Analytical 

#ifndef MASS_SPRING_DAMPER_IMPLEMENTATIONS 
    #include "MassSpringDamper_imp.hpp" 
#endif 

#endif // ANALYTICAL_MODEL_DECLARATIONS 

