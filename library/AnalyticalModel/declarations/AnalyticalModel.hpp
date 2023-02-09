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

        virtual ~Model () = 0;

        /**
          * @brief 
          * Compute mass matrix from structural dynamic model. @n 
          * Pure virtual function. 
          * 
          * @tparam R a type of floating number e.g. double 
          * 
          * @return vector suitable for row-major matrix wrapper 
          */
        virtual Vector<R> ComputeMassMatrix () const = 0; 

        /**
          * @brief 
          * Compute damping matrix from structural dynamic model. @n 
          * Pure virtual function. 
          * 
          * @tparam R a type of floating number e.g. double 
          * 
          * @return vector suitable for row-major matrix wrapper 
          */
        virtual Vector<R> ComputeDampingMatrix () const = 0; 

        /**
          * @brief 
          * Compute stiffness matrix from structural dynamic model. @n 
          * Pure virtual function. 
          * 
          * @tparam R a type of floating number e.g. double 
          * 
          * @return vector suitable for row-major matrix wrapper 
          */
        virtual Vector<R> ComputeStiffnessMatrix () const = 0; 

        /**
          * @brief 
          * Compute dynamic stiffness from structural dynamic model. @n 
          * Pure virtual function. 
          * 
          * @tparam R a type of floating number e.g. double 
          * @tparam C a type of floating complex number e.g. std::complex<float> 
          * 
          * @return vector suitable for row-major matrix wrapper 
          */
        virtual Vector<C> ComputeDynamicStiffness ( const R omega ) const = 0;

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
    class MassSpringDamper : public Model<Z,R,C> {

        Vector<R> Masses_;
        Vector<R> Dampers_;
        Vector<R> Springs_;

        Z Dim_;

        public:

        /**
          * @brief 
          * Generate a mass-spring-damper SD model. 
          * 
          * @tparam R a type of floating number e.g. double 
          * 
          * @param Masses a vector of point masses 
          * @param Dampers a vector of damping coefficient 
          * @param Springs a vector of spring stiffness
          */
        MassSpringDamper ( 

            const Vector<R> Masses, 
            const Vector<R> Dampers, 
            const Vector<R> Springs 

        );

        ~MassSpringDamper () override {};

        /**
          * @brief 
          * Compute mass matrix from Mass-Spring-Damper model. 
          * 
          * @tparam R a type of floating number e.g. double 
          * 
          * @return vector suitable for row-major matrix wrapper 
          */
        Vector<R> ComputeMassMatrix () const override;

        /**
          * @brief 
          * Compute damping matrix from Mass-Spring-Damper model.
          * 
          * @tparam R a type of floating number e.g. double 
          * 
          * @return vector suitable for row-major matrix wrapper 
          */
        Vector<R> ComputeDampingMatrix () const override;

        /**
          * @brief 
          * Compute stiffness matrix from Mass-Spring-Damper model. 
          * 
          * @tparam R a type of floating number e.g. double 
          * 
          * @return vector suitable for row-major matrix wrapper 
          */
        Vector<R> ComputeStiffnessMatrix () const override;

        /**
          * @brief 
          * Compute dynamic stiffness from Mass-Spring-Damper model. 
          * 
          * @tparam R a type of floating number e.g. double 
          * @tparam C a type of floating complex number e.g. std::complex<float> 
          * 
          * @return vector suitable for row-major matrix wrapper 
          */
        Vector<C> ComputeDynamicStiffness ( const R omega ) const override;

    };

} // Analytical 

#ifndef MASS_SPRING_DAMPER_IMPLEMENTATIONS 
    #include "MassSpringDamper_imp.hpp" 
#endif 

#endif // ANALYTICAL_MODEL_DECLARATIONS 

