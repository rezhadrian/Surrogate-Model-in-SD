/**
  * @file GaussSeidel.hpp 
  *
  * @brief declarations of templated Gauss Seidel iterative solver  
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezhadr@outlook.com 
  */

#ifndef GAUSS_SEIDEL_DECLARATIONS 
#define GAUSS_SEIDEL_DECLARATIONS 

#include "LibrariesLoader_LA.hpp" 

namespace linalg {


    template < class Matrix, class Vector, typename T, typename Function >
    /**
      * Iterative Gauss - Seidel solver.
      *
      * @tparam Matrix object with functions nRow(), nCol(), and operator().
      * @tparam Vector object with functions size() and operator[].
      * @tparam T type of real floating point 
      * @tparam Function callable object of type T(Vector).
      */
    struct GaussSeidel {

          size_t   MaxIter_;
               T       Tol_;
               T     omega_;
        Function      Norm_;

        /**
          * Initiate Gauss - Seidel solver.
          *
          * @param MaxIter max number of iterations.
          * @param     Tol allowable norm of error.
          * @param    Norm a measure of vector length.
          * @param   omega relaxation factor
          */
        GaussSeidel ( size_t MaxIter, T Tol, const Function& Norm, T omega = 1 );

        /**
          * Update approximate solution to Ax = b by one iteration.
          * 
          * @param A LHS matrix in equation Ax = b.
          * @param b RHS vector in equation Ax = b.
          * @param x approximate solution which will be updated.
          */
        void UpdateSolution (
            const Matrix& A, const Vector& b, Vector& x
        ) const;

        /**
          * Solve system of linear equations in the form of Ax = b.
          *
          * @param A LHS matrix in equation Ax = b.
          * @param b RHS vector in equation Ax = b.
          * @param x solution to the system of equations.
          */
        void Solve (
            const Matrix& A, const Vector& b, Vector& x
        ) const;

    }; // GaussSeidel 


} // linalg 

#ifndef GAUSS_SEIDEL_IMPLEMENTATIONS 
    #include "GaussSeidel_imp.hpp" 
#endif 

#endif // GAUSS_SEIDEL_DECLARATIONS 

