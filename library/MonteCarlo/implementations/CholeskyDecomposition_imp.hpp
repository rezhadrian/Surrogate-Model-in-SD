/**
  * @file CholeskyDecomposition_imp.hpp
  *
  * @brief implement function to perform Cholesky decomposition 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef CHOLESKY_DECOMPOSITION_IMPLEMENTATIONS 
#define CHOLESKY_DECOMPOSITION_IMPLEMENTATIONS 

#ifndef MONTE_CARLO_DECLARATIONS 
    #include "MonteCarlo.hpp" 
#endif 


namespace MonteCarlo {


    template < class LTriangularMatrix, class SymMatrix >
    /**
      * Using Cholesky-Crout Algorithm 
      *
      * Reference:
      * Rushcel, J. (2016). 
      * Parallel Implementation of The Cholesky Decomposition on CPUs and GPUs. 
      * Porto Alegre. Universidade Federal do Rio Grande do Sul Instituto de 
      * Informatica Curso de Ciencia da Computacao.
      */
    LTriangularMatrix CholeskyDecompose ( const SymMatrix& A ) {

        LTriangularMatrix L ( A.dimension() );

        for ( auto j = 0; j < A.dimension(); j++ ) {

            auto sum = 0.0;

            for ( auto k = 0; k < j; k++ ) {
                sum += L(j,k) * L(j,k);
            }

            L(j,j) = std::sqrt ( A(j,j) - sum );

            if ( !( L(j,j) > 0) ) {

                throw std::runtime_error (
                    "CholeskyDecompose: input matrix is not pos. definite"
                );

            }

            for ( auto i = j+1; i < A.dimension(); i++ ) {
                
                sum = 0.0;

                for ( auto k = 0; k < j; k++ ) {
                    sum += L(i,k) * L(j,k);
                }

                L(i,j) = 1.0 / L(j,j) * ( A(i,j) - sum );

            }

        }

        return L;

    }


} // MonteCarlo : CholeskyDecompose 


#endif // CHOLESKY_DECOMPOSITION_IMPLEMENTATIONS 

