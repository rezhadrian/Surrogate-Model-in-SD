/**
  * @file Cholesky_imp.hpp 
  *
  * @brief implementations of templated Cholesky decomposition 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezhadr@outlook.com 
  */

#ifndef CHOLESKY_IMPLEMENTATIONS 
#define CHOLESKY_IMPLEMENTATIONS 

#ifndef CHOLESKY_DECLARATIONS 
    #include "Cholesky.hpp" 
#endif 


namespace linalg {

    template < class LTMatrix, class SymMatrix > 
    /**
      * Using Cholesky-Crout Algorithm 
      *
      * Reference:
      * Rushcel, J. (2016). 
      * Parallel Implementation of The Cholesky Decomposition on CPUs and GPUs. 
      * Porto Alegre. Universidade Federal do Rio Grande do Sul Instituto de 
      * Informatica Curso de Ciencia da Computacao.
      */
    LTMatrix Cholesky ( const SymMatrix& A ) {

        if ( A.nRow() != A.nCol() ) {
            throw std::runtime_error (
                "Cholesky: Matrix is not symmetic"
            );
        }

        LTMatrix L ( A.nRow() );

        for ( auto j = 0; j < A.nRow(); j++ ) {

            auto sum = 0.0;

            for ( auto k = 0; k < j; k++ ) {
                sum += L(j,k) * L(j,k);
            }

            L(j,j) = std::sqrt ( A(j,j) - sum );

            if ( !( L(j,j) > 0) ) {

                throw std::runtime_error (
                    "Cholesky: input matrix is not pos. definite"
                );

            }

            for ( auto i = j+1; i < A.nRow(); i++ ) {
                
                sum = 0.0;

                for ( auto k = 0; k < j; k++ ) {
                    sum += L(i,k) * L(j,k);
                }

                L(i,j) = 1.0 / L(j,j) * ( A(i,j) - sum );

            }

        }

        return L;

    }

} // linalg : Cholesky 


#endif // CHOLESKY_IMPLEMENTATIONS 

