/**
  * This file is part of a Surrogate Model library for structural dynamics. 
  *
  * The library is free: you can distribute it and/or modify it.
  * The library is distributed in the hope that it can be useful and helpful,
  * particularly for students and self - learning individuals.
  * 
  * The library comes WITHOUT ANY WARRANTY: without even the implied warranty 
  * of MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE.
  *
  * The library is licensed under the GNU General Public License v3.0
  *
  *
  * The library is based on works by:
  * - Felix Schneider
  * - Iason Papaioannou
  * 
  * Chair of Structural Mechanics 
  * Engineering Risk Analysis Group 
  *
  * Technische Universitaet Muenchen
  *
  * www.cee.ed.tum.de/bm
  * www.cee.ed.tum.de/era 
  */ 

/**
  * Using Cholesky-Crout Algorithm but inverted to produce upper triangular
  *
  * Reference:
  * Rushcel, J. (2016). 
  * Parallel Implementation of The Cholesky Decomposition on CPUs and GPUs. 
  * Porto Alegre. Universidade Federal do Rio Grande do Sul Instituto de 
  * Informatica Curso de Ciencia da Computacao.
  */

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

#ifndef CHOLESKY_DECOMPOSITION_DECLARATIONS 
    #include "CholeskyDecomposition.hpp" 
#endif 


namespace MonteCarlo {

    template < class UTriangularMatrix, class SymMatrix >
    UTriangularMatrix CholeskyDecompose ( const SymMatrix& A ) {

        UTriangularMatrix L ( A.dimension() );

        for ( auto j = 0; j < A.dimension(); j++ ) {

            auto sum = 0.0;

            for ( auto k = 0; k < j; k++ ) {
                sum += L(k,j) * L(k,j);
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
                    sum += L(k,i) * L(k,j);
                }

                L(j,i) = 1.0 / L(j,j) * ( A(j,i) - sum );

            }

        }

        return L;

    }

} // MonteCarlo : CholeskyDecompose 


#endif // CHOLESKY_DECOMPOSITION_IMPLEMENTATIONS 

