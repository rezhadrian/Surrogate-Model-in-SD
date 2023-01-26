/**
  * @file LinAlg_imp.hpp 
  *
  * @brief implement linear algebra functionalities used in MonteCarlo 
  * 
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef LINALG_IMPLEMENTATIONS 
#define LINALG_IMPLEMENTATIONS 

#ifndef LINALG_DECLARATIONS 
    #include "LinAlg.hpp" 
#endif 


namespace MonteCarlo {


    template < typename Z, typename R >
    DenseSymMatrix<Z,R>::DenseSymMatrix ( const Z dim ) :
        dim_( dim ),
        dat_( Vector<R> ( (dim+1) * dim / 2, 0 ) )
    {}


    template < typename Z, typename R >
    R DenseSymMatrix<Z,R>::operator() ( const Z i, const Z j ) const {

        if ( i > dim_ - 1 || j > dim_ - 1 ) {

            throw std::runtime_error (
                "DenseSymMatrix: index exceeds dimension"
            );

        }

        if ( j > i ) {
            return operator() ( j, i );
        }

        return dat_[ (i+1) * i / 2 + j ];

    }


    template < typename Z, typename R >
    R& DenseSymMatrix<Z,R>::operator() ( const Z i, const Z j ) {

        if ( i > dim_ - 1 || j > dim_ - 1 ) {

            throw std::runtime_error (
                "DenseSymMatrix: index exceeds dimension"
            );

        }

        if ( j > i ) {
            return operator() ( j, i );
        }

        return dat_[ (i+1) * i / 2 + j ];

    }


} // MonteCarlo : DenseSymMatrix 


namespace MonteCarlo {


    template < typename Z, typename R >
    DenseLTriangularMatrix<Z,R>::DenseLTriangularMatrix ( const Z dim ) :
        dim_( dim ),
        dat_( Vector<R> ( (dim+1) * dim / 2, 0 ) )
    {}


    template < typename Z, typename R >
    R DenseLTriangularMatrix<Z,R>::operator() ( const Z i, const Z j ) const {

        if ( i > dim_ - 1 || j > dim_ - 1 ) {

            throw std::runtime_error (
                "DenseTriangularMatrix: index exceeds dimension"
            );

        }

        if ( j > i ) {
            return 0.0;
        }

        return dat_[ (i+1) * i / 2 + j ];

    }


    template < typename Z, typename R >
    R& DenseLTriangularMatrix<Z,R>::operator() ( const Z i, const Z j ) {

        if ( i > dim_ - 1 || j > dim_ - 1 ) {

            throw std::runtime_error (
                "DenseTriangularMatrix: index exceeds dimension"
            );

        }

        if ( j > i ) {

            throw std::runtime_error (
                "DenseTriangularMatrix: cannot write to upper triangular part"
            );

        }

        return dat_[ (i+1) * i / 2 + j ];

    }


} // MonteCarlo : DenseTriangularMatrix 


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


#endif // LINALG_IMPLEMENTATIONS 

