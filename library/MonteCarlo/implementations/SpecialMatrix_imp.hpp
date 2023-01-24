/**
  * @file SpecialMatrix_imp.hpp
  *
  * @brief implement special matrix class for Monte Carlo 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef SPECIAL_MATRIX_IMPLEMENTATIONS 
#define SPECIAL_MATRIX_IMPLEMENTATIONS 

#ifndef MONTE_CARLO_DECLARATIONS 
    #include "MonteCarlo.hpp" 
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
    DenseUTriangularMatrix<Z,R>::DenseUTriangularMatrix ( const Z dim ) :
        dim_( dim ),
        dat_( Vector<R> ( (dim+1) * dim / 2, 0 ) )
    {}


    template < typename Z, typename R >
    R DenseUTriangularMatrix<Z,R>::operator() ( const Z i, const Z j ) const {

        if ( i > dim_ - 1 || j > dim_ - 1 ) {

            throw std::runtime_error (
                "DenseUTriangularMatrix: index exceeds dimension"
            );

        }

        if ( i > j ) {
            return 0.0;
        }

        return dat_[ i*dim_ - (i-1)*i/2 + j-i ];

    }

    
    template < typename Z, typename R >
    R& DenseUTriangularMatrix<Z,R>::operator() ( const Z i, const Z j ) {

        if ( i > dim_ - 1 || j > dim_ - 1 ) {

            throw std::runtime_error (
                "DenseUTriangularMatrix: index exceeds dimension"
            );

        }

        if ( i > j ) {

            throw std::runtime_error (
                "DenseUTriangularMatrix: cannot write to lower triangular part"
            );

        }

        return dat_[ i*dim_ - (i-1)*i/2 + j-i ];

    } 


} // MonteCarlo : DenseUTriangularMatrix 


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


#endif // SPECIAL_MATRIX_IMPLEMENTATIONS 

