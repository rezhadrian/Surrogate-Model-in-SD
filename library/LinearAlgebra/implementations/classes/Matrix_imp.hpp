/**
  * @file Matrix_imp.hpp 
  *
  * @brief file contains implementations of templated Matrix classes.
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef MATRIX_IMPLEMENTATIONS 
#define MATRIX_IMPLEMENTATIONS 

#include <stdexcept>
#ifndef MATRIX_DECLARATIONS 
    #include "Matrix.hpp" 
#endif 


// Implementations of DenseMatrix class 

namespace linalg {


    template < typename T >
    DenseMatrix<T>::DenseMatrix ( size_t nRow, size_t nCol ) :
        nRow_( nRow ),
        nCol_( nCol ),
        data_( Vector<T> (nRow*nCol, 0.0) )
    {}


    template < typename T >
    DenseMatrix<T>::DenseMatrix ( size_t nRow, size_t nCol, Vector<T> data ) :
        nRow_( nRow ),
        nCol_( nCol ),
        data_( std::move(data) ) 
    {

        if ( nRow_ * nCol_ != data_.size() ) {
            throw std::runtime_error (
                "DenseMatrix: data size does not match row and column sizes"
            );
        }

    }


    template < typename T >
    T& DenseMatrix<T>::at ( size_t i, size_t j ) {

        if ( i > nRow_ || j > nCol_ ) {
            throw std::runtime_error (
                "DenseMatrix: attempt to access non - existent element"
            );
        }

        return data_[ i * nCol() + j ];
    }


    template < typename T >
    T DenseMatrix<T>::operator() ( size_t i, size_t j ) const {

        if ( i > nRow() || j > nCol() ) {
            throw std::runtime_error (
                "DenseMatrix: attempt to access non - existent element"
            );
        }

        return data_[ i * nCol() + j ];
    }


    template < typename T >
    Vector<T> DenseMatrix<T>::operator* ( const Vector<T>& v ) const {

        if ( nCol_ != v.size() ) {
            throw std::runtime_error (
                "Matrix multiplication: sizes don't match"
            );
        }

        Vector<T> result ( nRow_, 0.0 );

        for ( auto i = 0; i < nRow_; i++ ) {
        for ( auto j = 0; j < nCol_; j++ ) {

            result[i] += this->operator()(i,j) * v[j];

        }
        }

        return result;
    }


    template < typename T >
    Vector<T> DenseMatrix<T>::TransProd ( const Vector<T>& v ) const {

        if ( nRow_ != v.size() ) {
            throw std::runtime_error (
                "Matrix TransProd: sizes don't match"
            );
        }

        Vector<T> result ( nCol_, 0.0 );

        for ( auto j = 0; j < nCol_; j++ ) {
        for ( auto i = 0; i < nRow_; i++ ) {

            result[j] += this->operator()(i,j) * v[i];

        }
        }

        return result;
    }


    template < typename T >
    DenseMatrix<T> DenseMatrix<T>::TransProd ( 

        const DenseMatrix<T>& other 

    ) const {

        if ( nRow_ != other.nRow() ) {
            throw std::runtime_error (
                "Matrix TransProd: sizes don't match"
            );
        }

        DenseMatrix<T> result ( nCol_, other.nCol() );

        for ( auto i = 0; i < nCol_       ; i++ ) {
        for ( auto j = 0; j < other.nCol(); j++ ) {
        for ( auto k = 0; k < nRow_       ; k++ ) {

            result.at (i,j) += this->operator()(k,i) * other(k,j);

        }
        }
        }

        return result;
    }


    #ifdef LA_COMPLEX 

    template < typename T >
    Vector<T> DenseMatrix<T>::ConjTransProd ( const Vector<T>& v ) const {

        if ( nRow_ != v.size() ) {
            throw std::runtime_error (
                "Matrix ConjTransProd: sizes don't match"
            );
        }

        Vector<T> result ( nCol_, 0.0 );

        for ( auto j = 0; j < nCol_; j++ ) {
        for ( auto i = 0; i < nRow_; i++ ) {

            result[j] += std::conj(this->operator()(i,j)) * v[i];

        }
        }

        return result;
    }


    template < typename T >
    DenseMatrix<T> DenseMatrix<T>::ConjTransProd ( 

        const DenseMatrix<T>& other 

    ) const {

        if ( nRow_ != other.nRow() ) {
            throw std::runtime_error (
                "Matrix ConjTransProd: sizes don't match"
            );
        }

        DenseMatrix<T> result ( nCol_, other.nCol() );

        for ( auto i = 0; i < nCol_; i++ ) {
        for ( auto j = 0; j < other.nCol(); j++ ) {
        for ( auto k = 0; k < nRow_; k++ ) {

            result.at (i,j) += 
                std::conj( this->operator()(k,i) ) * other(k,j);

        }
        }
        }

        return result;
    }

    #endif // LA_COMPLEX 

} // linalg : DenseMatrix 


// Implementations of DenseSymMatrix class 

namespace linalg {


    template < typename T >
    DenseSymMatrix<T>::DenseSymMatrix ( const size_t dim ) :
        dim_( dim ),
        data_( Vector<T> ( (dim+1) * dim / 2, 0 ) )
    {}


    template < typename T >
    T DenseSymMatrix<T>::operator() ( const size_t i, const size_t j ) const {

        if ( i > dim_ - 1 || j > dim_ - 1 ) {

            throw std::runtime_error (
                "DenseSymMatrix: index exceeds dimension"
            );

        }

        if ( j > i ) {
            return operator() ( j, i );
        }

        return data_[ (i+1) * i / 2 + j ];

    }


    template < typename T >
    T& DenseSymMatrix<T>::operator() ( const size_t i, const size_t j ) {

        if ( i > dim_ - 1 || j > dim_ - 1 ) {

            throw std::runtime_error (
                "DenseSymMatrix: index exceeds dimension"
            );

        }

        if ( j > i ) {
            return operator() ( j, i );
        }

        return data_[ (i+1) * i / 2 + j ];

    }


} // linalg : DenseSymMatrix 


// Implementations of DenseLTriangularMatrix class 

namespace linalg {


    template < typename T >
    DenseLTriangularMatrix<T>::DenseLTriangularMatrix ( const size_t dim ) :
        dim_( dim ),
        data_( Vector<T> ( (dim+1) * dim / 2, 0 ) )
    {}


    template < typename T >
    T DenseLTriangularMatrix<T>::operator() ( 

        const size_t i, const size_t j 

    ) const {

        if ( i > dim_ - 1 || j > dim_ - 1 ) {

            throw std::runtime_error (
                "DenseLTriangularMatrix: index exceeds dimension"
            );

        }

        if ( j > i ) {
            return 0.0;
        }

        return data_[ (i+1) * i / 2 + j ];

    }


    template < typename T >
    T& DenseLTriangularMatrix<T>::operator() ( const size_t i, const size_t j ) {

        if ( i > dim_ - 1 || j > dim_ - 1 ) {

            throw std::runtime_error (
                "DenseLTriangularMatrix: index exceeds dimension"
            );

        }

        if ( j > i ) {

            throw std::runtime_error (
                "DenseLTriangularMatrix: cannot write to upper triangular part"
            );

        }

        return data_[ (i+1) * i / 2 + j ];

    }

} // linalg : DenseLTriangularMatrix 


// Implementations of DiagonalMatrix class 

namespace linalg {


    template < typename T >
    DiagonalMatrix<T>::DiagonalMatrix ( size_t nRow, size_t nCol ) :
        nRow_( nRow ),
        nCol_( nCol ),
        data_( Vector<T> (std::min(nRow,nCol), 0.0) )
    {}


    template < typename T >
    T& DiagonalMatrix<T>::at ( size_t i, size_t j ) {

        if ( i > nRow_ || j > nCol_ || i != j) {

            throw std::runtime_error (
                "DiagonalMatrix: attempt to access non - existent element"
            );

        }

        return data_[i];

    }


    template < typename T >
    T DiagonalMatrix<T>::operator() ( size_t i, size_t j ) const {

        if ( i > nRow_ || j > nCol_ ) {

            throw std::runtime_error (
                "DiagonalMatrix: attempt to access non - existent element"
            );

        }

        if ( i == j ) {
            return data_[i];
        }

        return 0.0;

    }


    template < typename T >
    Vector<T> DiagonalMatrix<T>::operator* ( const Vector<T>& v ) const {

        if ( nCol_ != v.size() ) {

            throw std::runtime_error (
                "Matrix multiplication: sizes don't match"
            );

        }

        Vector<T> result ( nRow_, 0.0 );

        for ( auto i = 0; i < nRow_; i++ ) {
        for ( auto j = 0; j < nCol_; j++ ) {

            result[i] += this->operator()(i,j) * v[j];

        }
        }

        return result;

    }


} // linalg : DiagonalMatrix 


#endif // MATRIX_IMPLEMENTATIONS 

