/**
  * @file Matrix_imp.hpp 
  *
  * @brief file contains implementations of templated Matrix classes.
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezhadr@outlook.com 
  */

#ifndef MATRIX_IMPLEMENTATIONS 
#define MATRIX_IMPLEMENTATIONS 

#ifndef MATRIX_DECLARATIONS 
    #include "Matrix.hpp" 
#endif 


namespace linalg {


    template < typename C >
    DenseMatrix<C>::DenseMatrix ( size_t nRow, size_t nCol ) :
        nRow_( nRow ),
        nCol_( nCol ),
        data_( Vector<C> (nRow*nCol, 0.0) )
    {}


    template < typename C >
    DenseMatrix<C>::~DenseMatrix () {}


    template < typename C >
    size_t DenseMatrix<C>::nRow () const { return nRow_; }


    template < typename C >
    size_t DenseMatrix<C>::nCol () const { return nCol_; }


    template < typename C >
    C& DenseMatrix<C>::at ( size_t IndexRow, size_t IndexCol ) {
        return data_[IndexRow * nCol_ + IndexCol];
    }


    template < typename C >
    C DenseMatrix<C>::operator() ( size_t IndexRow, size_t IndexCol ) const {
        return data_[IndexRow * nCol_ + IndexCol];
    }


} // linalg : DenseMatrix 


#endif // MATRIX_IMPLEMENTATIONS 

