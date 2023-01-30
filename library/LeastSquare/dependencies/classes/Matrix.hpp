/**
  * @file Matrix.hpp 
  *
  * @brief file contains declarations of templated Matrix classes.
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezhadr@outlook.com 
  */

#ifndef MATRIX_DECLARATIONS 
#define MATRIX_DECLARATIONS 

#include "LibrariesLoader_LA.hpp" 

template < typename T >
using Vector = std::vector<T>;


namespace linalg {

    template < typename T >
    using Vector = std::vector<T>;


    template < typename T >
    /**
      * A class of matrix that stores all of its elements, including zeros.
      */
    class DenseMatrix {

        protected: 

            size_t     nRow_, nCol_;
            Vector<T>  data_;

        public:

            /**
              * Construct matrix with given size and initialize elements to zero.
              */
             DenseMatrix ( size_t nRow, size_t nCol );
             DenseMatrix ( size_t nRow, size_t nCol, Vector<T> data );
            ~DenseMatrix () {}

            size_t nRow() const { return nRow_; }
            size_t nCol() const { return nCol_; }

            /**
              * Function to access matrix element. Read and write access.
              * @return reference to matrix element.
              */
            T& at ( size_t i, size_t j );

            /**
              * Function to access matrix element. Read only.
              * @return value of matrix element.
              */
            T operator() ( size_t i, size_t j ) const;

            /**
              * Multiplication with vector.
              * @return a new vector.
              */
            Vector<T> operator* ( const Vector<T>& v ) const;

            Vector<T> TransProd ( const Vector<T>& v ) const;

            #ifdef COMPLEX 
            Vector<T> ConjTransProd ( const Vector<T>& v ) const;
            DenseMatrix<T> ConjTransProd ( const DenseMatrix<T>& other ) const;
            #endif 

    }; // DenseMatrix 

} // linalg


namespace linalg {

    template < typename T >
    using Vector = std::vector<T>;


    template < typename T >
    /**
      * Special class of matrix that stores only main diagonal elements.
      * Can be non - square matrix.
      */
    class DiagonalMatrix {

        protected:

            size_t     nRow_, nCol_;
            Vector<T>  data_;

        public:

            /**
              * Construct matrix with given size and initialize elements to zero.
              */
             DiagonalMatrix ( size_t nRow, size_t nCol );
            ~DiagonalMatrix () {}

            size_t nRow() const { return nRow_; }
            size_t nCol() const { return nCol_; }

            /**
              * Function to access matrix element. Read and write access.
              * @return reference to matrix element.
              */
            T& at ( size_t i, size_t j );

            /**
              * Function to access matrix element. Read only.
              * @return value of matrix element.
              */
            T operator() ( size_t i, size_t j ) const;

            /**
              * Multiplication with vector.
              * @return a new vector.
              */
            Vector<T> operator* ( const Vector<T>& v ) const;

    }; // DiagonalMatrix 

}


#ifndef MATRIX_IMPLEMENTATIONS 
    #include "Matrix_imp.hpp" 
#endif 

#endif // MATRIX_DECLARATIONS 

