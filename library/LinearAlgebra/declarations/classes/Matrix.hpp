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


namespace linalg {

    template < typename T >
    using Vector = std::vector<T>;

}


namespace linalg {

    template < typename T >
    /**
      * A class of matrix that stores all of its elements, including zeros 
      * Not necessarily square nor symmetric 
      * 
      * @tparam T real or complex floating point e.g. std::complex<double> 
      */
    class DenseMatrix {

        size_t     nRow_, nCol_;
        Vector<T>  data_;

        public:

        /**
          * Construct matrix with given sizes and initialize elements to zero
          */
        DenseMatrix ( size_t nRow, size_t nCol );

        /**
          * Construct matrix with given sizes and use given vector as data 
          */
        DenseMatrix ( size_t nRow, size_t nCol, Vector<T> data );

        ~DenseMatrix () {}

        size_t nRow() const { return nRow_; }
        size_t nCol() const { return nCol_; }

        /**
          * Function to access matrix element. Read and write access 
          * @return reference to matrix element 
          */
        T& at ( size_t i, size_t j );

        /**
          * Function to access matrix element. Read only 
          * @return value of matrix element 
          */
        T operator() ( size_t i, size_t j ) const;

        /**
          * Multiplication with vector 
          * @return a new vector 
          */
        Vector<T> operator* ( const Vector<T>& v ) const;

        /**
          * Multiplication of self transpose with vector 
          * @return a new vector 
          */
        Vector<T> TransProd ( const Vector<T>& v ) const;

        /**
          * Multiplication of self transpose with other matrix 
          * @return a new matrix 
          */
        DenseMatrix<T> TransProd ( const DenseMatrix<T>& other ) const;

        #ifdef LA_COMPLEX 

        /**
          * Multiplication of self conjugate transpose with vector 
          * @return a new vector 
          */
        Vector<T> ConjTransProd ( const Vector<T>& v ) const;

        /**
          * Multiplication of self conjugate transpose with other matrix 
          * @return a new matrix 
          */
        DenseMatrix<T> ConjTransProd ( const DenseMatrix<T>& other ) const;

        #endif // LA_COMPLEX 

    };


} // linalg : DenseMatrix 


namespace linalg {

    template < typename T >
    /**
      * A subclass of dense matrix that is symmetric 
      * Only stores lower triangular and diagonal element 
      * 
      * @tparam T real or complex floating point e.g. std::complex<double> 
      */
    class DenseSymMatrix {

        size_t    dim_;
        Vector<T> data_;

        public: 

        /**
          * Construct matrix with given dim and initialize elements to zero
          */
        DenseSymMatrix ( const size_t dimension );

        size_t nRow() const { return dim_; }
        size_t nCol() const { return dim_; }

        /**
          * Function to access matrix element. Read and write access 
          * @return reference to matrix element 
          */
        T  operator() ( const size_t i, const size_t j ) const;

        /**
          * Function to access matrix element. Read only 
          * @return value of matrix element 
          */
        T& operator() ( const size_t i, const size_t j );

        size_t dimension () const { return dim_; };

    }; 


} // linalg : DenseSymMatrix


namespace linalg {

    template < typename T >
    /**
      * A subclass of dense matrix that is lower triangular and square 
      * Only stores lower triangular and diagonal element 
      * 
      * @tparam T real or complex floating point e.g. std::complex<double> 
      */
    class DenseLTriangularMatrix {

        size_t    dim_;
        Vector<T> data_;

        public: 

        /**
          * Construct matrix with given dim and initialize elements to zero
          */
        DenseLTriangularMatrix ( const size_t dimension );

        size_t nRow() const { return dim_; }
        size_t nCol() const { return dim_; }

        /**
          * Function to access matrix element. Read and write access 
          * @return reference to matrix element 
          */
        T  operator() ( const size_t i, const size_t j ) const;

        /**
          * Function to access matrix element. Read only 
          * @return value of matrix element 
          */
        T& operator() ( const size_t i, const size_t j );

        size_t dimension () const { return dim_; };

        const Vector<T>& data () { return data_; }

    }; 

} // linalg : DenseLTriangularMatrix 


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

