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


    template < typename C >
    /**
      * Dense matrix store all matrix elements, including zeros.
      * Use row major storage scheme.
      */
    class DenseMatrix {

        size_t nRow_, nCol_;
        Vector<C> data_;

        public:

            /**
              * Initialize DenseMatrix with zero elements.
              */
             DenseMatrix ( size_t nRow, size_t nCol );
            ~DenseMatrix ();

            size_t nRow () const;
            size_t nCol () const;

            /**
              * Function to access matrix element. Usable to read and write.
              * @return reference to matrix element at given indices.
              */
            C& at ( size_t IndexRow, size_t IndexCol );

            /**
              * Function to access matrix element. Read only.
              * @return value of matrix element at given indices.
              */
            C operator() ( size_t IndexRow, size_t IndexCol ) const;

    }; // DenseMatrix 


} // linalg 


#ifndef MATRIX_IMPLEMENTATIONS 
    #include "Matrix_imp.hpp" 
#endif 

#endif // MATRIX_DECLARATIONS 

