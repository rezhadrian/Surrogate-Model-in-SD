/**
  * @file Cholesky.hpp 
  *
  * @brief declarations of templated Cholesky decomposition 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef CHOLESKY_DECLARATIONS 
#define CHOLESKY_DECLARATIONS 

#include "LibrariesLoader_LA.hpp" 


namespace linalg {

    template < class LTMatrix, class SymMatrix > 
    /**
      * Cholesky decomposition of SPD matrix e.g. correlation matrix 
      *
      * @tparam LTMatrix lower triangular matrix class 
      * @tparam SymMatrix symmetric matrix class 
      * 
      * @return lower triangular part of the decomposition 
      */
    LTMatrix Cholesky ( const SymMatrix& A );

} // linagl : Cholesky 


#ifndef CHOLESKY_IMPLEMENTATIONS 
    #include "Cholesky_imp.hpp" 
#endif 

#endif // CHOLESKY_DECLARATIONS 

