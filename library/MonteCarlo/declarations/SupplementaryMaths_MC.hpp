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
  * @file SupplementaryMaths_MC.hpp
  *
  * @brief declare additional maths functions for Monte Carlo
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef SUPPLEMENTARY_MATHS_MC_DECLARATIONS 
#define SUPPLEMENTARY_MATHS_MC_DECLARATIONS 

#include "LibrariesLoader_MC.hpp" 

namespace MonteCarlo {

    template < typename Z, typename R >
    R Power ( const Z power, const R number );

} // MonteCarlo 

#ifndef SUPPLEMENTARY_MATHS_MC_IMPLEMENTATIONS 
    #include "SupplementaryMaths_imp_MC.hpp" 
#endif 

#endif // SUPPLEMENTARY_MATHS_MC_DECLARATIONS 

