/**
  * @file LatinHypercubeSampling_imp.hpp
  *
  * @brief implement functions to produce latin hypercube sample
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef LATIN_HYPERCUBE_SAMPLING_IMPLEMENTATIONS 
#define LATIN_HYPERCUBE_SAMPLING_IMPLEMENTATIONS 

#ifndef MONTE_CARLO_DECLARATIONS 
    #include "MonteCarlo.hpp" 
#endif 


namespace MonteCarlo {


    template < typename Z, typename R >
    /**
      * Generate random samples using latin hypercube sampling 
      * Variables from the same dimension are contiguous in memory 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      */
    Vector<R> UnsortedLHS ( const Z nPoints, const Z dim );


    template < typename Z >
    /**
      * Generate indices to sort UnsortedLHS 
      * Make variables from each sample point contigous in memory 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      */
    Vector<Z> TransposerIndices ( const Z nPoints, const Z dim );


    template < typename Z, typename R > 
    /**
      * Sort LHS result based on given indices 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      */
    Vector<R> SortLHS ( 

        const Vector<R>& UnsortedLHS, 
        const Vector<Z>& SorterIndices 

    );


} // MonteCarlo : LHS SubFunctions Declarations 


namespace MonteCarlo {

    template < typename Z, typename R >
    Vector<R> LHS ( const Z nPoints, const Z dim ) {

        auto Samples = UnsortedLHS<Z,R> ( nPoints, dim );
        auto SorterIndices = TransposerIndices<Z> ( nPoints, dim );

        return SortLHS<Z,R> ( Samples, SorterIndices );

    }

} // MonteCarlo : LHS 


namespace MonteCarlo {

    template < typename Z, typename R >
    Vector<R> UnsortedLHS ( const Z nPoints, const Z dim ) {

        if ( nPoints < 1 ) {

            throw std::runtime_error (
                "LHS: Number of points must be positive integer"
            );

        }

        if ( dim < 1 ) {

            throw std::runtime_error (
                "LHS: dimension must be positive integer"
            );

        }

        R range = 1.0 / nPoints;

        std::random_device device;

        std::default_random_engine generator ( device () );
        std::mt19937 shuffler ( device () );

        std::uniform_real_distribution <R> RandomVariable ( 0.0, range );


        Vector<Z> IntervalIndices ( nPoints );

        std::iota (

            IntervalIndices.begin(),
            IntervalIndices.end(),
            0

        );

        Vector<R> result ( nPoints * dim );

        for ( auto i = 0; i < dim; i++ ) {

            std::shuffle (

                IntervalIndices.begin(),
                IntervalIndices.end(),
                shuffler

            );

            std::transform (

                IntervalIndices.begin(),
                IntervalIndices.end(),

                result.begin() + i * nPoints,

                [range, &RandomVariable, &generator ]( auto m ) {
                    
                    return m * range + RandomVariable ( generator );
                    
                }

            );

        }

        return result;

    }

} // MonteCarlo : UnsortedLHS 


namespace MonteCarlo {

    template < typename Z >
    Vector<Z> TransposerIndices ( const Z nPoints, const Z dim ) {

        if ( nPoints < 1 ) {
            throw std::runtime_error (
                "TransposerIndices: num of points must be greater than unity"
            );
        }

        if ( dim < 1 ) {
            throw std::runtime_error (
                "TransposerIndices: dimension must be greater than unity"
            );
        }

        Vector<Z> FirstPointIndices ( dim ); 

        std::iota (

            FirstPointIndices.begin(), 
            FirstPointIndices.end(), 
            0

        );

        Vector<Z> SorterIndices ( nPoints * dim );

        for ( auto i = 0; i < nPoints; i++ ) {

            std::transform (

                FirstPointIndices.begin(), 
                FirstPointIndices.end(), 
                SorterIndices.begin() + i * dim, 

                [i,nPoints](const auto m) {
                    
                    return m * nPoints + i;
                    
                }

            );

        }

        return SorterIndices;

    }

} // MonteCarlo : TransposerIndices 


namespace MonteCarlo {

    template < typename Z, typename R > 
    Vector<R> SortLHS ( 

        const Vector<R>& UnsortedLHS, 
        const Vector<Z>& SorterIndices 

    ) {

        Vector<R> SortedLHS ( UnsortedLHS.size() );

        std::transform (

            SorterIndices.begin(), SorterIndices.end(), 

            SortedLHS.begin(), 

            [&UnsortedLHS]( const auto m) {

                return UnsortedLHS[m];

            }

        );

        return SortedLHS;

    }

} // MonteCarlo : SortLHS 


#endif // LATIN_HYPERCUBE_SAMPLING_IMPLEMENTATIONS 

