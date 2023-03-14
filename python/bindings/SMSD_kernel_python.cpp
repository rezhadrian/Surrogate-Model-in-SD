#include "pybind11/pybind11.h" 
#include "pybind11/stl.h" 
#include "pybind11/complex.h" 

#include "Surrogate_MassSpringDamper.hpp" 

typedef size_t Z; 
typedef double R;
typedef std::complex<R> C; 

typedef std::vector<R> VectorR;
typedef std::vector<C> VectorC;


PYBIND11_MODULE ( SMSD, m ) {

    m.doc() = "Surrogate Model for Mass Spring Damper System";

    m.def ( "RandomSampling", &MonteCarlo::RandomSampling<Z,R,C> );

    pybind11::class_ < Analytical::MassSpringDamper<Z,R,C> > 
    ( m, "MassSpringDamper" ) 

        .def ( 

            pybind11::init< const VectorR&,
                            const VectorR&, 
                            const VectorR& > () 

        ); 

    pybind11::class_ < MassSpringDamper::Surrogate::DirectMCS >
    ( m, "DirectMCS") 

        .def (

            pybind11::init< const VectorR&, 
                            const VectorR&, 
                            const VectorR&, 
                            const R, 
                            const Z > () 

        ) 

        .def (

            "SetIndices", 
            &MassSpringDamper::Surrogate::DirectMCS::SetIndices, 
            "something"

        )

        .def (

            "ComputeResponse", 
            &MassSpringDamper::Surrogate::DirectMCS::ComputeResponse, 
            "something"

        );

    pybind11::class_< MassSpringDamper::Surrogate::IntrusivePCE > 
    ( m, "IntrusivePCE" )

        .def( 

            pybind11::init< const Analytical::MassSpringDamper<Z,R,C>*, 
                            const R, 
                            const Z > () 

        )

        .def (

            "SetIndices", 
            &MassSpringDamper::Surrogate::IntrusivePCE::SetIndices, 
            "something"

        )

        .def (

            "ComputeResponse", 
            &MassSpringDamper::Surrogate::IntrusivePCE::ComputeResponse, 
            "something"

        ) 

        .def (

            "Train", 
            &MassSpringDamper::Surrogate::IntrusivePCE::Train, 
            "something"

        );

} 

