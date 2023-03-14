import SMSD 
import numpy as np 
import matplotlib.pyplot as plt 

# =============================================================================
# Define analytical model: Mass Spring Damper model 
# =============================================================================

M = [  1.0,  1.0 ] 
C = [  1.0,  0.5 ]
K = [ 20.0, 10.0 ]

F = [  1.0,  1.0 ]

Omega = 10.0 

AnalyticModel = SMSD.MassSpringDamper ( M, C, K ) 

# =============================================================================
# Define random coefficients for model parameters 
# =============================================================================

Dim = 1 

MassBasisCoeffs    = [ 0.000, 0.000, 0.020, 0.000 ]
DampingBasisCoeffs = [ 0.000, 0.000, 0.002, 0.005 ]
SpringBasisCoeffs  = [ 0.000, 0.000, 0.200, 0.300 ]
ForceBasisCoeffs   = [ 0.000, 0.000, 0.000, 0.000 ]

# =============================================================================
# Create and train intrusive PCE surrogate model 
# =============================================================================

dMCS = SMSD.DirectMCS ( AnalyticModel, Omega, Dim )

dMCS.SetIndices ( 1, 1 ) 



# =============================================================================
# Create and train intrusive PCE surrogate model 
# =============================================================================


iPCE = SMSD.IntrusivePCE ( AnalyticModel, Omega, Dim )

iPCE.SetIndices ( 

    4,  # max sum of indices 
    4   # max index

)

iPCE.Train ( 

    F, 
    MassBasisCoeffs, 
    DampingBasisCoeffs, 
    SpringBasisCoeffs, 
    ForceBasisCoeffs 

)

# =============================================================================
# Perform Monte Carlo Simulation 
# =============================================================================

RandomInput = SMSD.RandomSampling ( 50, 1 ) 

dMCSResponse   =   dMCS.ComputeResponse ( 
                        RandomInput, 
                        F, 
                        MassBasisCoeffs, 
                        DampingBasisCoeffs, 
                        SpringBasisCoeffs, 
                        ForceBasisCoeffs 
                    )

iPCEResponse   =   iPCE.ComputeResponse ( RandomInput )

# =============================================================================
# Post - Processing 
# =============================================================================

dMCSu1  = np.array (  dMCSResponse[ ::2])
iPCEu1  = np.array (  iPCEResponse[ ::2] )

plt.plot (

        dMCSu1.real, 
        dMCSu1.imag, 
        '.', color="blue", markersize=3, alpha=0.2, 
        label="dMCS"

)

plt.plot (

        iPCEu1.real, 
        iPCEu1.imag, 
        '.', color="red", markersize=3, alpha=0.2, 
        label="iPCE"

)

plt.legend(loc='best')

plt.show() 

