from scipy import optimize
import math
import matplotlib.pyplot as plt
import numpy as np

x0 = 0.5e-3  # [m]

tWall = 378.15  # [K]
tSat = 373.15  # [K]
L = 2260.0e3  # [J/kg]
cpV = 2030.  # [J/(kg K)]
rhoV = 0.597  # [kg/m**3]
k =  0.0248  # [W/(m K)]
aV = k/(rhoV*cpV)  # [m**2/s]


def beta():
    """Solve the transcendental equation."""
    def f(_beta):
        return _beta*math.exp(_beta**2)*math.erf(_beta) \
               -cpV*(tWall - tSat)/(np.sqrt(np.pi)*L)

    solution = optimize.fsolve(f, 1.)
    return solution[0]

t0 = x0**2/(4.*beta()**2*aV)
print('startTime = %f s for x0 = %f m' %  (t0, x0))
