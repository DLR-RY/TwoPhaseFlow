import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import optimize
import math

t0 = 1.36  # [s]
tWall = 378.15  # [K]
tSat = 373.15  # [K]
L = 2260.0e3  # [J/kg]
cpV = 2030.  # [J/(kg K)]
rhoV = 0.581  # [kg/m**3]
k = 0.025  # [W/(m K)]
aV = k/(rhoV*cpV)  # [m**2/s]
A = 1e-7  # [m**2]


def beta():
    """Solve the transcendental equation."""
    def f(_beta):
        return _beta*math.exp(_beta**2)*math.erf(_beta) \
               - cpV*(tWall - tSat)/(np.sqrt(np.pi)*L)

    solution = optimize.fsolve(f, 1.)
    return solution[0]


def x(t):
    """Interface position"""
    return 2*beta()*np.sqrt(aV*t)


def dx(t):
    """Interface velocity"""
    return beta()*np.sqrt(aV*t)/t

# plot analytic solution
t = np.linspace(1.36, 100.0, 1000)

plt.plot(t, x(t)*1e3, label='analytic', color='k')
plt.show()
