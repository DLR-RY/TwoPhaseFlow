# -*- coding: utf-8 -*-
import pandas as pd
import math
import matplotlib.pyplot as plt
import casefoam
from casefoam import postFunctions
import seaborn as sns
from sinwave_prosperetti import sinwave_prosperetti
from scipy import integrate

solutionDir = 'surfaces'
file = 'alpha.water_constantIso.raw'

rho1 = 1  # Density liquid 1
rho2 = 1  # Density liquid 2
wavelength = 0.003  # Wavelength
H0 = 3e-5  # Initial height
nu = 0.001  # Kinematic viscosity
sigma = 1       # Surface tension
omega0 = math.sqrt(sigma*(2*math.pi/wavelength)**3/(2*rho1))
ana = sinwave_prosperetti(rho1, rho2, wavelength, H0, nu, sigma)


postFunction = postFunctions.getFreeSurfaceWallAndCentre

sol = casefoam.posField_to_timeSeries(
    solutionDir, file, postFunction, axis=1)
sol = sol.reset_index()
sol.columns = ['time', 'min', 'mean', 'max',
               'interfaceType']
sol['max'] /= 3e-5
sol['time'] /= 1/omega0
sol = sol.sort_values('time')
sol['analytical'] = sol.apply(lambda x: abs(ana.a(x['time'])), axis=1)
sol['error'] = H0*(sol['max'] - sol['analytical'])**2

print(sol)
int_err = integrate.trapz(sol['error'], x=sol['time'])
err = 1/wavelength*math.sqrt(1/25*int_err)
print("err",err)
print("omega0",omega0,1/omega0)

ax = sol.plot(x='time',y='max')
sol.plot(x='time',y='analytical',ax=ax)

plt.show()
