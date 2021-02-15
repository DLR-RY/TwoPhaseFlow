# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np
import pandas as pd

rhov = 0.5850
rhol = 958.4
cpv = 2030
cpl = 4216
lambdav = 0.0203
lambdal = 0.67135
L = 2.26e6

Tsat = 373.15
Tinf = Tsat + 1

RadiusInit = 50e-6
RadiusDomain = 250e-06

Ja = cpl*(Tinf-Tsat)/L


def betaf(x, betaG, rhov, rhol):
    return np.exp(- betaG**2*(((1 - x)**-2) - 2*(1 - rhov/rhol)*x - 1))


betaGValue = 0.0

q, err = quad(betaf, 0, 1, args=(betaGValue, rhov, rhol,))
# residuum
res = 2*betaGValue**2*q - (rhol*cpl*(Tinf - Tsat)/(rhov*(L + (cpl - cpv)
                                                         * (Tinf - Tsat))))

i = 0
dx = 1

while abs(res) > 1e-10 and i < 1000:
    i += 1
    q, err = quad(betaf, 0, 1, args=(betaGValue, rhov, rhol,))

    res = 2*betaGValue**2*q - (rhol*cpl*(Tinf - Tsat)/(rhov*(L + (cpl - cpv)
                                                             * (Tinf - Tsat))))
    if res > 0:
        betaGValue = betaGValue - dx
        dx = dx/10
    else:
        betaGValue = betaGValue + dx

print("betaGValue = ", betaGValue)
print("residium = ", res)

# calculate scriven bubble radius
alpha_v = lambdal/(cpl*rhol)
ti = RadiusInit**2/(4*betaGValue**2*alpha_v)
t = np.arange(ti, 0.5e-2 + 1e-6, 1e-6)
RadiusScriven = 2*betaGValue*(alpha_v*t)**0.5

sol = pd.DataFrame()
sol['t'] = t
sol['r'] = RadiusScriven
sol.set_index('t', inplace=True)

# calculate temperature profile which corresponds to initial time step
N = 1e4
dr = RadiusDomain/N
drInit = RadiusDomain/N
r = np.arange(RadiusInit, RadiusDomain + dr, dr)
T = pd.DataFrame(columns=['r', 'T'])

for i, ri in enumerate(r):
    qT, err = quad(betaf, 1 - RadiusInit/ri, 1,
                   args=(betaGValue, rhov, rhol,))
    Ti = Tinf - 2*betaGValue**2*((rhov*(L + (cpl - cpv)*(Tinf - Tsat)))
                                 / (rhol*cpl))*qT
    T.loc[i, 'r'] = ri -RadiusInit
    T.loc[i, 'T'] = Ti
    if Ti < Tsat:
        T.loc[i, 'T'] = Tsat


# write data to file
T.set_index('r', inplace=True)
T.loc[-100] = Tsat
T.loc[100] = Tinf
T.sort_index(inplace=True)

T.to_csv('T01.csv', sep=',', header=False)
sol.to_csv('data.dat', sep=' ')

# radius over time
print("betaG ", betaGValue)
print("alphal ", alpha_v)

plt.title('Radius')
plt.plot(t, RadiusScriven)
plt.xlabel('Time $t$ [s]')
plt.ylabel('Radius $r$ [m]')
plt.tight_layout()
plt.show()

T.plot(legend=False, title='Temperature profile at $t_0$')
plt.xlabel('Radius $r$ [m]')
plt.ylabel('Temperature $T$ [K]')
plt.tight_layout()
plt.show()
