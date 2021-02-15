import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt

test = 1

# constants
rhov = 0.581
rhol = 958.4
cpv = 2030.0
cpl = 4216.0
lambdav = 0.025
lambdal = 0.671
L = 2.260e6  # verdampfungsenthalpie

alpha_l = lambdal/(cpl*rhol)

beta = rhov/rhol
C = lambdal/(rhov*L)
B = alpha_l/(C*beta)

# bei 0.1 s


def dphi(y, t, phidot_0):
    """DGL System
    s paper:
    A Volume of Fluid Based Method for Fluid Flows with Phase Change
    von
    Samuel W. J. Welch and John Wilson 2000
    """
    dphi = np.zeros(2)
    dphi[0] = y[1]  # dy/dt = v
    dphi[1] -= (t + phidot_0)*y[1]  # dv/dt = -10
    return dphi


# eta im paper
entDimx = np.arange(0, 10.001, 1e-3)

N = 100

t_start = 0.1 
t_end = t_start + 0.5
T_inf = 378.15  # K überhitzung
dt = 1e-3
t = np.arange(t_start, t_end + dt, dt)
# ist ein bisschen overkill man müsste die DGL nur einmal lösen und
# dann mit eta skalieren
etaCoeff = np.sqrt(0.5/alpha_l)/np.sqrt(t[0])
a = 0.001
b = 0.1
# lösung basiert auf der shooting method kombiniert secant method

for n in range(0, N):
    c = 0.5*(a+b)
    phidot_0 = c

    y0 = [373.15/B, phidot_0]
    Y = odeint(dphi, y0, entDimx, args=(phidot_0,))

    F = Y[-1, 0]*B - T_inf
    if abs(F) < 1e-12:
        print('phidot_0=%f | n=%i | F=%f' % (phidot_0, n, F))
        break  # we are done
    elif F <= 0.0:
        a = c
    else:
        b = c
vs = np.zeros(len(t))

for i in range(0, len(t)):
    etaCoeff = np.sqrt(0.5/alpha_l)/np.sqrt(t[i])
    vs[i] = phidot_0*B*etaCoeff*C

s = vs
s = [0]

for i in range(0, len(t) - 1):
    s.append(s[i] + vs[i]*dt)

# write results in DataFrame
data = pd.DataFrame()
data['t'] = t
data.set_index('t', inplace=True)
data['vs'] = vs
data['s'] = s

plt.figure()
data['s'].plot()

plt.figure()
data['vs'].plot()
# Zeit Geschwindigkeit Position
data.to_csv('data.dat')


# save temperature profile at 0.1
etaCoeff = np.sqrt(0.5/alpha_l)/np.sqrt(t[0])
a = 0.001
b = 0.1

for n in range(0, N):
    c = 0.5*(a+b)
    phidot_0 = c
    # evaluate F(c) = (estimate of u(3) using u(0)=c)
    y0 = [373.15/B, phidot_0]
    Y = odeint(dphi, y0, entDimx, args=(phidot_0,))

    F = Y[-1, 0]*B - T_inf
    if abs(F) < 1e-12:
        break  # we are done
    elif F <= 0.0:
        a = c
    else:
        b = c

plt.figure()
plt.plot(entDimx/etaCoeff*1000, Y[:, 0]*B)


plt.figure()

etaCoeff1_1 = np.sqrt(0.5/alpha_l)/np.sqrt(t[-1])

Tdata = np.zeros([len(entDimx), 2])
Tdata[:, 0] = entDimx/etaCoeff
Tdata[:, 1] = Y[:, 0]*B

# data.iloc[-1,1]
Tdata06 = np.zeros([len(entDimx), 2])
Tdata06[:, 0] = entDimx/etaCoeff1_1 + data.iloc[-1,1] 
Tdata06[:, 1] = Y[:, 0]*B

# temperatur profile
# np.savetxt('T01.dat', Tdata, delimiter=' ')
# np.savetxt('T11.dat', Tdata11, delimiter=' ')

np.savetxt('T01.csv', Tdata, delimiter=',')
np.savetxt('T06.csv', Tdata06, delimiter=',')
#Tdata.to_csv('data.dat')
#Tdata11.to_csv('data.dat')

# plot results
#plt.plot(Tdata[:, 0], Tdata[:, 1])
plt.plot(Tdata06[:, 0], Tdata06[:, 1])
plt.show()

# print('x=%f' % (entDimx[-1]/etaCoeff))
# plt.plot(entDimx, Y[:, 0]*B)
