import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.special import erfc


class sinwave_prosperetti:
    # Input
    def __init__(self, rho1, rho2, wavelength,H0,nu, sigma):
        self.rho1 = rho1  # Density liquid 1
        self.rho2 = rho2  # Density liquid 2
        self.wavelength = wavelength  # Wavelength
        self.H0 = H0  # Initial height
        self.nu = nu  # Kinematic viscosity
        self.sigma = sigma       # Surface tension

        # Constantes
        self.waveNb = (2*math.pi)/self.wavelength
        self.omega0 = abs(
            np.sqrt((self.sigma*self.waveNb**3)/(self.rho1+self.rho2)))
        self.epsilon = (self.nu*self.waveNb**2)/self.omega0
        self.beta = (self.rho1*self.rho2)/(self.rho1+self.rho2)**2

        self.p = [1,
                  -4*self.beta*(self.epsilon*self.omega0)**0.5,
                  2*(1-6*self.beta)*self.epsilon*self.omega0,
                  4*(1-3*self.beta)*(self.epsilon*self.omega0)**(3/2),
                  (1-4*self.beta)*(self.epsilon*self.omega0)**2+self.omega0**2]
        self.r = np.roots(self.p)

        self.Z = [0, 0, 0, 0]
        self.Z[0] = (self.r[1]-self.r[0]) * \
            (self.r[2]-self.r[0])*(self.r[3]-self.r[0])
        self.Z[1] = (self.r[0]-self.r[1]) * \
            (self.r[2]-self.r[1])*(self.r[3]-self.r[1])
        self.Z[2] = (self.r[0]-self.r[2]) * \
            (self.r[1]-self.r[2])*(self.r[3]-self.r[2])
        self.Z[3] = (self.r[0]-self.r[3]) * \
            (self.r[1]-self.r[3])*(self.r[2]-self.r[3])

    def f(self, t, i):
        return np.exp(((self.r[i]**2-self.epsilon*self.omega0)*t)/ \
            self.omega0)*erfc((self.r[i]*t**0.5)/(self.omega0**0.5))


    def a(self, t):
        A1 = (4*(1-4*self.beta)*self.epsilon**2) / \
            (8*(1-4*self.beta)*self.epsilon**2+1)*erfc(np.sqrt(self.epsilon*t))

        sum = 0
        for i in range(0, 4):
            sum += (self.r[i]*self.omega0**2)/(self.Z[i] * \
                (self.r[i]**2-self.epsilon*self.omega0))*self.f(t, i)
        return np.real(A1 + sum)

    def plot(self, duration, N=1001):

        self.result = []
        self.tau = np.linspace(0, duration, N)
        for t in self.tau:
            self.result.append(self.a(t))

        plt.plot(self.tau, self.result)
        plt.show()

    def tabulatedData(self, duration, N=1001):

        self.result = []
        self.tau = np.linspace(0, duration, N)
        for t in self.tau:
            self.result.append(self.a(t))

        data = {'t': self.tau, 'analytical': self.result}
        analytical = pd.DataFrame(data)
        return analytical

# rho1 = 1  # Density liquid 1
# rho2 = 1  # Density liquid 2
# wavelength = 0.003  # Wavelength
# H0 = 3e-5  # Initial height
# nu = 0.001  # Kinematic viscosity
# sigma = 1       # Surface tension
# ana = sinwave_prosperetti(rho1, rho2, wavelength,H0, nu, sigma)
# # ana.plot(25)
# print(ana.a(1))
