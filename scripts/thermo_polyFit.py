import CoolProp.CoolProp as cp
import numpy as np
import matplotlib.pyplot as plt

Temp = np.linspace(45,350,100)
capa = np.linspace(0.5e5,10e5,100)
kappa = np.linspace(0.5e5,10e5,100)
mu = np.linspace(0.5e5,10e5,100)
rho = np.linspace(0.5e5,10e5,100)
rho_perf = np.linspace(0.5e5,10e5,100)
p = 11.1e5

def perfectGas(p,T,R):
    return p/(R*T)

for idx,T in np.ndenumerate(Temp):
    capa[idx] = cp.PropsSI('C','P',p,'T',T,'H2')
    kappa[idx] = cp.PropsSI('CONDUCTIVITY','P',p,'T',T,'H2')
    mu[idx] = cp.PropsSI('VISCOSITY','P',p,'T',T,'H2')
    rho[idx] = cp.PropsSI('D','P',p,'T',T,'H2')
    rho_perf[idx] = perfectGas(p,T,4124.2)


fit_cp = np.polyfit(Temp,capa, 7)
plt.plot(Temp,capa,Temp,np.polyval(fit_cp,Temp))
print("fitcp ",fit_cp[::-1])
plt.figure()
fit_kappa = np.polyfit(Temp,kappa, 7)
print("fit kappa ",fit_kappa[::-1])
plt.plot(Temp,kappa,Temp,np.polyval(fit_kappa,Temp))

plt.figure()
fit_mu = np.polyfit(Temp,mu, 7)
print("fit mu ",fit_mu[::-1])
plt.plot(Temp,mu,Temp,np.polyval(fit_mu,Temp))

plt.figure()
fit_rho = np.polyfit(Temp,rho, 7)
print("fit rho ",fit_rho[::-1])
plt.plot(Temp,rho,Temp,rho_perf)

plt.show()
