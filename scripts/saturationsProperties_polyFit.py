import CoolProp.CoolProp as cp
import numpy as np
import matplotlib.pyplot as plt

pressure = np.linspace(0.5e5,10e5,100)
TSat = np.linspace(0.5e5,10e5,100)
L = np.linspace(0.5e5,10e5,100)

for idx,p in np.ndenumerate(pressure):
    H_V = cp.PropsSI('H','P',p,'Q',1,'N2')
    H_L = cp.PropsSI('H','P',p,'Q',0,'N2')
    L[idx[0]] = H_V - H_L

for idx,p in np.ndenumerate(pressure):
    TSat[idx[0]] = cp.PropsSI('T','P',p,'Q',0,'N2')

fit_Tsat = np.polyfit(pressure,TSat, 5)
plt.plot(pressure,TSat,pressure,np.polyval(fit_Tsat,pressure))
plt.draw()
print(fit_Tsat[::-1])
plt.figure()
fit_L = np.polyfit(pressure,L, 5)
print(fit_L[::-1])
plt.plot(pressure,L,pressure,np.polyval(fit_L,pressure))
plt.show()

