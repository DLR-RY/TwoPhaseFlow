# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import optimize, integrate
import math
import numpy as np

t0 = 1.36  # [s]
tWall = 378.15  # [K]
tSat = 373.15  # [K]
L = 2260.0e3  # [J/kg]
cpV = 2030.  # [J/(kg K)]
rhoV = 0.581  # [kg/m**3]
k = 0.025  # [W/(m K)]
aV = k/(rhoV*cpV)  # [m**2/s]
A = 1e-7  # [m**2]

class stefanProblem():

    def beta(self):
        """Solve the transcendental equation."""
        def f(_beta):
            return _beta*math.exp(_beta**2)*math.erf(_beta) \
                - cpV*(tWall - tSat)/(np.sqrt(np.pi)*L)

        solution = optimize.fsolve(f, 1.)
        return solution[0]


    def x(self,t):
        """Interface position"""
        return 2*self.beta()*np.sqrt(aV*t)*1000


    def dx(self,t):
        """Interface velocity"""
        return self.beta()*np.sqrt(aV*t)/t

surfPos = pd.read_csv("results/stefanProblem.csv")

analytical_data = pd.read_csv('init/data.csv')
analytical_data.columns = ['t','x(t)']
analytical_data['x(t)'] *= 1000

ana = stefanProblem()
surfPos_plicRDF = surfPos[surfPos.interFaceType == 'plicRDF']

ax = sns.lineplot(x='time', y='x(t)', hue='Method',
                  style='Resolution', data=surfPos_plicRDF,
                  style_order=['grid1','grid2','grid3'],
                  hue_order=['implicitGrad', 'explicitGrad', 'Schrage'])
                  
ax = analytical_data.iloc[::50, :].plot(x='t',y='x(t)',label='analytical',color='black',ax=ax,marker='o', linestyle='None')
ax.set(ylim=(0, 5))
ax.set_ylabel('x(t) [mm]')
ax.set_xlabel('t [s]')
plt.savefig('results/stefan_plicRDF.pdf')
plt.savefig('results/stefan_plicRDF.png')

surfPos_isoSurface = surfPos[surfPos.interFaceType == 'isoSurface']

plt.figure()
ax = sns.lineplot(x='time', y='x(t)', hue='Method',
                  style='Resolution', data=surfPos_isoSurface,
                  style_order=['grid1','grid2','grid3'],
                  hue_order=['implicitGrad', 'explicitGrad', 'Schrage'])
                  
ax = analytical_data.iloc[::50, :].plot(x='t',y='x(t)',label='analytical',color='black',ax=ax,marker='o', linestyle='None')
ax.set(ylim=(0, 5))
ax.set_ylabel('x(t) [mm]')
ax.set_xlabel('t [s]')
plt.savefig('results/stefan_isoSurface.pdf')
plt.savefig('results/stefan_isoSurface.png')


plt.show()
