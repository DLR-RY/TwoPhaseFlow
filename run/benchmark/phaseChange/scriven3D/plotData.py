# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import interpolate

sns.set_style("ticks")

class scriven3D():

    def __init__(self):
        dir_name = os.path.dirname(os.path.abspath(__file__))
        data = pd.read_csv(dir_name + '/init/data.dat')
        data.columns = ['t','r(t)']
        data['r(t)'] *= 1000
        self.posInterp = interpolate.interp1d(data['t'], data['r(t)'])


    def x(self,t):
        """Interface position"""
        return self.posInterp(t)

surfPos = pd.read_csv("results/scriven3D.csv")

analytical_data = pd.read_csv('init/data.dat')
analytical_data.columns = ['t','r(t)']
analytical_data['r(t)'] *= 1000

ana = scriven3D()
surfPos_plicRDF = surfPos[surfPos.interFaceType == 'plicRDF']

ax = sns.lineplot(x='time', y='r(t)', hue='Method',
                  style='Resolution', data=surfPos_plicRDF,
                  style_order=['grid1','grid2','grid3','grid4'],
                  hue_order=['implicitGrad','interfaceRes','embeddedPhaseChange'])
ax = analytical_data.iloc[::50, :].plot(
    x='t', y='r(t)', label='analytical', color='black', ax=ax, marker='o')

# ax.set(ylim=(0, 4))
ax.set_ylabel('r(t) [mm]')
ax.set_xlabel('t [s]')
plt.savefig('results/scriven3D_plicRDF.pdf')
plt.savefig('results/scriven3D_plicRDF.png')

surfPos_isoSurface = surfPos[surfPos.interFaceType == 'isoSurface']

plt.figure()
ax = sns.lineplot(x='time', y='r(t)', hue='Method',
                  style='Resolution', data=surfPos_isoSurface,
                  style_order=['grid1','grid2','grid3','grid4'],
                  hue_order=['implicitGrad','interfaceRes','embeddedPhaseChange'])
ax = analytical_data.iloc[::50, :].plot(
    x='t', y='r(t)', label='analytical', color='black', ax=ax, marker='o')

# ax.set(ylim=(0, 4))
ax.set_ylabel('r(t) [mm]')
ax.set_xlabel('t [s]')
plt.savefig('results/scriven3D_isoSurface.pdf')
plt.savefig('results/scriven3D_isoSurface.png')


plt.show()
