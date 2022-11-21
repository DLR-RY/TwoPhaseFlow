# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import interpolate

sns.set_style("ticks")


class SuckingInterface():

    def __init__(self):
        dir_name = os.path.dirname(os.path.abspath(__file__))
        data = pd.read_csv(dir_name + '/init/data.dat')
        data.columns = ['t', 'v(t)', 'x(t)']
        data['x(t)'] *= 1000
        self.posInterp = interpolate.interp1d(data['t'], data['x(t)'])

    def x(self, t):
        """Interface position"""
        return self.posInterp(t)


surfPos = pd.read_csv("results/suckingInterface.csv")

analytical_data = pd.read_csv('init/data.dat')
analytical_data.columns = ['t', 'v(t)', 'x(t)']
analytical_data['x(t)'] *= 1000

ana = SuckingInterface()

surfPos_plicRDF = surfPos[surfPos.interFaceType == 'plicRDF']

mask = (surfPos_plicRDF["Resolution"] == 'grid1') & (
    surfPos_plicRDF["Method"] == 'Schrage') & (surfPos_plicRDF['time'] > 0.2)
surfPos_plicRDF.loc[mask, :] = 200

# cut off Schrage models for better visulatiuzation
surfPos_plicRDF = surfPos_plicRDF[surfPos_plicRDF['x(t)'] < 4]
ax = sns.lineplot(x='time', y='x(t)', hue='Method',
                  style='Resolution', data=surfPos_plicRDF,
                  style_order=['grid1', 'grid2', 'grid3'],
                  hue_order=['implicitGrad', 'explicitGrad', 'Schrage'])
ax = analytical_data.iloc[::50, :].plot(
    x='t', y='x(t)', label='analytical', color='black', ax=ax, marker='o')
ax.set(ylim=(0, 4))
ax.set_ylabel('x(t) [mm]')
ax.set_xlabel('t [s]')
plt.savefig('results/suckingInterface_plicRDF.pdf')
plt.savefig('results/suckingInterface_plicRDF.png')

surfPos_isoSurface = surfPos[surfPos.interFaceType == 'isoSurface']

# cut off Schrage models for better visulatiuzation
mask = (surfPos_isoSurface["Resolution"] == 'grid1') & (
    surfPos_isoSurface["Method"] == 'Schrage') & (surfPos_isoSurface['time'] > 0.2)
surfPos_isoSurface.loc[mask, :] = 200

surfPos_isoSurface = surfPos_isoSurface[surfPos_isoSurface['x(t)'] < 4]
plt.figure()
ax = sns.lineplot(x='time', y='x(t)', hue='Method',
                  style='Resolution', data=surfPos_isoSurface,
                  style_order=['grid1', 'grid2', 'grid3'],
                  hue_order=['implicitGrad', 'explicitGrad', 'Schrage'])
ax = analytical_data.iloc[::50, :].plot(
    x='t', y='x(t)', label='analytical', color='black', ax=ax, marker='o')
ax.set(ylim=(0, 4))
ax.set_ylabel('x(t) [mm]')
ax.set_xlabel('t [s]')
plt.savefig('results/suckingInterface_isoSurface.pdf')
plt.savefig('results/suckingInterface_isoSurface.png')


plt.show()
