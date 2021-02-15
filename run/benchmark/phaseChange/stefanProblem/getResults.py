# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
from casefoam import postFunctions
import seaborn as sns

sns.set_style("ticks")

caseStructure = [['isoSurface', 'plicRDF'],
                 ['implicitGrad', 'explicitGrad', 'Schrage'],
                 ['grid1', 'grid2', 'grid3']]

baseCase = 'Cases'
solutionDir = 'surfaces'
file = 'alpha.water_constantIso.raw'
postFunction = postFunctions.getFreeSurfaceWallAndCentre
surfPos = casefoam.posField_to_timeSeries(
    solutionDir, file, postFunction, caseStructure, baseCase)

surfPos.columns = [
    'min', 'x(t)', 'max', 'interFaceType', 'Method', 'Resolution']

surfPos = surfPos.sort_index()
surfPos = surfPos.reset_index('time')
surfPos[['min', 'x(t)', 'max']] *= 1000
surfPos = surfPos.drop(columns=['min',  'max'])

# plot analytical solution
analytical_data = pd.read_csv('init/data.csv')
analytical_data['x(t)'] *= 1000
ax = analytical_data.plot(x='t', y='x(t)', label='analytical', color='black')
surfPos_plicRDF = surfPos[surfPos.interFaceType == 'plicRDF']


ax = sns.lineplot(x='time', y='x(t)', hue='Method',
                  style='Resolution', data=surfPos_plicRDF,
                  ax=ax, style_order=['grid1', 'grid2', 'grid3'],
                  hue_order=['implicitGrad', 'explicitGrad', 'Schrage'])

ax.set_ylabel('x(t) [mm]')
ax.set_xlabel('t [s]')
plt.savefig('Stefan_plicRDF.pdf')
surfPos_plicRDF.to_csv("Stefan_plicRDF.csv", index=False)

ax = analytical_data.plot(x='t', y='x(t)', label='analytical', color='black')
surfPos_isoSurface = surfPos[surfPos.interFaceType == 'isoSurface']

ax = sns.lineplot(x='time', y='x(t)', hue='Method',
                  style='Resolution', data=surfPos_isoSurface,
                  ax=ax, style_order=['grid1', 'grid2', 'grid3'],
                  hue_order=['implicitGrad', 'explicitGrad', 'Schrage'])

ax.set_ylabel('x(t) [mm]')
ax.set_xlabel('t [s]')
plt.savefig('Stefan_isoSurface.pdf')
surfPos_isoSurface.to_csv("Stefan_isoSurface.csv", index=False)

plt.show()
