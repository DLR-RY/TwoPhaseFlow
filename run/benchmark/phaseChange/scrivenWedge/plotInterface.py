# -*- coding: utf-8 -*-
#%%
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
from casefoam import postFunctions
import seaborn as sns
#output_notebook()
sns.set_style("ticks")

caseStructure = [['isoSurface', 'plicRDF'],
        ['implicitGrad' ,'explicitGrad','interfaceRes'],
        ['grid1', 'grid2', 'grid3', 'grid4']]

baseCase = 'Cases' # 'Cases/isoSurface/explicitGrad/grid1/'
solutionDir = 'surfaces'
file = 'interfaceEnergyFluxLiquid_freeSurf.raw'
postFunction = postFunctions.getRadius

surfPos = casefoam.posField_to_timeSeries(solutionDir, file,postFunction, caseStructure, baseCase)

surfPos.columns = ['min', 'r', 'max', 'interFaceType', 'Method', 'Resolution']
surfPos = surfPos.sort_index()
surfPos = surfPos.reset_index('time')
surfPos[['min', 'r', 'max']] *= 1000
#%%
# plot analytical solution
analytical_data = pd.read_csv('init/data.dat',delim_whitespace=True)
#%%

analytical_data['r'] *= 1000
ax = analytical_data.plot(x='t',y='r',label='analytical',color='black')
surfPos_plicRDF = surfPos[surfPos.interFaceType == 'plicRDF']

ax = sns.lineplot(x='time', y='r', hue='Method',
                  style='Resolution', data=surfPos_plicRDF,
                  ax=ax,style_order=['grid1', 'grid2', 'grid3', 'grid4'],
                  hue_order=['implicitGrad','explicitGrad','interfaceRes'])
ax.set_ylabel('radius [mm]')
ax.set_xlabel('t [s]')
plt.savefig('Scriven_plicRDF.pdf')
surfPos_plicRDF.to_csv("Scriven_plicRDF.csv",index=False)
#%%
ax = analytical_data.plot(x='t',y='r',label='analytical',color='black')
surfPos_isoSurface = surfPos[surfPos.interFaceType == 'isoSurface']

ax = sns.lineplot(x='time', y='r', hue='Method',
                  style='Resolution', data=surfPos_isoSurface,
                  ax=ax,style_order=['grid1', 'grid2', 'grid3', 'grid4'],
                  hue_order=['implicitGrad','explicitGrad','interfaceRes'])
ax.set_ylabel('radius [mm]')
ax.set_xlabel('t [s]')
plt.savefig('Scriven_isoSurface.pdf')
surfPos_isoSurface.to_csv("Scriven_isoSurface.csv",index=False)
plt.show()

# %%


# %%
