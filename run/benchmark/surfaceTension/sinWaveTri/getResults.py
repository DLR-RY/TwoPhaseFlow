# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
from casefoam import postFunctions
import seaborn as sns

cases = [['isoSurface', 'plicRDF'],
         ['gradAlpha', 'RDF', 'fitParaboloid'],
         ['coarse', 'mid', 'fine']]

baseCase = 'Cases'
solutionDir = 'reconSurfaces'
file = 'alpha.water_freeSurf.raw'

analytical = pd.read_csv("analyticalWave.dat",
                         delim_whitespace=True, header=None)
analytical.columns = ['time', "analytical"]
analytical["analytical"] *= 100

postFunction = postFunctions.getFreeSurfaceWallAndCentre

sol = casefoam.posField_to_timeSeries(
    solutionDir, file, postFunction, cases, baseCase, axis=1)
sol = sol.reset_index()
sol.columns = ['time', 'min', 'mean', 'max',
               'interfaceType', 'Method', 'nCells']
sol = sol.replace('coarse', 32)
sol = sol.replace('mid', 64)
sol = sol.replace('fine', 128)
sol['max'] /= 3e-5
sol['time'] /= 1.475e-5
ax = analytical.plot(x='time', y='analytical', c='black', marker='o')
sns.set_style("ticks")
plicRDF = sol[sol['interfaceType'] == 'plicRDF']
plicRDF.to_csv("sinwaveTri.csv", index=False)


sns.lineplot(x='time', y='max', hue='Method', style="nCells", data=plicRDF,
             hue_order=['gradAlpha', 'fitParaboloid', 'RDF'], ax=ax).set(ylim=(0, 3.5))
plt.ylabel('Relative amplitude')
plt.xlabel('Non-dimensional time')
plt.savefig("sinWaveTri_plicRDF.pdf")

plt.show()
