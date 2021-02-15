# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
from casefoam import postFunctions
import seaborn as sns
import numpy as np


def userFunction(caseComb, time, currentDataFrame, axis):
    t = time
    err = (currentDataFrame.iloc[:, axis]-5)/5
    L1 = ((currentDataFrame.iloc[:, axis]-5)/5).mean()
    L2 = ((currentDataFrame.iloc[:, axis]-5)**2/5).mean()
    LMax = err.abs().max()
    df = pd.DataFrame(np.array([time, L1, L2, LMax], ndmin=2),
                      columns=['time', 'L1', 'L2', 'LMax'])
    df = df.set_index('time')
    return df


case = [['plicRDF'],
        ['gradAlpha', 'heightFunction', 'RDF', 'fitParaboloid'],
        ['Grid1', 'Grid2', 'Grid3', 'Grid4', 'Grid5', 'Grid6', 'Grid7', 'Grid8', 'Grid9', 'Grid10']]

baseCase = 'Cases'
solutionDir = 'reconSurfaces'
file = 'K__freeSurf.raw'


sol = casefoam.posField_to_timeSeries(
    solutionDir, file, userFunction, case, baseCase, axis=3)
sol = sol.reset_index()

res = np.linspace(10, 128, 10)
res = res.astype(int)
sol = sol.replace('Grid1', res[0])
sol = sol.replace('Grid2', res[1])
sol = sol.replace('Grid3', res[2])
sol = sol.replace('Grid4', res[3])
sol = sol.replace('Grid5', res[4])
sol = sol.replace('Grid6', res[5])
sol = sol.replace('Grid7', res[6])
sol = sol.replace('Grid8', res[7])
sol = sol.replace('Grid9', res[8])
sol = sol.replace('Grid10', res[9])
print(sol.head())
L1_sol = sol.groupby(['var_1', 'var_2']).mean()
L1_sol = L1_sol.reset_index()
L1_sol = L1_sol[['var_1', 'var_2', 'L1']]
L1_sol.columns = ['Method', 'Res', 'value']
L1_sol['Error'] = 'L1'

L2_sol = sol.groupby(['var_1', 'var_2']).mean()
L2_sol = L2_sol.reset_index()
L2_sol = L2_sol[['var_1', 'var_2', 'L2']]
L2_sol.columns = ['Method', 'Res', 'value']
L2_sol['Error'] = 'L2'

LMax_sol = sol.groupby(['var_1', 'var_2']).max()
LMax_sol = LMax_sol.reset_index()
LMax_sol = LMax_sol[['var_1', 'var_2', 'LMax']]
LMax_sol.columns = ['Method', 'Res', 'value']
LMax_sol['Error'] = 'LMax'

stat_sol = pd.DataFrame()
stat_sol = stat_sol.append(L2_sol)
stat_sol = stat_sol.append(LMax_sol)

sns.lineplot(x='Res', y='value',
             hue='Method', style='Error',
             data=stat_sol, markers=True,
            hue_order=['gradAlpha', 'fitParaboloid', 'RDF', 'heightFunction']).set(xscale='log', yscale='log', ylim=(1e-5, 5))


plt.ylabel('curvature error')
plt.xlabel('Resolution per Radius')
plt.savefig('advectCircle_hex.pdf')
stat_sol.to_csv("curv_error_hex.csv", index=False)
# tri
case = [['plicRDF'],  # [['isoSurface', 'plicRDF'],
        ['gradAlpha', 'RDF', 'fitParaboloid'],
        ['Grid1', 'Grid2', 'Grid3', 'Grid4', 'Grid5', 'Grid6', 'Grid7', 'Grid8', 'Grid9', 'Grid10']]


baseCase = 'CasesTri'
solutionDir = 'reconSurfaces'
file = 'K__freeSurf.raw'

sol = casefoam.posField_to_timeSeries(
    solutionDir, file, userFunction, case, baseCase, axis=3)
sol = sol.reset_index()

res = np.linspace(10, 128, 10)
res = res.astype(int)
sol = sol.replace('Grid1', res[0])
sol = sol.replace('Grid2', res[1])
sol = sol.replace('Grid3', res[2])
sol = sol.replace('Grid4', res[3])
sol = sol.replace('Grid5', res[4])
sol = sol.replace('Grid6', res[5])
sol = sol.replace('Grid7', res[6])
sol = sol.replace('Grid8', res[7])
sol = sol.replace('Grid9', res[8])
sol = sol.replace('Grid10', res[9])


L1_sol = sol.groupby(['var_1', 'var_2']).mean()
L1_sol = L1_sol.reset_index()
L1_sol = L1_sol[['var_1', 'var_2', 'L1']]
L1_sol.columns = ['Method', 'Res', 'value']
L1_sol['Error'] = 'L1'

L2_sol = sol.groupby(['var_1', 'var_2']).mean()
L2_sol = L2_sol.reset_index()
L2_sol = L2_sol[['var_1', 'var_2', 'L2']]
L2_sol.columns = ['Method', 'Res', 'value']
L2_sol['Error'] = 'L2'

LMax_sol = sol.groupby(['var_1', 'var_2']).max()
LMax_sol = LMax_sol.reset_index()
LMax_sol = LMax_sol[['var_1', 'var_2', 'LMax']]
LMax_sol.columns = ['Method', 'Res', 'value']
LMax_sol['Error'] = 'LMax'

stat_sol = pd.DataFrame()
stat_sol = stat_sol.append(L2_sol)
stat_sol = stat_sol.append(LMax_sol)

plt.figure()
sns.lineplot(x='Res', y='value',
             hue='Method', style='Error',
             data=stat_sol, markers=True,
            hue_order=['gradAlpha', 'fitParaboloid', 'RDF']).set(xscale='log', yscale='log', ylim=(5e-5, 20))
plt.ylabel('curvature error')
plt.xlabel('Resolution per Radius')
plt.savefig('advectCircle_tri.pdf')
stat_sol.to_csv("curv_error_tri.csv", index=False)

plt.show()
