# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
import seaborn as sns
# output_notebook()
sns.set_style("ticks")

caseStructure = [['gradAlpha', 'RDF', 'heightFunction', 'fitParaboloid'],
                 ['hex'],  # ,'tri'],
                 ['Grid1', 'Grid2', 'Grid3', 'Grid4', 'Grid5', 'Grid6', 'Grid7', 'Grid8', 'Grid9', 'Grid10']]


baseCase = 'Cases'  # 'Cases/isoSurface/explicitGrad/grid1/'
solutionDir = 'reconstructionError/0'
file = 'error.dat'

curv = casefoam.time_series(solutionDir, file, caseStructure, baseCase)
curv = curv[[5, 6, 9, 'var_0', 'var_1', 'var_2']]
curv.columns = ['Curv_max', 'Curv_mean', 'Length', 'Method', 'gridType', 'Res']
curv['DeltaX'] = 0.5/curv['Length']
curv = pd.melt(
    curv, id_vars=curv.columns[2:], value_vars=curv.columns[:2], value_name='curv', var_name='error')


ax = sns.lineplot(x='DeltaX', y='curv', hue='Method', style='error',
                data=curv, markers=True,
                hue_order=['gradAlpha', 'fitParaboloid', 'RDF', 'heightFunction']).set(xscale='log', yscale='log', ylim=(1e-5, 10))
# ax.set()
plt.ylabel('curvature error')
plt.xlabel('Resolution per Radius')
plt.savefig('curv2d_hex.pdf')
curv.to_csv("curv2d_hex.csv", index=False)

caseStructureTri = [['gradAlpha', 'RDF', 'heightFunction', 'fitParaboloid'],
                    ['tri'],  # ,'tri'],
                    ['Grid1', 'Grid2', 'Grid3', 'Grid4', 'Grid5', 'Grid6', 'Grid7', 'Grid8', 'Grid9', 'Grid10']]
plt.figure()
baseCase = 'CasesTri'  # 'Cases/isoSurface/explicitGrad/grid1/'
solutionDir = 'reconstructionError/0'
file = 'error.dat'


curvTri = casefoam.time_series(solutionDir, file, caseStructureTri, baseCase)
curvTri = curvTri[[5, 6, 9, 'var_0', 'var_1', 'var_2']]
curvTri.columns = ['Curv_max', 'Curv_mean',
                   'Length', 'Method', 'gridType', 'Res']
curvTri['DeltaX'] = 0.5/curvTri['Length']
curvTri = pd.melt(
    curvTri, id_vars=curvTri.columns[2:], value_vars=curvTri.columns[:2], value_name='curv', var_name='error')


sns.lineplot(x='DeltaX', y='curv', hue='Method', style='error', data=curvTri,
             markers=True,
            hue_order=['gradAlpha', 'fitParaboloid', 'RDF']).set(xscale='log', yscale='log', ylim=(None, 1000))
plt.ylabel('curvature error')
plt.xlabel('Resolution per Radius')
plt.savefig('curv2d_tri.pdf')
curvTri.to_csv("curv2d_tri.csv", index=False)
plt.show()
