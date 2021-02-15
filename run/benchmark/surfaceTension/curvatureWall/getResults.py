# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
import seaborn as sns
sns.set_style("ticks")

caseStructure = [['gradAlpha', 'RDF', 'fitParaboloid'],
                 ['15', '30', '45', '60', '75'],
                 ['Grid1', 'Grid2', 'Grid3', 'Grid4', 'Grid5']]


baseCase = 'Cases'
solutionDir = 'reconstructionError/0'
file = 'error.dat'

curv = casefoam.time_series(solutionDir, file, caseStructure, baseCase)
curv = curv[[5, 6, 9, 'var_0', 'var_1', 'var_2']]
curv.columns = ['Curv_max', 'Curv_mean',
                'Length', 'Method', 'ContactAngle', 'Res']
curv['DeltaX'] = 0.5/curv['Length']
curv = curv.replace('15', 'alpha=15')
curv = curv.replace('30', 'alpha=30')
curv = curv.replace('45', 'alpha=45')
curv = curv.replace('60', 'alpha=60')
curv = curv.replace('75', 'alpha=75')

curv_gradAlpha = curv[curv.Method == 'gradAlpha'].copy()
curv_gradAlpha = pd.melt(
    curv_gradAlpha, id_vars=curv.columns[2:], value_vars=curv.columns[:2], value_name='curv', var_name='error')
curv_RDF = curv[curv.Method == 'RDF'].copy()
curv_RDF = pd.melt(
    curv_RDF, id_vars=curv.columns[2:], value_vars=curv.columns[:2], value_name='curv', var_name='error')
curv_fitParaboloid = curv[curv.Method == 'fitParaboloid'].copy()
curv_fitParaboloid = pd.melt(
    curv_fitParaboloid, id_vars=curv.columns[2:], value_vars=curv.columns[:2], value_name='curv', var_name='error')

ax = sns.lineplot(x='DeltaX', y='curv', hue='ContactAngle', style='error',
                  data=curv_gradAlpha, markers=True).set(xscale='log', yscale='log')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.ylabel('curvature error')
plt.xlabel('Resolution per Radius')
plt.savefig('curv2d_wall_gradAlpha_hex.pdf')
curv_gradAlpha.to_csv("curv2d_wall_gradAlpha_hex.csv", index=False)

plt.figure()
ax = sns.lineplot(x='DeltaX', y='curv', hue='ContactAngle', style='error',
                  data=curv_RDF, markers=True).set(xscale='log', yscale='log')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.ylabel('curvature error')
plt.xlabel('Resolution per Radius')
plt.savefig('curv2d_wall_RDF_hex.pdf')
curv_RDF.to_csv("curv2d_wall_RDF_hex.csv", index=False)

plt.figure()
ax = sns.lineplot(x='DeltaX', y='curv', hue='ContactAngle', style='error',
                  data=curv_fitParaboloid, markers=True).set(xscale='log', yscale='log')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.ylabel('curvature error')
plt.xlabel('Resolution per Radius')
plt.savefig('curv2d_wall_fitParaboloid_hex.pdf')
curv_fitParaboloid.to_csv("curv2d_wall_fitParaboloid_hex.csv", index=False)

caseStructureTri = [['gradAlpha', 'RDF', 'fitParaboloid'],
                    ['15', '30', '45', '60', '75'],
                    ['Grid1', 'Grid2', 'Grid3', 'Grid4', 'Grid5']]


baseCase = 'CasesTri'
solutionDir = 'reconstructionError/0'
file = 'error.dat'

curv = casefoam.time_series(solutionDir, file, caseStructureTri, baseCase)
curv = curv[[5, 6, 9, 'var_0', 'var_1', 'var_2']]
curv.columns = ['Curv_max', 'Curv_mean',
                'Length', 'Method', 'ContactAngle', 'Res']
curv['DeltaX'] = 0.5/curv['Length']
curv = curv.replace('15', 'alpha=15')
curv = curv.replace('30', 'alpha=30')
curv = curv.replace('45', 'alpha=45')
curv = curv.replace('60', 'alpha=60')
curv = curv.replace('75', 'alpha=75')

curv_gradAlpha = curv[curv.Method == 'gradAlpha'].copy()
curv_gradAlpha = pd.melt(
    curv_gradAlpha, id_vars=curv.columns[2:], value_vars=curv.columns[:2], value_name='curv', var_name='error')
curv_RDF = curv[curv.Method == 'RDF'].copy()
curv_RDF = pd.melt(
    curv_RDF, id_vars=curv.columns[2:], value_vars=curv.columns[:2], value_name='curv', var_name='error')
curv_fitParaboloid = curv[curv.Method == 'fitParaboloid'].copy()
curv_fitParaboloid = pd.melt(
    curv_fitParaboloid, id_vars=curv.columns[2:], value_vars=curv.columns[:2], value_name='curv', var_name='error')


plt.figure()
ax = sns.lineplot(x='DeltaX', y='curv', hue='ContactAngle', style='error',
                  data=curv_gradAlpha, markers=True).set(xscale='log', yscale='log')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.ylabel('curvature error')
plt.xlabel('Resolution per Radius')
plt.savefig('curv2d_wall_gradAlpha_tri.pdf')
curv_gradAlpha.to_csv("curv2d_wall_gradAlpha_tri.csv", index=False)


plt.figure()
ax = sns.lineplot(x='DeltaX', y='curv', hue='ContactAngle', style='error',
                  data=curv_RDF, markers=True).set(xscale='log', yscale='log')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.ylabel('curvature error')
plt.xlabel('Resolution per Radius')
plt.savefig('curv2d_wall_RDF_tri.pdf')
curv_RDF.to_csv("curv2d_wall_RDF_tri.csv", index=False)

plt.figure()
ax = sns.lineplot(x='DeltaX', y='curv', hue='ContactAngle', style='error',
                  data=curv_fitParaboloid, markers=True).set(xscale='log', yscale='log')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.ylabel('curvature error')
plt.xlabel('Resolution per Radius')
plt.savefig('curv2d_wall_fitParaboloid_tri.pdf')
curv_fitParaboloid.to_csv("curv2d_wall_fitParaboloid_tri.csv", index=False)
plt.show()
