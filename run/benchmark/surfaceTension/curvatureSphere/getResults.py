# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
import seaborn as sns
#output_notebook()
sns.set_style("ticks")

caseStructure = [['gradAlpha', 'RDF','fitParaboloid','heightFunction'],
                 ['hex'],
                 ['Grid1','Grid2','Grid3','Grid4','Grid5','Grid6','Grid7','Grid8','Grid9','Grid10']]


baseCase = 'Cases'
solutionDir = 'reconstructionError/0'
file = 'error.dat'

curv = casefoam.time_series(solutionDir, file,caseStructure, baseCase)
curv = curv[[5,6,9,'var_0','var_2']]
curv.columns = ['Curv_max','Curv_mean','Length','Method','Res']
curv['DeltaX'] = 0.75/curv['Length']
curv2 =  pd.melt(curv,id_vars=curv.columns[2:],value_vars=curv.columns[:2],value_name='curv',var_name='error')

sns.lineplot(x='DeltaX', y='curv', hue='Method',style='error', data=curv2,markers=True,
            hue_order=['gradAlpha', 'fitParaboloid', 'RDF', 'heightFunction']).set(xscale='log',yscale='log',xlim=(None, 120))
plt.ylabel('curvature error')
plt.xlabel('Resolution per Radius')
plt.savefig('curv3d_hex.pdf')
curv2.to_csv("curv3d_hex.csv",index=False)

caseStructure = [['gradAlpha', 'RDF','fitParaboloid'],
                 ['tri'],
                 ['Grid1','Grid2','Grid3','Grid4','Grid5','Grid6','Grid7','Grid8','Grid9','Grid10']]


baseCase = 'CasesTri'
solutionDir = 'reconstructionError/0'
file = 'error.dat'

curvTri = casefoam.time_series(solutionDir, file,caseStructure, baseCase)
curvTri = curvTri[[5,6,9,'var_0','var_2']]
curvTri.columns = ['Curv_max','Curv_mean','Length','Method','Res']
curvTri['DeltaX'] = 0.75/curvTri['Length']
curvTri2 =  pd.melt(curvTri,id_vars=curvTri.columns[2:],value_vars=curvTri.columns[:2],value_name='curv',var_name='error')


plt.figure()
sns.lineplot(x='DeltaX', y='curv', hue='Method',style='error', data=curvTri2,markers=True,
            hue_order=['gradAlpha', 'fitParaboloid', 'RDF']).set(xscale='log',yscale='log',xlim=(13, None))
plt.ylabel('curvature error')
plt.xlabel('Resolution per Radius')
plt.savefig('curv3d_tri.pdf')
curvTri2.to_csv("curv3d_tri.csv",index=False)

plt.show()