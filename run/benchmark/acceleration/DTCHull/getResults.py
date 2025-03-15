import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import casefoam
import numpy as np



caseStructure = [['IFLocalEuler','IFEuler']]

baseCase = 'Cases'


solutionDir = 'forces/0'
file = 'force.dat'

forces = casefoam.time_series(solutionDir,file,caseStructure,baseCase)

forces = forces.reset_index()
forces = forces.drop(columns=[4,5,6,7,8,9])
forces.columns = ['t','Fx', 'Fy','Fz','Case']
forces = forces.groupby('Case').apply(lambda group: group.iloc[10:], include_groups=False).reset_index() # remove the first 10 timesteps
forces = forces
forces['index'] = forces.groupby('Case').cumcount()
# sns.lineplot(data=forces, x='index', y='Fx', hue='Case')  # Plot over the index
print(forces)
sns.lineplot(data=forces, x='index', y='Fx', hue='Case')
plt.show()
# sol = sol.drop(columns='minU')
# table = sol[sol['t'] == 0.01]
# table = table.drop(columns='t')

# table = table.set_index(['GridType','Resolution'])

# table.to_html("results.html")
# table.to_csv("results.csv")