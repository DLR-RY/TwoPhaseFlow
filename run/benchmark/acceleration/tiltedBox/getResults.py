import pandas as pd
import matplotlib.pyplot as plt
import casefoam
import numpy as np



caseStructure = [['hex'],
                 ['Grid1','Grid2','Grid3']]

baseCase = 'Cases'


solutionDir = 'maxU/0'
file = 'fieldMinMax.dat'

sol = casefoam.time_series(solutionDir,file,caseStructure,baseCase)

sol.reset_index(inplace=True)
sol.columns = ['t','minU', 'maxU','GridType','Resolution']
sol = sol.drop(columns='minU')
table = sol[sol['t'] == 0.01]
table = table.drop(columns='t')

table = table.set_index(['GridType','Resolution'])

table.to_html("results.html")
table.to_csv("results.csv")