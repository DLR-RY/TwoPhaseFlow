# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
from casefoam import postFunctions
import seaborn as sns
#output_notebook()
sns.set_style("ticks")

case = [['isoSurface', 'plicRDF'],
        ['implicitGrad', 'explicitGrad', 'Schrage'],
        ['grid1', 'grid2', 'grid3']]


baseCase = 'Cases'
solutionDir = 'surfaces'
file = 'alpha.water_constantIso.raw'
postFunction = postFunctions.getFreeSurfaceWallAndCentre
surfPos = casefoam.posField_to_timeSeries(solutionDir, file,postFunction, case, baseCase)

surfPos.columns = ['min', 'x(t)', 'max', 'interFaceType', 'Method', 'Resolution']

surfPos = surfPos.sort_index()
surfPos = surfPos.reset_index('time')
surfPos[['min', 'x(t)', 'max']] *= 1000
surfPos = surfPos.drop(columns=['min',  'max'])

surfPos.to_csv("results/stefanProblem.csv",index=False)
