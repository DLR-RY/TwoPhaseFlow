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
file = 'interfaceEnergyFluxLiquid_freeSurf.raw'
postFunction = postFunctions.getRadius

surfPos = casefoam.posField_to_timeSeries(solutionDir, file,postFunction, case, baseCase)

surfPos.columns = ['min', 'r(t)', 'max', 'interFaceType', 'Method', 'Resolution']

surfPos = surfPos.sort_index()
surfPos = surfPos.reset_index('time')
surfPos[['min', 'r(t)', 'max']] *= 1000

surfPos.to_csv("results/scriven3D.csv",index=False)
