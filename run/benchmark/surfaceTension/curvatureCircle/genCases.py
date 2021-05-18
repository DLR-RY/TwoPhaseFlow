import casefoam
import numpy as np
import math
# directory of the base case
baseCase = 'curvature2D'

# list of parent, child and grandchild names
caseStructure = [['gradAlpha', 'RDF','heightFunction','fitParaboloid'],
                 ['hex'],
                 ['Grid1','Grid2','Grid3','Grid4','Grid5','Grid6','Grid7','Grid8','Grid9','Grid10']]

# dictionarys with data for the caseData dictionary
update_gradAlpha = {
    'constant/transportProperties': {
        'surfaceForces': {'surfaceTensionForceModel': 'gradAlpha'}}}

update_RDF = {
    'constant/transportProperties': {
        'surfaceForces': {'surfaceTensionForceModel': 'RDF'}}}


update_heightFunction = {
    'constant/transportProperties': {
        'surfaceForces': {'surfaceTensionForceModel': 'heightFunction'}}}

update_fitParaboloid = {
    'constant/transportProperties': {
        'surfaceForces': {'surfaceTensionForceModel': 'fitParaboloid'}}}

update_hex = dict() #{'system/blockMeshDict': {}} # do nothing

# this function does the same as update_coarse etc but is more elegant
def changeBlockMesh(Nx):
    return {
    'system/blockMeshDict': {
            'blocks': ['hex',
                       [0, 1, 2, 3, 4, 5, 6, 7],
                       '(%s 1 %s)' % (Nx, Nx),
                       'simpleGrading',
                       '(1 1 1)']}}

res = np.logspace(math.log10(25),math.log10(600),10)
res = res.astype(int)
# dictionary of data to update
caseData = {'gradAlpha': update_gradAlpha,
            'RDF': update_RDF,
            'fitParaboloid': update_fitParaboloid,
            'heightFunction': update_heightFunction,
            'hex': update_hex,
            'Grid1': changeBlockMesh(res[0]),
            'Grid2': changeBlockMesh(res[1]),
            'Grid3': changeBlockMesh(res[2]),
            'Grid4': changeBlockMesh(res[3]),
            'Grid5': changeBlockMesh(res[4]),
            'Grid6': changeBlockMesh(res[5]),
            'Grid7': changeBlockMesh(res[6]),
            'Grid8': changeBlockMesh(res[7]),
            'Grid9': changeBlockMesh(res[8]),
            'Grid10': changeBlockMesh(res[9])
            }

# generate cases
casefoam.mkCases(baseCase,caseStructure, caseData, hierarchy='tree',writeDir='Cases')
