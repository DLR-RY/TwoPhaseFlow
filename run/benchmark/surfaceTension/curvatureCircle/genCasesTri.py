import casefoam
import numpy as np
import math
# directory of the base case
baseCase = 'curvature2DTri'

# list of parent, child and grandchild names
caseStructure = [['gradAlpha', 'RDF','fitParaboloid'],
                 ['tri'],
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

update_tri = dict() #{'system/blockMeshDict': {}} # do nothing

# this function does the same as update_coarse etc but is more elegant

def replaceString(Nx):
    return {
    'triSquare.geo': {'#!stringManipulation': {'nx=10': f"nx={Nx}"}}}

res = np.logspace(math.log10(25),math.log10(600),10)
res = res.astype(int)
# dictionary of data to update
caseData = {'gradAlpha': update_gradAlpha,
            'RDF': update_RDF,
            'fitParaboloid': update_fitParaboloid,
            'heightFunction': update_heightFunction,
            'tri': update_tri,
            'Grid1': replaceString(res[0]),
            'Grid2': replaceString(res[1]),
            'Grid3': replaceString(res[2]),
            'Grid4': replaceString(res[3]),
            'Grid5': replaceString(res[4]),
            'Grid6': replaceString(res[5]),
            'Grid7': replaceString(res[6]),
            'Grid8': replaceString(res[7]),
            'Grid9': replaceString(res[8]),
            'Grid10': replaceString(res[9])
            }

# generate cases
casefoam.mkCases(baseCase,caseStructure, caseData, hierarchy='tree',writeDir='CasesTri')
