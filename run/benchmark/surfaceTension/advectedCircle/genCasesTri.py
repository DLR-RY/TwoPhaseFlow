import casefoam
import numpy as np

case = [['plicRDF'], #[['isoSurface', 'plicRDF'],
        ['gradAlpha','RDF','fitParaboloid'],
        ['Grid1','Grid2','Grid3','Grid4','Grid5','Grid6','Grid7','Grid8','Grid9','Grid10']]

update_isoSurface = {
    'system/fvSolution': {
        'solvers': {'alpha.water': {'advectionScheme': 'MULESScheme',
                                    'reconstructionScheme': 'isoSurface'}}}}

update_plicRDF = {
    'system/fvSolution': {
        'solvers': {'alpha.water': {'advectionScheme': 'isoAdvection',
                                    'reconstructionScheme': 'plicRDF'}}}}

update_gradAlpha = {
    'constant/transportProperties': {
        'surfaceForces': {'surfaceTensionForceModel': 'gradAlpha'}}}

update_RDF = {
    'constant/transportProperties': {
        'surfaceForces': {'surfaceTensionForceModel': 'RDF'}}}

update_fitParaboloid = {
    'constant/transportProperties': {
        'surfaceForces': {'surfaceTensionForceModel': 'fitParaboloid'}}}

def replaceString(Nx):
    return {
    'triSquare.geo': {'#!stringManipulation': {'replaceNx': '%s' % (4*Nx),
                                              'replaceNz': '%s' % (Nx) }}}

res = np.linspace(10,128,10)
res = res.astype(int)

data = {#'isoSurface': update_isoSurface,
        'plicRDF': update_plicRDF,
        'gradAlpha': update_gradAlpha,
        'RDF': update_RDF,
        'fitParaboloid': update_fitParaboloid,
        'Grid1': replaceString(res[0]),
        'Grid2': replaceString(res[1]),
        'Grid3': replaceString(res[2]),
        'Grid4': replaceString(res[3]),
        'Grid5': replaceString(res[4]),
        'Grid6': replaceString(res[5]),
        'Grid7': replaceString(res[6]),
        'Grid8': replaceString(res[7]),
        'Grid9': replaceString(res[8]),
        'Grid10': replaceString(res[9])}

casefoam.mkCases('advectedCircleTri', case, data, 'tree',writeDir='CasesTri')
