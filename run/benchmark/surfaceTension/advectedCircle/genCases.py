import casefoam
import numpy as np

case = [['plicRDF'], #[['isoSurface', 'plicRDF'],
        ['gradAlpha','RDF','heightFunction','fitParaboloid'],
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

update_heightFunction = {
    'constant/transportProperties': {
        'surfaceForces': {'surfaceTensionForceModel': 'heightFunction'}}}

def changeBlockMesh(Nx):
    return {
    'system/blockMeshDict': {
            'blocks': ['hex',
                       [0, 1, 2, 3, 4, 5, 6, 7],
                       '(%s 1 %s)' % (4*Nx, Nx),
                       'simpleGrading',
                       '(1 1 1)']}}

res = np.linspace(10,128,10)
res = res.astype(int)

data = {#'isoSurface': update_isoSurface,
        'plicRDF': update_plicRDF,
        'gradAlpha': update_gradAlpha,
        'RDF': update_RDF,
        'heightFunction': update_heightFunction,
        'fitParaboloid': update_fitParaboloid,
        'Grid1': changeBlockMesh(res[0]),
        'Grid2': changeBlockMesh(res[1]),
        'Grid3': changeBlockMesh(res[2]),
        'Grid4': changeBlockMesh(res[3]),
        'Grid5': changeBlockMesh(res[4]),
        'Grid6': changeBlockMesh(res[5]),
        'Grid7': changeBlockMesh(res[6]),
        'Grid8': changeBlockMesh(res[7]),
        'Grid9': changeBlockMesh(res[8]),
        'Grid10': changeBlockMesh(res[9])}

casefoam.mkCases('advectedCircle', case, data, 'tree',writeDir='Cases')
