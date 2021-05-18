import casefoam

case = [['isoSurface', 'plicRDF'],
        ['gradAlpha','RDF','fitParaboloid'],
        ['coarse', 'mid', 'fine']]

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

def updateBlockMesh(x):
    return  {'system/blockMeshDict': {
            'blocks': ['hex',
                   [0, 1, 2, 3, 4, 5, 6, 7],
                   '(%s %s 1)' % (x , 3*x),
                   'simpleGrading',
                   '(1 1 1)']}}


data = {'isoSurface': update_isoSurface,
        'plicRDF': update_plicRDF,
        'gradAlpha': update_gradAlpha,
        'RDF': update_RDF,
        'fitParaboloid': update_fitParaboloid,
        'coarse': updateBlockMesh(32),
        'mid': updateBlockMesh(64),
        'fine': updateBlockMesh(128)}

casefoam.mkCases('sinWave', case, data, 'tree',writeDir='Cases')
