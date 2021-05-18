import casefoam

case = [['plicRDF'], #[['isoSurface', 'plicRDF'],
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

def replaceString(Nx):
    return {
    'triSquare.geo': {'#!stringManipulation': {'replaceNx': '%s' % Nx,
                                              'replaceNz': '%s' % (3*Nx) }}}


data = {#'isoSurface': update_isoSurface,
        'plicRDF': update_plicRDF,
        'gradAlpha': update_gradAlpha,
        'RDF': update_RDF,
        'fitParaboloid': update_fitParaboloid,
        'coarse': replaceString(32),
        'mid': replaceString(64),
        'fine': replaceString(128)}

casefoam.mkCases('sinWaveTri', case, data, 'tree',writeDir='Cases')
