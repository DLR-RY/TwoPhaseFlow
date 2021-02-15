import casefoam

case = [['isoSurface', 'plicRDF'],
        ['implicitGrad', 'explicitGrad', 'Schrage'],
        ['grid1', 'grid2', 'grid3']]

update_isoSurface = {
    'system/fluid/fvSolution': {
        'solvers': {'alpha.water': {'advectionScheme': 'MULESScheme',
                                    'reconstructionScheme': 'isoSurface'}}}}

update_plicRDF = {
    'system/fluid/fvSolution': {
        'solvers': {'alpha.water': {'advectionScheme': 'isoAdvection',
                                    'reconstructionScheme': 'plicRDF'}}}}

update_implicitGrad = {
    'constant/fluid/phaseChangeProperties': {
         'energySourceTermModel': 'implicitGrad'}}
update_explicitGrad = {
    'constant/fluid/phaseChangeProperties': {
        'energySourceTermModel': 'selectedGradExplicit'}}
update_Schrage = {
    'constant/fluid/phaseChangeProperties': {
        'energySourceTermModel': 'Schrage'}}

def changeBlockMesh(Nx):
    return {
            'system/fluid/blockMeshDict': {
                'blocks': ['hex',
                          [0, 1, 2, 3, 4, 5, 6, 7],
                          '(%s %s %s)' % (Nx, Nx, Nx),
                          'simpleGrading',
                          '(1 1 1)']}
            }

data = {'isoSurface': update_isoSurface,
        'plicRDF': update_plicRDF,
        'implicitGrad': update_implicitGrad,
        'explicitGrad': update_explicitGrad,
        'Schrage': update_Schrage,
        'grid1': changeBlockMesh(50),
        'grid2': changeBlockMesh(100),
        'grid3': changeBlockMesh(200)}

casefoam.mkCases('scriven3D', case, data, 'tree', writeDir='Cases')
