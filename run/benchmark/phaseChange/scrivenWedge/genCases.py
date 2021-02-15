import casefoam

case = [[ 'plicRDF','isoSurface'],
        ['implicitGrad', 'explicitGrad', 'interfaceRes'],
        ['grid1', 'grid2', 'grid3', 'grid4']]

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
update_interfaceRes = {
    'constant/fluid/phaseChangeProperties': {
        'energySourceTermModel': 'Schrage'}}

def changeBlockMesh(Nx):
    return {
            'system/fluid/blockMeshDict': {
                'blocks': ['hex',
                          [0, 4, 1, 0, 3, 5, 2, 3],
                          '(%s 1 %s)' % (Nx, Nx),
                          'simpleGrading',
                          '(1 1 1)']}
            }

data = {'isoSurface': update_isoSurface,
        'plicRDF': update_plicRDF,
        'implicitGrad': update_implicitGrad,
        'explicitGrad': update_explicitGrad,
        'interfaceRes': update_interfaceRes,
        'grid1': changeBlockMesh(50),
        'grid2': changeBlockMesh(100),
        'grid3': changeBlockMesh(200),
        'grid4': changeBlockMesh(400)}

casefoam.mkCases('scrivenWedge', case, data, 'tree', writeDir='Cases')
