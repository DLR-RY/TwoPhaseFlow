import casefoam

# directory of the base case
baseCase = 'tiltedBox'

# list of parent, child and grandchild names
caseStructure = [['hex'],
                 ['Grid1','Grid2','Grid3']]


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


# dictionary of data to update
caseData = {'hex': update_hex,
            'Grid1': changeBlockMesh(32),
            'Grid2': changeBlockMesh(64),
            'Grid3': changeBlockMesh(128),
            }

# generate cases
casefoam.mkCases(baseCase,caseStructure, caseData, hierarchy='tree',writeDir='Cases')
