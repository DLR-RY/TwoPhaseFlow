import casefoam

# directory of the base case
baseCase = 'DTCHull'

# list of parent, child and grandchild names
caseStructure = [['IFLocalEuler','IFEuler']]

# this function does the same as update_coarse etc but is more elegant
# DDTSCHEME         localEuler;
# CFL               10;
# AlphaCFL          5;  
# SOLVER            interFoam;
# ENDTIME           4000;
# DELTAT            1;
# WRITEINTERVAL     100;

def parameters(ddtScheme, cfl, alphaCFL, solver, endTime, deltaT, writeInterval):
    return {
        'system/simulationParameters': {
            'DDTSCHEME': ddtScheme,
            'CFL': cfl,
            'AlphaCFL': alphaCFL,
            'SOLVER': solver,
            'ENDTIME': endTime,
            'DELTAT': deltaT,
            'WRITEINTERVAL': writeInterval
        }
    }


# dictionary of data to update
caseData = {'IFLocalEuler': parameters('localEuler',5.0,2.5,'interFoam', 4000, 1, 500),
            'IFEuler': parameters('Euler',5.0,2.5,'interFoam', 200, 1e-3, 500)
            }

# generate cases
casefoam.mkCases(baseCase,caseStructure, caseData, hierarchy='tree',writeDir='Cases')
