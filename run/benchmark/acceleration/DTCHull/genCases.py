import casefoam

# directory of the base case
baseCase = "DTCHull"

# list of parent, child and grandchild names
caseStructure = [
    [
        "InterFoam_gravity_LocalEuler",
        "InterFoam_gravity_Euler",
        "InterFlow_gravity_LocalEuler",
        "InterFlow_gravity_Euler",
        "InterFlow_improvedGravity_LocalEuler",
        "InterFlow_improvedGravity_Euler",
    ]
]

# this function does the same as update_coarse etc but is more elegant
# DDTSCHEME         localEuler;
# CFL               10;
# AlphaCFL          5;
# SOLVER            interFoam;
# ENDTIME           4000;
# DELTAT            1;
# WRITEINTERVAL     100;
# GRAVITYMODEL      improvedGravity;
# NALPHASUBCYCLES   5;

def parameters(
    ddtScheme, cfl, alphaCFL, solver, endTime, deltaT, writeInterval, gravityModel, nAlphaSubCycles
):
    return {
        "system/simulationParameters": {
            "DDTSCHEME": ddtScheme,
            "CFL": cfl,
            "AlphaCFL": alphaCFL,
            "SOLVER": solver,
            "ENDTIME": endTime,
            "DELTAT": deltaT,
            "WRITEINTERVAL": writeInterval,
            "GRAVITYMODEL": gravityModel,
            "NALPHASUBCYCLES": nAlphaSubCycles
        }
    }


# dictionary of data to update
caseData = {
    "InterFoam_gravity_LocalEuler": parameters(
        "localEuler", 0.5, 0.5, "interFoam", 4000, 1, 500, "gravity",1
    ),
    "InterFoam_gravity_Euler": parameters(
        "Euler", 0.5, 0.5, "interFoam", 200, 1e-3, 500, "gravity",1 
    ),
    "InterFlow_gravity_LocalEuler": parameters(
        "localEuler", 0.5, 0.5, "interFlow", 4000, 1, 500, "gravity", 5
    ),
    "InterFlow_gravity_Euler": parameters(
        "Euler", 0.5, 0.5, "interFlow", 200, 1e-3, 500, "gravity", 5
    ),
    "InterFlow_improvedGravity_LocalEuler": parameters(
        "localEuler", 0.5, 0.5, "interFlow", 4000, 1, 500, "improvedGravity",5 
    ),
    "InterFlow_improvedGravity_Euler": parameters(
        "Euler", 0.5, 0.5, "interFlow", 200, 1e-3, 500, "improvedGravity",5 
    ),
}

# generate cases
casefoam.mkCases(baseCase, caseStructure, caseData, hierarchy="tree", writeDir="Cases")
