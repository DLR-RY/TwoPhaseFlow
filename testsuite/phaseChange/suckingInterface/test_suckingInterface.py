import os
import pytest
import oftest
import pandas as pd
from oftest import run_case, run_reset_case
import casefoam
from casefoam import postFunctions
# from testsuite.surfaceTension.sinWave.sinwave_prosperetti import sinwave_prosperetti
# from .sinwave_prosperetti import sinwave_prosperetti
from scipy import optimize, integrate
import math
from scipy import interpolate

import numpy as np



class SuckingInterface():

    def __init__(self):
        dir_name = os.path.dirname(os.path.abspath(__file__))
        data = pd.read_csv(dir_name + "/analyticalResults.csv")
        self.posInterp = interpolate.interp1d(data['t'], data['s'])


    def x(self,t):
        """Interface position"""
        return self.posInterp(t)


@pytest.fixture()
def expected():
    dir_name = os.path.dirname(os.path.abspath(__file__))
    expected = pd.read_csv(os.path.join(dir_name, "expected_results.csv"))
    return expected



@pytest.fixture  # descruct all tests
def load_errorfile():
    
    solutionDir = 'surfaces'
    file = 'alpha.water_constantIso.raw'
    dir_name = os.path.dirname(os.path.abspath(__file__))
    ana = SuckingInterface()
    # dir_name = oftest.base_dir()
    postFunction = postFunctions.getFreeSurfaceWallAndCentre

    sol = casefoam.posField_to_timeSeries(
        solutionDir, file, postFunction,baseCase=dir_name, axis=0)
    sol = sol.reset_index()
    sol.columns = ['time', 'min', 'mean', 'max',
                'interfaceType']
    sol = sol.sort_values('time')
    sol['analytical'] = sol.apply(lambda x: abs(ana.x(x['time'])), axis=1)
    sol['error'] = (sol['max'] - sol['analytical'])**2

    int_err = integrate.trapz(sol['error'], x=sol['time'])
    err = 1/0.5*int_err
    return err


def simMod(PCM):
    dir_name = os.path.dirname(os.path.abspath(__file__))
    filemod = {
            "system/simulationParameter": [
            ("PCM", PCM)
        ]
    }
    # case_mod = oftest.Case_modifiers(filemod, dir_name)
    meta_data = {"PCM": PCM}
    case_mod = oftest.Case_modifiers(filemod, dir_name, meta_data)
    return case_mod


explicitGrad = simMod("selectedGradExplicit")
implicitGrad = simMod("implicitGrad")
Schrage = simMod("Schrage")

parameters = [
    explicitGrad,
    implicitGrad,
    Schrage
]

results = {
    "err": [],
    "PCM": [],
    "runtime": [],
}


@pytest.mark.parametrize(
    "run_reset_case",
    parameters,
    indirect=["run_reset_case"],
    ids=[
        "explicitGrad",
        "implicitGrad",
        "Schrage"
    ]
)
def test_SuckingInterface(run_reset_case, expected, load_errorfile):

    log = oftest.path_log()
    err = load_errorfile
    scheme = run_reset_case.meta_data["PCM"]

    assert oftest.case_status(log) == "completed"  # checks if run completes

    results["err"].append(err)
    results["runtime"].append(oftest.run_time(log)[0])
    results["PCM"].append(scheme)

    exp_res = oftest.expected_results([1],(scheme,))
    assert err == pytest.approx(exp_res["err"],rel=0.001)



def test_results():
    print("results", results)
    dir_name = os.path.dirname(os.path.abspath(__file__))
    res = pd.DataFrame(results)
    res.to_csv(os.path.join(dir_name, "results_suckingInterface.csv"),index=False)