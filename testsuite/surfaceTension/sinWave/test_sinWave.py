import os
import pytest
import oftest
import pandas as pd
from oftest import run_case, run_reset_case
import casefoam
from casefoam import postFunctions
# from testsuite.surfaceTension.sinWave.sinwave_prosperetti import sinwave_prosperetti
# from .sinwave_prosperetti import sinwave_prosperetti
from scipy import integrate
import math
import numpy as np
from scipy.special import erfc


class sinwave_prosperetti:
    # Input
    def __init__(self, rho1, rho2, wavelength,H0,nu, sigma):
        self.rho1 = rho1  # Density liquid 1
        self.rho2 = rho2  # Density liquid 2
        self.wavelength = wavelength  # Wavelength
        self.H0 = H0  # Initial height
        self.nu = nu  # Kinematic viscosity
        self.sigma = sigma       # Surface tension

        # Constantes
        self.waveNb = (2*math.pi)/self.wavelength
        self.omega0 = abs(
            np.sqrt((self.sigma*self.waveNb**3)/(self.rho1+self.rho2)))
        self.epsilon = (self.nu*self.waveNb**2)/self.omega0
        self.beta = (self.rho1*self.rho2)/(self.rho1+self.rho2)**2

        self.p = [1,
                  -4*self.beta*(self.epsilon*self.omega0)**0.5,
                  2*(1-6*self.beta)*self.epsilon*self.omega0,
                  4*(1-3*self.beta)*(self.epsilon*self.omega0)**(3/2),
                  (1-4*self.beta)*(self.epsilon*self.omega0)**2+self.omega0**2]
        self.r = np.roots(self.p)

        self.Z = [0, 0, 0, 0]
        self.Z[0] = (self.r[1]-self.r[0]) * \
            (self.r[2]-self.r[0])*(self.r[3]-self.r[0])
        self.Z[1] = (self.r[0]-self.r[1]) * \
            (self.r[2]-self.r[1])*(self.r[3]-self.r[1])
        self.Z[2] = (self.r[0]-self.r[2]) * \
            (self.r[1]-self.r[2])*(self.r[3]-self.r[2])
        self.Z[3] = (self.r[0]-self.r[3]) * \
            (self.r[1]-self.r[3])*(self.r[2]-self.r[3])

    def f(self, t, i):
        return np.exp(((self.r[i]**2-self.epsilon*self.omega0)*t)/ \
            self.omega0)*erfc((self.r[i]*t**0.5)/(self.omega0**0.5))

    def a(self, t):
        A1 = (4*(1-4*self.beta)*self.epsilon**2) / \
            (8*(1-4*self.beta)*self.epsilon**2+1)*erfc(np.sqrt(self.epsilon*t))

        sum = 0
        for i in range(0, 4):
            sum += (self.r[i]*self.omega0**2)/(self.Z[i] * \
                (self.r[i]**2-self.epsilon*self.omega0))*self.f(t, i)
        return np.real(A1 + sum)

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
    rho1 = 1  # Density liquid 1
    rho2 = 1  # Density liquid 2
    wavelength = 0.003  # Wavelength
    H0 = 3e-5  # Initial height
    nu = 0.001  # Kinematic viscosity
    sigma = 1       # Surface tension
    omega0 = math.sqrt(sigma*(2*math.pi/wavelength)**3/(2*rho1))
    ana = sinwave_prosperetti(rho1, rho2, wavelength, H0, nu, sigma)
    # dir_name = oftest.base_dir()
    postFunction = postFunctions.getFreeSurfaceWallAndCentre

    sol = casefoam.posField_to_timeSeries(
        solutionDir, file, postFunction,baseCase=dir_name, axis=1)
    sol = sol.reset_index()
    sol.columns = ['time', 'min', 'mean', 'max',
                'interfaceType']
    sol['max'] /= 3e-5
    sol['time'] /= 1/omega0
    sol = sol.sort_values('time')
    sol['analytical'] = sol.apply(lambda x: abs(ana.a(x['time'])), axis=1)
    sol['error'] = H0*(sol['max'] - sol['analytical'])**2

    int_err = integrate.trapz(sol['error'], x=sol['time'])
    err = 1/wavelength*math.sqrt(1/25*int_err)
    return err


def simMod(res, STM):
    dir_name = os.path.dirname(os.path.abspath(__file__))
    filemod = {
            "system/simulationParameter": [
            ("STM", STM),
            ("nx", res),
            ("ny", 3*res),
            ("nz", 1),
        ]
    }
    # case_mod = oftest.Case_modifiers(filemod, dir_name)
    meta_data = {"STM": STM, "Res": res}
    case_mod = oftest.Case_modifiers(filemod, dir_name, meta_data)
    return case_mod


c_gradAlpha32 = simMod(32, "gradAlpha")
c_gradAlpha64 = simMod(64, "gradAlpha")
c_gradAlpha128 = simMod(128, "gradAlpha")

c_RDF32 = simMod(32, "RDF")
c_RDF64 = simMod(64, "RDF")
c_RDF128 = simMod(128, "RDF")

c_fitParaboloid32 = simMod(32, "fitParaboloid")
c_fitParaboloid64 = simMod(64, "fitParaboloid")
c_fitParaboloid128 = simMod(128, "fitParaboloid")

parameters = [
    c_gradAlpha32,
    c_gradAlpha64,
    pytest.param(c_gradAlpha128, marks=pytest.mark.slow),
    c_RDF32,
    c_RDF64,
    pytest.param(c_RDF128, marks=pytest.mark.slow),
    c_fitParaboloid32,
    c_fitParaboloid64,
    pytest.param(c_fitParaboloid128, marks=pytest.mark.slow),
]

results = {
    "err": [],
    "STM": [],
    "Res": [],
    "runtime": [],
}


@pytest.mark.parametrize(
    "run_reset_case",
    parameters,
    indirect=["run_reset_case"],
    ids=[
        "gradAlpha32",
        "gradAlpha64",
        "gradAlpha128",
        "RDF32",
        "RDF64",
        "RDF128",
        "fitParaboloid32",
        "fitParaboloid64",
        "fitParaboloid128"
    ]
)
def test_sinWave(run_reset_case, expected, load_errorfile):

    log = oftest.path_log()
    err = load_errorfile
    scheme = run_reset_case.meta_data["STM"]
    res = run_reset_case.meta_data["Res"]


    results["err"].append(err)
    results["runtime"].append(oftest.run_time(log)[0])
    results["STM"].append(scheme)
    results["Res"].append(res)

    assert oftest.case_status(log) == "completed"  # checks if run completes


def test_results():
    print("results", results)
    dir_name = os.path.dirname(os.path.abspath(__file__))
    res = pd.DataFrame(results)
    res.to_csv(os.path.join(dir_name, "results_sinWave.csv"),index=False)