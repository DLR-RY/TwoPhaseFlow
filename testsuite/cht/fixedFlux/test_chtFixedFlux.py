import os
import pytest
import oftest
import pandas as pd
from oftest import run_case, run_reset_case
import casefoam

@pytest.fixture  
def load_errorfile():
    dir_name = os.path.dirname(os.path.abspath(__file__))

    heatFlux = casefoam.time_series(
        'fluid/energyFluxes/0', 'interfaceEnergyFluxes.dat',baseCase=dir_name)

    massFlux = casefoam.time_series(
        'fluid/massFlux/0', 'volFieldValue.dat',baseCase=dir_name)
    return heatFlux, massFlux

def test_chtFixedHeatFlux(run_reset_case, load_errorfile):

    log = oftest.path_log()
    assert oftest.case_status(log) == "completed"  # checks if run completes
    
    heatFlux, massFlux = load_errorfile

    print(heatFlux.tail(1))
    print(massFlux.tail(1))
    q = heatFlux.tail(1)[1]
    mdot = massFlux.tail(1)[1]

    assert pytest.approx(q, rel=1e-4) == 1e-6
    assert pytest.approx(mdot, rel=1e-4) == 1e-12


