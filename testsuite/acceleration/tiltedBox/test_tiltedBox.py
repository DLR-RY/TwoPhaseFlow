import os
from numpy import logspace
import pytest
import oftest
from oftest import run_reset_case

@pytest.fixture  # descruct all tests
def load_magU():
    dir_name = os.path.dirname(os.path.abspath(__file__))
    # dir_name = oftest.base_dir()
    magUFile = os.path.join(dir_name,"postProcessing", "maxU","0","fieldMinMax.dat")
    magU = oftest.read_functionObject(magUFile)
    magU.columns = ['time','minU','maxU']
    return magU

def test_tiltedBox(run_reset_case,load_magU):
    logs = oftest.path_logs()
    assert len(logs) > 0
    for log in logs:
        assert oftest.case_status(log) == 'completed' # checks if run completes
    maxU = load_magU.max()["maxU"]
    assert maxU < 1e-7
    


