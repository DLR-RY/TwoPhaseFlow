import os
import pytest
import oftest
from oftest import run_case, clean_case

@pytest.fixture(scope='class') # descruct all tests
def load_errorfiles():
    dir_name = os.path.dirname(os.path.abspath(__file__))
    # dir_name = oftest.base_dir()
    ferr_serial = os.path.join(dir_name,"volumeFractionError.dat")
    ferr_par = os.path.join(dir_name,"volumeFractionErrorPar.dat")
    err_serial = oftest.read_functionObject(ferr_serial)
    err_parrallel = oftest.read_functionObject(ferr_par)
    return err_serial, err_parrallel


class TestDiscUniFlow:

    def test_completed(self,run_case):
        log = oftest.path_log()
        assert oftest.case_status(log) == 'completed' # checks if run completes

    def test_parallel(self,run_case,load_errorfiles):
        err_serial, err_parrallel = load_errorfiles
        # exlude last two lines
        error = abs(err_parrallel.iloc[:,:-2]) - abs(err_serial.iloc[:,:-2])
        max_error = max(error.abs().max().values)

        assert max_error < 1e-14

    def test_accuracy(self,run_case,load_errorfiles):
        err_serial, err_parrallel = load_errorfiles

        max_err_shape =  abs(err_serial.iloc[:,1].max())
        max_err_mass =  abs(err_serial.iloc[:,2].max())
        max_err_bound =  abs(err_serial.iloc[:,3].max())
        assert  max_err_shape <= 0.0122
        assert  max_err_mass <= 1e-13
        assert  max_err_bound <= 1e-8

    def test_clean(clean_case):
        pass