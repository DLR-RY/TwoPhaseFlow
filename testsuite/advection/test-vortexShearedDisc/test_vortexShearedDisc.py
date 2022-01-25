import os
import pytest
import oftest
import pandas as pd
from oftest import run_case, run_reset_case


@pytest.fixture(scope="class")  # descruct all tests
def load_errorfiles():
    dir_name = os.path.dirname(os.path.abspath(__file__))
    # dir_name = oftest.base_dir()
    ferr_serial = os.path.join(dir_name, "volumeFractionError.dat")
    ferr_par = os.path.join(dir_name, "volumeFractionErrorPar.dat")
    err_serial = oftest.read_functionObject(ferr_serial)
    err_parrallel = oftest.read_functionObject(ferr_par)
    return err_serial, err_parrallel

@pytest.fixture()
def expected():
    dir_name = os.path.dirname(os.path.abspath(__file__))
    expected = pd.read_csv(os.path.join(dir_name, "expected_results.csv"))
    return expected



@pytest.fixture  # descruct all tests
def load_errorfile():
    dir_name = os.path.dirname(os.path.abspath(__file__))
    # dir_name = oftest.base_dir()
    ferr = os.path.join(dir_name, "volumeFractionError.dat")
    err = oftest.read_functionObject(ferr)
    return err


class TestVortexShearedDisc:
    def test_completed(self, run_reset_case):
        log = oftest.path_log()
        assert oftest.case_status(log) == "completed"  # checks if run completes

    def test_parallel(self, run_reset_case, load_errorfiles):
        err_serial, err_parrallel = load_errorfiles
        error = abs(err_parrallel.iloc[:, 0:4]) - abs(err_serial.iloc[:, 0:4])
        max_error = max(error.abs().max().values)

        assert max_error < 1e-8

    def test_accuracy(self, run_reset_case, load_errorfiles):
        err_serial, err_parrallel = load_errorfiles
        err = err_serial[err_serial[0] == 8]

        max_err_shape = abs(err.iloc[:, 1].max())
        max_err_mass = abs(err.iloc[:, 2].max())
        max_err_bound = abs(err.iloc[:, 3].max())
        assert max_err_shape <= 2.35e-03*1.01
        assert max_err_mass <= 1e-13
        assert max_err_bound <= 5e-5


def simMod(res, reconScheme):
    dir_name = os.path.dirname(os.path.abspath(__file__))
    filemod = {
        "system/simulationParameter": [
            ("RECONSCHEME", reconScheme),
            ("nx", res),
            ("ny", 1),
            ("nz", res),
        ]
    }
    case_mod = oftest.Case_modifiers(filemod, dir_name)
    meta_data = {"reconScheme": reconScheme, "Res": res, "script": "Allrun.paraStudy"}
    case_mod = oftest.Case_modifiers(filemod, dir_name, meta_data)
    return case_mod


c_isoAlpha32 = simMod(32, "isoAlpha")
c_isoAlpha64 = simMod(64, "isoAlpha")
c_isoAlpha128 = simMod(128, "isoAlpha")
c_isoAlpha256 = simMod(256, "isoAlpha")
c_isoAlpha512 = simMod(512, "isoAlpha")

c_plicRDF32 = simMod(32, "plicRDF")
c_plicRDF64 = simMod(64, "plicRDF")
c_plicRDF128 = simMod(128, "plicRDF")
c_plicRDF256 = simMod(256, "plicRDF")
c_plicRDF512 = simMod(512, "plicRDF")

parameters = [
    c_isoAlpha32,
    c_isoAlpha64,
    c_isoAlpha128,
    c_isoAlpha256,
    pytest.param(c_isoAlpha512, marks=pytest.mark.slow),
    c_plicRDF32,
    c_plicRDF64,
    c_plicRDF128,
    c_plicRDF256,
    pytest.param(c_plicRDF512, marks=pytest.mark.slow),
]

results = {
    "err_shape": [],
    "err_mass": [],
    "err_bound": [],
    "reconScheme": [],
    "Res": [],
    "runtime": [],
}


@pytest.mark.parametrize(
    "run_reset_case",
    parameters,
    indirect=["run_reset_case"],
    ids=[
        "isoAlpha32",
        "isoAlpha64",
        "isoAlpha128",
        "isoAlpha256",
        "isoAlpha512",
        "plicRDF32",
        "plicRDF64",
        "plicRDF128",
        "plicRDF256",
        "plicRDF512"
    ]
)
def test_vortexShearedDisc(run_reset_case, expected, load_errorfile):

    log = oftest.path_log()
    assert oftest.case_status(log) == "completed"  # checks if run completes
    err = load_errorfile
    max_err_shape = abs(err.iloc[1, 1].max())
    max_err_mass = abs(err.iloc[:, 2].max())
    max_err_bound = abs(err.iloc[:, 3].max())
    scheme = run_reset_case.meta_data["reconScheme"]
    res = run_reset_case.meta_data["Res"]

    mask = (expected["reconScheme"] == scheme) & (expected["Res"] == res)

    expected_err_shape = expected[mask]["err_shape"].max()

    results["err_shape"].append(max_err_shape)
    results["err_mass"].append(max_err_mass)
    results["err_bound"].append(max_err_bound)
    results["runtime"].append(oftest.run_time(log)[0])
    results["reconScheme"].append(scheme)
    results["Res"].append(res)

    assert max_err_shape < expected_err_shape*1.05
    assert max_err_mass <= 1e-13
    assert max_err_bound <= 5e-5


def test_results():
    print("results", results)
    dir_name = os.path.dirname(os.path.abspath(__file__))
    res = pd.DataFrame(results)
    res['err_shape'] = res['err_shape'].apply(lambda x: '%.3E' % x)
    res['err_mass'] = res['err_mass'].apply(lambda x: '%.3E' % x)
    res['err_bound'] = res['err_bound'].apply(lambda x: '%.3E' % x)
    res.to_csv(os.path.join(dir_name, "results_shearedDisc.csv"),index=False)
