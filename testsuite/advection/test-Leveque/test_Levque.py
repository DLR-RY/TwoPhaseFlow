import os
import pytest
import oftest
import pandas as pd
from oftest import run_case
from oftest import run_reset_case


@pytest.fixture  # descruct all tests
def load_errorfiles():
    dir_name = os.path.dirname(os.path.abspath(__file__))
    # dir_name = oftest.base_dir()
    ferr = os.path.join(dir_name, "volumeFractionError.dat")
    err = oftest.read_functionObject(ferr)
    return err


def simMod(res, reconScheme):
    dir_name = os.path.dirname(os.path.abspath(__file__))
    filemod = {
        "system/simulationParameter": [
            ("RECONSCHEME", reconScheme),
            ("nx", res),
            ("ny", res),
            ("nz", res),
        ]
    }
    case_mod = oftest.Case_modifiers(filemod, dir_name)
    meta_data = {"reconScheme": reconScheme, "Res": res}
    case_mod = oftest.Case_modifiers(filemod, dir_name, meta_data)
    return case_mod

c_isoAlpha32 = simMod(32, "isoAlpha")
c_isoAlpha64 = simMod(64, "isoAlpha")
c_isoAlpha128 = simMod(128, "isoAlpha")

c_plicRDF32 = simMod(32, "plicRDF")
c_plicRDF64 = simMod(64, "plicRDF")
c_plicRDF128 = simMod(128, "plicRDF")

parameters = [
    (c_isoAlpha32,  8.661e-03),
    (c_isoAlpha64,  3.391E-03),
    (c_isoAlpha128, 5.448E-04),
    (c_plicRDF32,   7.743E-03),
    (c_plicRDF64,   3.079E-03),
    (c_plicRDF128,  6.969E-04),
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
    "run_reset_case, expected",
    parameters,
    indirect=["run_reset_case"],
    ids=[
        "isoAlpha32",
        "isoAlpha64",
        "isoAlpha128",
        "plicRDF32",
        "plicRDF64",
        "plicRDF128",
    ]
)
def test_deformation(run_reset_case, expected, load_errorfiles):

    print("runcase", run_reset_case)
    log = oftest.path_log()
    err = load_errorfiles
    max_err_shape = abs(err.iloc[1, 1].max())
    max_err_mass = abs(err.iloc[:, 2].max())
    max_err_bound = abs(err.iloc[:, 3].max())

    results["err_shape"].append(max_err_shape)
    results["err_mass"].append(max_err_mass)
    results["err_bound"].append(max_err_bound)
    results["runtime"].append(oftest.run_time(log)[0])
    results["reconScheme"].append(run_reset_case.meta_data["reconScheme"])
    results["Res"].append(run_reset_case.meta_data["Res"])


    assert max_err_shape < expected*1.05
    assert max_err_mass <= 1e-13
    assert max_err_bound <= 5e-6
    assert oftest.case_status(log) == "completed"  # checks if run completes


def test_results():
    print("results", results)
    dir_name = os.path.dirname(os.path.abspath(__file__))
    res = pd.DataFrame(results)
    res['err_shape'] = res['err_shape'].apply(lambda x: '%.3E' % x)
    res['err_mass'] = res['err_mass'].apply(lambda x: '%.3E' % x)
    res['err_bound'] = res['err_bound'].apply(lambda x: '%.3E' % x)
    res.to_csv(os.path.join(dir_name, "results_leveque.csv"),index=False)
