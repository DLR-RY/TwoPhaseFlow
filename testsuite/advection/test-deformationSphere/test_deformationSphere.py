import os
import pytest
import oftest
import pandas as pd
from oftest import run_case
from oftest import run_reset_case


@pytest.fixture  # descruct all tests
def load_errorfiles():
    dir_name = os.path.dirname(os.path.abspath(__file__))
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
            ("nz", 2 * res),
        ]
    }
    meta_data = {"reconScheme": reconScheme, "Res": res}
    case_mod = oftest.Case_modifiers(filemod, dir_name, meta_data)
    return case_mod


c_isoAlpha32 = simMod(32, "isoAlpha")  # 0.004244
c_isoAlpha64 = simMod(64, "isoAlpha")  # 0.001237
c_isoAlpha128 = simMod(128, "isoAlpha")  # 0.001237

c_plicRDF32 = simMod(32, "plicRDF")  # 0.00558
c_plicRDF64 = simMod(64, "plicRDF")  # 0.001494
c_plicRDF128 = simMod(128, "plicRDF")  # 0.001494


parameters = [
    (c_isoAlpha32,  4.466E-03),
    (c_isoAlpha64,  1.339E-03),
    (c_isoAlpha128, 3.543E-04),
    (c_plicRDF32,   4.471E-03),
    (c_plicRDF64,   1.270E-03),
    (c_plicRDF128,  3.176E-04),
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

    assert max_err_shape < (expected*1.05)
    assert max_err_mass <= 1e-13
    assert max_err_bound <= 5e-5
    assert oftest.case_status(log) == "completed"  # checks if run completes


def test_write_results():
    print("results", results)
    dir_name = os.path.dirname(os.path.abspath(__file__))
    res = pd.DataFrame(results)
    res['err_shape'] = res['err_shape'].apply(lambda x: '%.3E' % x)
    res['err_mass'] = res['err_mass'].apply(lambda x: '%.3E' % x)
    res['err_bound'] = res['err_bound'].apply(lambda x: '%.3E' % x)
    res.to_csv(os.path.join(dir_name, "results_deformation.csv"),index=False)
