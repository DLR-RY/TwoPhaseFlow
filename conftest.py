import sys
import pytest

# PyFoam 2023.7 bundles its own ThirdParty/six.py which lacks a `moves`
# submodule.  Python 3.12 is strict about submodule imports, so
# `from PyFoam.ThirdParty.six.moves import X` raises ModuleNotFoundError.
# Pre-register the real six.moves under that name so the import succeeds.
try:
    import six
    sys.modules.setdefault('PyFoam.ThirdParty.six.moves', six.moves)
except ImportError:
    pass


def pytest_addoption(parser):
    parser.addoption(
        "--writeNSteps", action="store", default=0, help="only perform specified number of timestep"
    )
    parser.addoption(
        "--no-clean-up", action='store_false',default=True ,help="do not clean case after run"
    )
