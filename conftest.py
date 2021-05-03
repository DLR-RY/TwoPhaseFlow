import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--writeNSteps", action="store", default=0, help="only perform specified number of timestep"
    )
    parser.addoption(
        "--no-clean-up", action='store_false',default=True ,help="do not clean case after run"
    )
