import os
from numpy import logspace
import pytest
import oftest
from oftest import run_reset_case

def test_tiltedBox(run_reset_case):
    logs = oftest.path_logs()
    assert len(logs) > 0
    for log in logs:
        assert oftest.case_status(log) == 'completed' # checks if run completes

