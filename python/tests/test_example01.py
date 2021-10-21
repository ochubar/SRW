# Imports from the example:
from __future__ import print_function #Python 2.7 compatibility
from srwpy.srwlib import *
from srwpy.uti_plot import *
import os

# Imports for tests:
import os
import pytest
from tests.conftest import filter_source, get_example_from_test


@pytest.mark.fast
def test_example01(examples_dir):
    current_dir = os.getcwd()
    os.chdir(examples_dir)

    code = filter_source(get_example_from_test())
    exec(code, globals(), globals())

    # end of example, start testing
    assert len(partTraj.arX) == npTraj
    assert len(partTraj.arY) == npTraj
    assert ctMesh == [0, 2.46, npTraj]

    os.chdir(current_dir)


@pytest.mark.fast
def test_example01_kick_matr(examples_dir):
    current_dir = os.getcwd()
    os.chdir(examples_dir)

    code = filter_source(get_example_from_test())
    exec(code, globals(), globals())

    # end of example, start testing
    assert len(partTraj.arX) == 10001
    assert len(partTraj.arY) == 10001
    assert ctMesh == [0, 2.5, 10001]

    os.chdir(current_dir)
