# Imports from the example:
from __future__ import print_function  # Python 2.7 compatibility
from srwpy.srwlib import *
from srwpy.uti_plot import *
import os

# Imports for tests:
import pytest


@pytest.mark.fast
def test_example01(example_code):
    exec(example_code, globals(), globals())

    # end of example, start testing
    assert len(partTraj.arX) == npTraj
    assert len(partTraj.arY) == npTraj
    assert ctMesh == [0, 2.46, npTraj]


@pytest.mark.fast
def test_example01_kick_matr(example_code):
    exec(example_code, globals(), globals())

    # end of example, start testing
    assert len(partTraj.arX) == 10001
    assert len(partTraj.arY) == 10001
    assert ctMesh == [0, 2.5, 10001]
