# Imports from the example:
from __future__ import print_function  # Python 2.7 compatibility
from srwpy.srwlib import *
from srwpy.uti_plot import *
import os

# Imports for tests:
import pytest


@pytest.mark.fast
def test_example02(example_code):
    exec(example_code, globals(), globals())

    # end of example, start testing
    assert len(partTraj.arX) == npTraj
    assert len(partTraj.arY) == npTraj
    assert ctMesh == [0, 12.2, npTraj]
