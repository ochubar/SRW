# Imports from the example:
from __future__ import print_function  # Python 2.7 compatibility
from srwpy.srwlib import *
from srwpy.uti_plot import *
import os
import sys

# Imports for tests:
import pytest


@pytest.mark.fast
def test_example06(example_code):
    exec(example_code, globals(), globals())

    # end of example, start testing
    assert len(powDenVsX) == 101
    assert len(powDenVsY) == 101
    for obj in [stkF, stkP]:
        assert hasattr(obj, "mesh")
        for param in ["zStart",
                      "eStart", "eFin", "ne",
                      "xStart", "xFin", "nx",
                      "yStart", "yFin", "ny"]:
            assert hasattr(obj.mesh, param)


@pytest.mark.fast
def test_example06_PETRA(example_code):
    exec(example_code, globals(), globals())

    # end of example, start testing
    assert len(arIS) == 20000
    assert len(arI) == 40401
    for obj in [stkF, stkP]:
        assert hasattr(obj, "mesh")
        for param in ["zStart",
                      "eStart", "eFin", "ne",
                      "xStart", "xFin", "nx",
                      "yStart", "yFin", "ny"]:
            assert hasattr(obj.mesh, param)
