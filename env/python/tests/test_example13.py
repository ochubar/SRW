# Imports from the example:
from __future__ import print_function  # Python 2.7 compatibility
from srwpy.srwlib import *
from srwpy.uti_plot import *  # required for plotting
import time

# Imports for tests:
import pytest


@pytest.mark.fast
def test_example13(example_code):
    exec(example_code, globals(), globals())

    # end of example, start testing
    assert len(arI0) == 102900
    assert len(arI1s) == 6720
    assert len(arI1m) == 6720
    assert hasattr(wfr, "mesh")
    for param in ["zStart",
                  "eStart", "eFin", "ne",
                  "xStart", "xFin", "nx",
                  "yStart", "yFin", "ny"]:
        assert hasattr(wfr.mesh, param)
