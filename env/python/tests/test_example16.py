# Imports from the example:
from __future__ import print_function  # Python 2.7 compatibility
from srwpy.srwlib import *
from srwpy.uti_math import fwhm
from srwpy.uti_plot import *
from scipy.special import jv

# Imports for tests:
import pytest


@pytest.mark.fast
def test_example16(example_code):
    exec(example_code, globals(), globals())

    # end of example, start testing
    assert len(sIn) == 1120
    assert len(arIinY) == 1120
    assert len(sOut) == 1152
    assert len(arI1y) == 1152
    assert len(analyticalIntens) == 1152

    assert hasattr(wfr, "mesh")
    for param in ["zStart",
                  "eStart", "eFin", "ne",
                  "xStart", "xFin", "nx",
                  "yStart", "yFin", "ny"]:
        assert hasattr(wfr.mesh, param)
