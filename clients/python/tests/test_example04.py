# Imports from the example:
from __future__ import print_function  # Python 2.7 compatibility
from srwpy.srwlib import *
from srwpy.uti_plot import *
import os

# Imports for tests:
import pytest


@pytest.mark.fast
def test_example04(example_code):
    exec(example_code, globals(), globals())

    # end of example, start testing
    assert len(arI1x) == 120
    assert len(arI1y) == 112
    assert hasattr(wfr, "mesh")
    assert hasattr(wfr, "partBeam")
    for param in ["zStart",
                  "eStart", "eFin", "ne",
                  "xStart", "xFin", "nx",
                  "yStart", "yFin", "ny"]:
        assert hasattr(wfr.mesh, param)
