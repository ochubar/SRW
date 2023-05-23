# Imports from the example:
from __future__ import print_function  # Python 2.7 compatibility
from srwpy.srwlib import *
from srwpy.uti_plot import *
import os

# Imports for tests:
import pytest


@pytest.mark.fast
def test_example03(example_code):
    exec(example_code, globals(), globals())

    # end of example, start testing
    assert len(arI2x) == 101
    assert len(arI2y) == 101
    for obj in [wfr1, wfr2]:
        assert hasattr(obj, "mesh")
        assert hasattr(obj, "partBeam")
        for param in ["zStart",
                      "eStart", "eFin", "ne",
                      "xStart", "xFin", "nx",
                      "yStart", "yFin", "ny"]:
            assert hasattr(obj.mesh, param)
