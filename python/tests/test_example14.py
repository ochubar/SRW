# Imports from the example:
from __future__ import print_function  # Python 2.7 compatibility
from srwpy.srwlib import *
from srwpy.uti_math import *
from srwpy.uti_plot import *  # required for plotting
import os
import time

# Imports for tests:
import pytest


@pytest.mark.fast
def test_example14(example_code):
    exec(example_code, globals(), globals())

    # end of example, start testing
    assert arF1vsXY is None
    assert arPI1vsT is None
    assert arP1vsT is None
    assert hasattr(wfr, "mesh")
    for param in ["zStart",
                  "eStart", "eFin", "ne",
                  "xStart", "xFin", "nx",
                  "yStart", "yFin", "ny"]:
        assert hasattr(wfr.mesh, param)
