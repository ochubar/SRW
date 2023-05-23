# Imports from the example:
from __future__ import print_function  # Python 2.7 compatibility
from srwpy.srwlib import *
from srwpy.uti_plot import *  # required for plotting
import os
import sys
import time

# Imports for tests:
import pytest


@pytest.mark.fast
def test_example09(example_code):
    exec(example_code, globals(), globals())

    # end of example, start testing
    assert len(arI1x) == 504
    assert len(arI1y) == 416
    assert len(arI1) == len(arI1x) * len(arI1y)
    assert hasattr(wfr, "mesh")
    for param in ["zStart",
                  "eStart", "eFin", "ne",
                  "xStart", "xFin", "nx",
                  "yStart", "yFin", "ny"]:
        assert hasattr(wfr.mesh, param)
