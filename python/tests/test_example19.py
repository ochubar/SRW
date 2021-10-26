# Imports from the example:
from __future__ import print_function  # Python 2.7 compatibility
from srwpy.srwlib import *
from srwpy.srwl_uti_smp import *
from srwpy import srwl_uti_smp_rnd_obj3d
from srwpy.uti_plot import *  # required for plotting
import os
import time
import numpy as np

# Imports for tests:
import pytest


@pytest.mark.fast
def test_example19(example_code):
    exec(example_code, globals(), globals())

    # end of example, start testing
    assert len(arI0) == 6084
    assert len(arI1) == 4194304
    assert len(arLogI1) == 4194304

    assert hasattr(wfr, "mesh")
    for param in ["zStart",
                  "eStart", "eFin", "ne",
                  "xStart", "xFin", "nx",
                  "yStart", "yFin", "ny"]:
        assert hasattr(wfr.mesh, param)
