# Imports from the example:
from __future__ import print_function  # Python 2.7 compatibility
import srwpy.uti_plot as uti_plot
from srwpy.srwlib import *
from srwpy.uti_math import fwhm
from scipy.special import jv

# Imports for tests:
import os
import pytest
from tests.conftest import filter_source, get_example_from_test


@pytest.mark.fast
def test_example16(examples_dir):
    current_dir = os.getcwd()
    os.chdir(examples_dir)
    os.environ["DISPLAY"] = ""

    code = filter_source(get_example_from_test())
    exec(code, globals(), globals())

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

    os.chdir(current_dir)
