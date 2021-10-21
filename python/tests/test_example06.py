# Imports from the example:
from __future__ import print_function #Python 2.7 compatibility
from srwpy.srwlib import *
from srwpy.uti_plot import *
import os
import sys

# Imports for tests:
import os
import pytest
from tests.conftest import filter_source, get_example_from_test


@pytest.mark.fast
def test_example06(examples_dir):
    current_dir = os.getcwd()
    os.chdir(examples_dir)

    code = filter_source(get_example_from_test())
    exec(code, globals(), globals())

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

    os.chdir(current_dir)


@pytest.mark.fast
def test_example06_PETRA(examples_dir):
    current_dir = os.getcwd()
    os.chdir(examples_dir)

    code = filter_source(get_example_from_test())
    exec(code, globals(), globals())

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

    os.chdir(current_dir)
