# Imports from the example:
from __future__ import print_function
import srwpy.uti_plot as uti_plot
from srwpy.srwlib import *
from srwpy.uti_math import matr_prod, fwhm

# Imports for tests:
import os
import pytest
from tests.conftest import filter_source, get_example_from_test


@pytest.mark.fast
def test_example15(examples_dir):
    current_dir = os.getcwd()
    os.chdir(examples_dir)
    os.environ["DISPLAY"] = ""

    code = filter_source(get_example_from_test())
    exec(code, globals(), globals())

    # end of example, start testing
    assert type(intensitiesToPlot) == dict
    assert len(intensitiesToPlot["intensity"]) == 3
    assert len(intensitiesToPlot["intensity"][0]) == 235872
    assert len(intensitiesToPlot["intensity"][1]) == 4212
    assert len(intensitiesToPlot["intensity"][2]) == 3024

    assert len(intensitiesToPlot["mesh_x"]) == 3
    assert len(intensitiesToPlot["mesh_x"][0]) == 3
    assert len(intensitiesToPlot["mesh_x"][1]) == 3
    assert len(intensitiesToPlot["mesh_x"][2]) == 3

    assert len(intensitiesToPlot["mesh_y"]) == 3
    assert len(intensitiesToPlot["mesh_y"][0]) == 3
    assert len(intensitiesToPlot["mesh_y"][1]) == 3
    assert len(intensitiesToPlot["mesh_y"][2]) == 3

    assert hasattr(wfr, "mesh")
    for param in ["zStart",
                  "eStart", "eFin", "ne",
                  "xStart", "xFin", "nx",
                  "yStart", "yFin", "ny"]:
        assert hasattr(wfr.mesh, param)

    os.chdir(current_dir)
