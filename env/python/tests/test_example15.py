# Imports from the example:
from __future__ import print_function
from srwpy.srwlib import *
from srwpy.uti_math import matr_prod, fwhm
from srwpy.uti_plot import *

# Imports for tests:
import pytest


@pytest.mark.fast
def test_example15(example_code):
    exec(example_code, globals(), globals())

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
