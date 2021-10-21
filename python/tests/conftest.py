import os
import pytest


@pytest.fixture
def examples_dir():
    _examples_dir = os.path.abspath(os.path.join(
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(__file__)
            )
        ), "examples", "python")
    )
    return _examples_dir
