import os
import pytest

REASON = "The test is executed via a separate step in CI"

@pytest.mark.skip(reason=REASON)
@pytest.mark.slow
def test_example08(examples_dir):
    current_dir = os.getcwd()
    os.chdir(examples_dir)
    os.environ["DISPLAY"] = ""

    print(REASON)

    os.chdir(current_dir)
