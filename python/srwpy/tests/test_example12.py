import os
import subprocess
import pytest


@pytest.mark.fast
@pytest.mark.skipif(sys.platform == "linux", reason="runs only on Linux")
def test_example12(examples_dir):
    current_dir = os.getcwd()
    os.chdir(examples_dir)
    os.environ["DISPLAY"] = ""
    subprocess.run("python SRWLIB_Example12.py".split())
    os.chdir(current_dir)
