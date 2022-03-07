import glob
import pytest
import subprocess
import sys


@pytest.mark.slow
@pytest.mark.parametrize("temp_example_file", [20], indirect=True)
@pytest.mark.skipif(sys.platform not in ["linux", "darwin"], reason="runs only on Unix")
def test_example10(temp_example_file):
    subprocess.run(f"time mpirun --oversubscribe -n 4 python {temp_example_file}".split(),
                   check=True)

    log_files = glob.glob("__srwl_logs__/srwl_stat_wfr_emit_prop_multi_e_*.log",
                          recursive=True)

    last_log = sorted(log_files)[-1]
    with open(last_log, "r") as f:
        content = f.read()
    print(f"Log file:\n{content}")

    assert content
