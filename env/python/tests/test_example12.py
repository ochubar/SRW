import glob
import pytest
import subprocess
import sys


@pytest.mark.slow
@pytest.mark.interactive
@pytest.mark.skipif(sys.platform != "linux", reason="runs only on Linux")
# @pytest.mark.skipif(sys.platform not in ["linux", "darwin"], reason="runs only on Unix")
def test_example12(temp_example_file):
    subprocess.run(f"time mpirun --oversubscribe -n 4 python {temp_example_file}".split(),
                   check=True)

    # The log files are not produced by this example.
    log_files = glob.glob("__srwl_logs__/srwl_stat_wfr_emit_prop_multi_e_*.log",
                          recursive=True)

    with pytest.raises(IndexError):
        last_log = sorted(log_files)[-1]
