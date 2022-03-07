import glob
import pytest
import subprocess
import sys


@pytest.mark.slow
@pytest.mark.skipif(sys.platform not in ["darwin"], reason="runs only on macOS")
def test_example20(temp_example_file):
    subprocess.run(f"time mpirun --oversubscribe -n 4 python {temp_example_file} "
                   f"--wm --wm_nop --wm_ch=6 --wm_nm=10 --wm_nmm=2 --wm_na=5 "
                   f"--wm_ns=5 --wm_ff=h5 --wm_fni=ex20_res".split(),
                   check=True)

    # This step takes too long, commenting it out for now:
    # time python SRWLIB_Example20.py --wm_ch=7 --wm_fnmi=ex20_res_mi.h5 --wm_ncm=10 --wm_fni=ex20_res

    log_files = glob.glob("__srwl_logs__/srwl_stat_wfr_emit_prop_multi_e_*.log",
                          recursive=True)

    last_log = sorted(log_files)[-1]
    with open(last_log, "r") as f:
        content = f.read()
    print(f"Log file:\n{content}")

    assert content
