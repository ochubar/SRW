#!/bin/bash
#
# Run an example (SRWLIB_Example10.py) to verify build
#
set -e

# Overridable parameters, e.g. timeout_s=2 bash basic-test.sh
: ${timeout_s:=20}
: ${ex_n:=10}
data_d=data_example_$ex_n
example_py=SRWLIB_Example$ex_n.py
out=$data_d/$example_py.log

# No color if stdout is not a tty
if [[ ! -t 1 ]]; then
    tput() { return; }
fi

# cleanup: kill hanging process, remove temporary directory, output files
pid=
to_remove=()
cleanup() {
    set +e
    if [[ -n $pid ]]; then
        kill -9 "$pid" && wait "$pid" >& /dev/null
        pid=
    fi
    local f
    for f in "${to_remove[@]}"; do
        rm -rf "$f"
    done
    to_remove=
}
trap cleanup EXIT ERR

# Start in examples directory
cd "$(dirname "$0")"
if [[ ! -d $data_d ]]; then
    mkdir -p "$data_d"
    to_remove+=( $data_d )
fi

# Run example with timeout
python -u "$example_py" >& "$out" &
pid=$!
status='PASS'
msg=' (timeout)'
color=2 # green
i=0
while (( $i < $timeout_s )); do
    sleep 1
    i=$((i + 1))
    if kill -0 "$pid" >& /dev/null; then
        continue
    fi
    code=$(wait "$pid" >&/dev/null || echo $?)
    pid=
    if (( $code )); then
        # Let the user know what went wrong
        tail "$out"
        status=FAIL
        msg=" (code=$code)"
        color=1 # red
        exit=1
    else
        status=PASS
        msg=
    fi
    break
done
echo "$example_py: $(tput setaf $color)$(tput bold)$status$(tput sgr0)$msg"
#TODO(robnagler) this should be generalized
to_remove+=( $data_d/ex${ex_n}_res_{int_se,prop_se,prop_me}.dat )
exit $code
