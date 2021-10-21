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


def filter_source(file):
    with open(file, "r") as f:
        lines = f.readlines()
        lines = [line.rstrip() for line in lines]

    for i, line in enumerate(lines):
        lstripped_line = line.lstrip()
        if lstripped_line.startswith("uti_plot_show"):
            lines[i] = f"# {line}"
        elif lstripped_line.startswith(("from ", "import ")):
            lines[i] = f"# {line}"
        else:
            pass

    code = "\n".join(lines)
    # print(code)

    return code


def get_example_from_test():
    test_name = os.environ.get("PYTEST_CURRENT_TEST").split(":")[-1].split(" ")[0]
    example_file = f"{test_name.replace('test_e', 'SRWLIB_E')}.py"
    return example_file
