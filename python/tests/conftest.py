import os
import pytest


@pytest.fixture(scope="function")
def examples_dir():
    _examples_dir = os.path.abspath(os.path.join(
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(__file__)
            )
        ), "examples", "python")
    )

    current_dir = os.getcwd()
    os.chdir(_examples_dir)

    yield _examples_dir

    os.chdir(current_dir)


def filter_source(file, filter=True, as_string=True):
    with open(file, "r") as f:
        lines = f.readlines()
        lines = [line.rstrip() for line in lines]

    if filter:
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
    if as_string:
        return code
    else:
        return lines


@pytest.fixture(scope="function")
def example_file():
    test_name = os.environ.get("PYTEST_CURRENT_TEST").split(":")[-1].split(" ")[0]
    _example_file = f"{test_name.replace('test_e', 'SRWLIB_E')}.py"
    return _example_file


@pytest.fixture(scope="function")
def no_display():
    os.environ["DISPLAY"] = ""


@pytest.fixture(scope="function")
def example_code(examples_dir, example_file, no_display):
    code = filter_source(example_file, filter=True, as_string=True)
    return code


def get_example_from_test():
    test_name = os.environ.get("PYTEST_CURRENT_TEST").split(":")[-1].split(" ")[0]
    example_file = f"{test_name.replace('test_e', 'SRWLIB_E')}.py"
    return example_file
