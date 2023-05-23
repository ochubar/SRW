import os
import pytest
import shutil


@pytest.fixture(scope="function")
def examples_dir():
    _examples_dir = os.path.abspath(os.path.join(
        os.path.dirname(
            os.path.dirname(__file__) #../
        ), "srwpy", "examples")
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
            startwhites = line[:len(line)-len(lstripped_line)]
            if lstripped_line.startswith("uti_plot_show"):
                lines[i] = f"{startwhites}pass # {line}"
            elif lstripped_line.startswith(("from ", "import ")):
                lines[i] = f"{startwhites}pass # {line}"
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
    current_test = os.environ.get("PYTEST_CURRENT_TEST")
    # This var will look like:
    #   tests/test_example10.py::test_example10[20] (setup)
    test_name = current_test.split(":")[-1].split(" ")[0].split("[")[0]
    _example_file = f"{test_name.replace('test_e', 'SRWLIB_E')}.py"
    return _example_file


@pytest.fixture(scope="function")
def no_display():
    os.environ["DISPLAY"] = ""


@pytest.fixture(scope="function")
def example_code(examples_dir, example_file, no_display):
    code = filter_source(example_file, filter=True, as_string=True)
    return code


@pytest.fixture(scope="function")
def raw_example_lines(examples_dir, example_file, no_display):
    lines = filter_source(example_file, filter=False, as_string=False)
    return lines


@pytest.fixture(scope="function")
def temp_example_file(tmp_path, example_file, raw_example_lines, request):
    output_file = tmp_path / example_file
    print(f"Output file: {output_file}")

    lines = raw_example_lines

    if hasattr(request, "param"):
        param_name = "nMacroElec"
        param_old_value = 50000
        param_new_value = request.param

        search_str = f"{param_name} = {param_old_value}"
        replace_str = f"{param_name} = {param_new_value}"

        for i, line in enumerate(lines):
            if search_str in line:
                print(f"\nReplacing {param_name} from {param_old_value} --> {param_new_value}\n")
                lines[i] = line.replace(search_str, replace_str)

    code = "\n".join(lines)
    output_file.write_text(code)

    yield output_file

    shutil.rmtree("__srwl_logs__", ignore_errors=True)
