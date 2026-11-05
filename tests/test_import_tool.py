"""tools/import_pipeline.py rewrites path literals and nothing else."""
import ast
import sys
import textwrap
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
import import_pipeline as ip  # noqa: E402

P = ip.SERVER_PREFIX


def patch(src, rel="pipeline/08_build/x.py"):
    return ip.patch_source(textwrap.dedent(src).lstrip("\n"), rel)


def test_single_and_double_quoted_literals():
    out, n = patch(f'''
        R = "{P}/outputs/a/results/"
        S = '{P}/raw/patent-data/'
        ''')
    assert n == 2
    assert 'R = _ws + "/outputs/a/results/"' in out and "S = _ws + '/raw/patent-data/'" in out
    ast.parse(out)


def test_ws_defined_after_docstring_and_future_import():
    out, _ = patch(f'''
        """Module docstring."""
        from __future__ import annotations
        import os
        R = "{P}/x"
        ''')
    lines = out.splitlines()
    assert lines[0] == '"""Module docstring."""'
    assert lines[1].startswith("from __future__")
    assert lines[2].startswith("_ws = ")
    assert out.count("_ws = ") == 1


def test_docstring_and_comments_are_left_alone():
    src = f'''
        """Reads {P}/outputs/x/ and "{P}/y"."""
        # also {P}/z
        a = 1
        '''
    out, n = patch(src)
    assert n == 0 and out == textwrap.dedent(src).lstrip("\n")


def test_text_without_the_prefix_is_unchanged():
    src = "x = '/data/other/path'\n"
    assert ip.patch_source(src, "pipeline/08_build/x.py") == (src, 0)


def test_expression_forms_stay_valid():
    out, n = patch(f'''
        from pathlib import Path
        a = Path("{P}/raw/")
        b = os.path.join("{P}/o", "q")
        c = "{P}/a" if flag else "{P}/b"
        d = [("{P}/m", 1)]
        ''')
    assert n == 5
    ast.parse(out)
    assert 'Path(_ws + "/raw/")' in out


def test_v1_scripts_directory_becomes_the_directory_of_the_file():
    out, _ = patch(f"import sys\nsys.path.insert(0,'{ip.V1_SCRIPTS}')\n", "pipeline/05_target_map/v1/a.py")
    assert "__file__" in out and "'..', 'v1'" not in out


def test_v2_files_reach_v1_through_a_sibling_directory():
    out, _ = patch(f"import sys\nsys.path.insert(0,'{ip.V1_SCRIPTS}')\n", "pipeline/05_target_map/v2/b.py")
    assert "'..', 'v1'" in out


def test_v2_assemble_pattern():
    out, _ = patch(f"V1='{P}/outputs/target-uniprot-map-20261016/'\nsys.path.insert(0,V1+'scripts')\n",
                   "pipeline/05_target_map/v2/v2_assemble.py")
    assert "V1=_ws + '/outputs/target-uniprot-map-20261016/'" in out
    assert "V1+'scripts'" not in out and "'..', 'v1'" in out


@pytest.mark.skipif(sys.version_info < (3, 12), reason="f-string tokens differ before Python 3.12")
def test_fstring_prefix():
    out, n = patch(f'x = f"{P}/outputs/{{name}}/y"\n')
    assert n == 1 and 'f"{_ws}/outputs/{name}/y"' in out


def test_the_edited_file_differs_only_by_the_literal_and_the_ws_line():
    src = f'"""d."""\nimport sys\nD = "{P}/outputs/q/"\nprint(D)\n'
    out, _ = ip.patch_source(src, "pipeline/08_build/x.py")
    old = src.splitlines()
    new = [l for l in out.splitlines() if not l.startswith("_ws = ")]
    assert [l.replace('_ws + "', '"').replace('"/outputs/q/"', f'"{P}/outputs/q/"') for l in new] == old


def test_public_database_prefix_is_rewritten_to_the_raw_data_directory():
    out, n = patch(f'''
        from pathlib import Path
        ROOT = Path("{ip.PUBLIC_PREFIX}")
        STD = Path("{ip.PUBLIC_PREFIX}/standardized")
        ''')
    assert n == 2
    assert 'ROOT = Path(_ws + "/raw/patent-data")' in out
    assert 'STD = Path(_ws + "/raw/patent-data/standardized")' in out
    ast.parse(out)
