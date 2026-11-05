"""Documentation hygiene: relative links resolve, and numbers stated in the docs match the manifest."""
import re

import pytest

from patbench import REPO, manifest

DOCS = sorted(p for p in REPO.rglob("*.md") if ".git" not in p.parts and ".pytest_cache" not in p.parts)


@pytest.mark.parametrize("doc", DOCS, ids=lambda p: str(p.relative_to(REPO)))
def test_relative_links_resolve(doc):
    text = doc.read_text(encoding="utf-8")
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    bad = []
    for target in re.findall(r"\]\(([^)\s]+)\)", text):
        if target.startswith(("http://", "https://", "#", "mailto:")):
            continue
        path = target.split("#")[0]
        if path and not (doc.parent / path).exists():
            bad.append(target)
    assert not bad, bad


def test_manifest_counts_quoted_in_the_docs():
    rows = manifest.read_manifest()
    n = {s: sum(1 for r in rows if r["status"] == s) for s in ("verbatim", "patched-paths", "authored")}
    prov = (REPO / "docs/CODE_PROVENANCE.md").read_text(encoding="utf-8")
    assert f"{n['verbatim']} verbatim" in prov and f"{n['patched-paths']} patched" in prov and f"{n['authored']} authored" in prov
    assert f"{len(rows)} files" in prov


CJK = re.compile(r"[\u3000-\u9fff\uff00-\uffef]")


@pytest.mark.parametrize("doc", DOCS + [REPO / "pyproject.toml", REPO / "Makefile"], ids=lambda p: str(p.relative_to(REPO)))
def test_documentation_is_english_only(doc):
    assert not CJK.search(doc.read_text(encoding="utf-8")), f"{doc.name} contains Chinese characters"
