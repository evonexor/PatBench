"""Documentation hygiene: relative links resolve and repository text files are English-only."""
import re

import pytest

from patbench import REPO

DOCS = sorted(
    p for p in REPO.rglob("*.md")
    if ".git" not in p.parts
    and ".pytest_cache" not in p.parts
)


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


CJK = re.compile(r"[\u3000-\u9fff\uff00-\uffef]")


@pytest.mark.parametrize(
    "doc",
    DOCS + [REPO / "pyproject.toml"],
    ids=lambda p: str(p.relative_to(REPO)),
)
def test_documentation_is_english_only(doc):
    assert not CJK.search(doc.read_text(encoding="utf-8")), f"{doc.name} contains non-English CJK characters"
