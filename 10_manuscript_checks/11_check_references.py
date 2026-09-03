"""
11_check_references.py -- the two-way cross-check of the reference apparatus

Purpose  Two directions, both must balance:
           (a) every reference in the list is cited in the text
           (b) every citation in the text has an entry in the list

Citation style in the draft is author-year in square brackets, several
sources separated by semicolons:  [Unger et al. 2023]  -  [Cameron et al.
2015; Rajpar et al. 2009]. The conversion into the numbered Cell Press form
happens at package build time (21_build_submission.py), not here.

Input    manuscript/MANUSCRIPT.md
Output   results/reference_check.txt (and return code 1 if anything is open)
Runtime  seconds
"""
from __future__ import annotations

import os
import pathlib
import re
import sys

_env = os.environ.get("PAPER_V2_ROOT")
ROOT = (pathlib.Path(_env) if _env
        else pathlib.Path(__file__).resolve().parents[1])
MS = ROOT / "manuscript" / "MANUSCRIPT.md"
RES = ROOT / "results"
RES.mkdir(parents=True, exist_ok=True)

OUTPUT: list[str] = []


try:  # the Windows console is cp1252; beta and rho would break it otherwise
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass


def say(s: str = "") -> None:
    print(s)
    OUTPUT.append(s)


# Bracket expressions that are NOT citations: references to our own figures,
# tables, files, status and significance marks.
# The trailing $ anchors EVERY alternative to the full string: a bare "[A-Z]"
# must never swallow "Unger et al. 2023" just because it starts with a
# capital letter.
NOT_A_CITATION = re.compile(
    r"^(?:[A-Z]|C|E|F|\d+"
    r"|(?:Fig\.?|Figure|Table|TS|S)\s?\d[\w., ]*"
    r"|(?:z|P|OR|n|rho|\u03c1)\s.*"
    r"|sic)$")


def key_of(z: str) -> str:
    """Reduce a citation to its comparison key: surname + year."""
    z = z.strip()
    z = re.sub(r"\s+", " ", z)
    return z.lower()


def citations_from_text(text: str) -> list[str]:
    """All author-year citations of the running text, without the list."""
    head = text.split("\n## References")[0]
    # Code blocks and inline code out -- file names live there.
    head = re.sub(r"```.*?```", " ", head, flags=re.S)
    head = re.sub(r"`[^`]*`", " ", head)
    hits: list[str] = []
    for raw in re.findall(r"\[([^\[\]]{4,1400})\]", head):
        if NOT_A_CITATION.match(raw.strip()):
            continue
        if not re.search(r"\b(19|20)\d{2}\b", raw):
            continue
        for single in raw.split(";"):
            single = single.strip()
            if single:
                hits.append(single)
    return hits


def references_from_list(text: str) -> list[tuple[int, str]]:
    """The numbered entries of the reference list."""
    if "\n## References" not in text:
        return []
    part = text.split("\n## References", 1)[1]
    entries: list[tuple[int, str]] = []
    for m in re.finditer(r"^(\d+)\.\s+(.+?)(?=^\d+\.\s|\Z)", part,
                         flags=re.M | re.S):
        entries.append((int(m.group(1)), re.sub(r"\s+", " ",
                                                m.group(2)).strip()))
    return entries


def fits(citation: str, reference: str) -> bool:
    """Citation and reference belong together when the year matches and one
    carrying word of the citation appears in the reference.

    Three characters as the lower bound, not four: "Hao" and "Liu" are
    surnames, and with four characters exactly these two references failed.
    """
    yc = set(re.findall(r"((?:19|20)\d{2})", citation))
    yr = set(re.findall(r"((?:19|20)\d{2})", reference))
    if yc and yr and not (yc & yr):
        return False
    filler = {"et", "al", "and", "the", "review", "for", "with", "von"}
    words = [w for w in re.findall(r"[A-Za-zÄÖÜäöüß.\-]{3,}", citation)
             if w.lower().strip(".-") not in filler]
    ref = reference.lower()
    if any(w.lower().strip(".-") in ref for w in words):
        return True
    # Two-letter surnames ("He", "Wu") fall through the word length. The
    # first token of a citation is always the surname; it is searched at a
    # word boundary of the reference.
    first = citation.strip().split()[0].strip(".,-")
    return bool(first) and re.search(r"\b" + re.escape(first.lower()) + r"\b",
                                     ref) is not None


def main() -> int:
    ms = MS.read_text(encoding="utf-8")
    citations = citations_from_text(ms)
    refs = references_from_list(ms)

    say("11_check_references.py -- two-way cross-check")
    say(f"  references in the list    : {len(refs)}")
    say(f"  citations in running text : {len(citations)} "
        f"({len(set(map(key_of, citations)))} distinct)")
    say("")

    # (a) Reference without a citation. Two series carry no publication and
    #     are cited via their accession -- the Cell Press policy explicitly
    #     allows this. They count as cited when their accession appears in
    #     the running text.
    head = ms.split("\n## References")[0]

    def cited(r: str) -> bool:
        if any(fits(z, r) for z in citations):
            return True
        akz = re.findall(r"\b(GSE\d{4,}|E-MTAB-\d+)\b", r)
        return bool(akz) and all(a in head for a in set(akz))

    without_citation = [(n, r) for n, r in refs if not cited(r)]
    say("(a) References NOT cited in the text")
    if not without_citation:
        say("    none -- every reference is cited.")
    for n, r in without_citation:
        say(f"    [{n:2d}] {r[:96]}")
    say("")

    # (b) Citation without a reference
    without_ref = sorted({z for z in citations
                          if not any(fits(z, r) for _, r in refs)},
                         key=key_of)
    say("(b) Citations in the text WITHOUT an entry in the list")
    if not without_ref:
        say("    none -- every citation has an entry.")
    for z in without_ref:
        say(f"    {z}")
    say("")

    open_items = len(without_citation) + len(without_ref)
    say("")
    say(f"=== {open_items} open items in the cross-check ===")

    p = RES / "reference_check.txt"
    p.write_text("\n".join(OUTPUT) + "\n", encoding="utf-8")
    print("->", p)
    return 1 if open_items else 0


if __name__ == "__main__":
    sys.exit(main())
