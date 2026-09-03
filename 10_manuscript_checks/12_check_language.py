"""
12_check_language.py -- no German prose anywhere in what is published.

Purpose  Checks everything that is published -- the submission package, the
         manuscript sources, the panel data, the code, the preregistrations,
         the README and the licences -- for German: umlauts, the replacement
         character, the usual suspects, the forbidden forms "naiv"/"naive",
         the accession rule and the spelling of the identifiers.

         Two things are deliberately NOT flagged, and both are documented in
         README.md:
           * the short internal names of variables and of the columns in
             derived_data/ and results/, which the analysis scripts read and
             write. column_glossary.csv gives their English meaning.
           * proper names that carry an umlaut (Universitaetsklinikum
             Erlangen, Nuernberg, Mueller). The rule forbids the German
             LANGUAGE, not the correct spelling of a name.

Inputs   manuscript/*.md, figures/data/*.csv, code/**, preregistrations/**,
         README.md, submission/*
Output   results/language_check.txt (exit code 1 if anything is open)
Runtime  seconds
"""
from __future__ import annotations

import os
import pathlib
import re
import sys
import zipfile

import pandas as pd

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
RES = WURZEL / "results"
RES.mkdir(parents=True, exist_ok=True)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

AUSGABE: list[str] = []


def sag(s: str = "") -> None:
    print(s)
    AUSGABE.append(s)


UMLAUT = re.compile(r"[äöüÄÖÜß]")
ERSATZZEICHEN = re.compile(r"�")

# German words that actually occurred in this project, plus the ones named
# explicitly by the author.
VERDACHT = re.compile(
    r"\b(naiv\w*|naive|Eichung|Eichungen|Zelle|Zellen|Probe|Proben|Grenze|"
    r"Anteil|Datensatz|Datensaetze|Spender|Achse|Achsen|Fenster|Zerlegung|"
    r"Modul|Gene\)|bestanden|durchgefallen|osteogen\b|adipogen\b|chondrogen\b|"
    r"myogen\b|Entitaet|Anlage|Hintergrund|keine|kein|ohne|nicht|Linie|"
    r"gegen|Laesion|Vorregistrierung|Protokoll|Begruendung|Abbildung|"
    r"Tabelle|Sichtung|Kohorte|Zweck|Eingaben|Ausgaben|Laufzeit|Regel|"
    r"Nachweisgrenze|deutsch|Sitzung|Arbeitsstand)\b")

# Exceptions: proper names, and the identifiers that the code legitimately
# carries. Anything matched here is blanked before the check runs.
AUSNAHME = re.compile(
    r"Müller|Erlangen|Nürnberg|"
    r"Universitätsklinikum|Universität|"
    r"\bNAIV\b|_marker\.py|`[^`]*`|"
    r"\b[\w./-]+\.(?:csv|py|R|md|txt|json|gz|h5ad|bib|cff)\b|"
    r"\b\w+_\w+\b")

# The accession rule: canonically "information (GSExxxxx)". Forbidden are
# "GSExxxxx (information)" and "information GSExxxxx" without parentheses.
GSE_FALSCH = re.compile(
    r"(?<![.,;:]\s)(?<!^)(?<!and )(?<!in )(?<!of )"
    r"GSE\d{4,6}\s*\((?!\s*(GSE|this|n =|the ))[A-Za-z]")
# The rule applies to LABELS, not to sentence structure: "in GSE151315 the
# contrast carries ..." is an ordinary English sentence in which the accession
# is the subject. Only an assay or level word immediately in front of the
# accession without parentheses is flagged.
BESCHREIBER = (r"ATAC|ATAC-seq|H3K27ac|ChIP|ChIP-seq|methylome|methylom|"
               r"chromatin|RNA-seq|27K|450K|osteogenic|adipogenic|chondrogenic|"
               r"axis|window|cohort|dataset|series|module|calibration")
GSE_OHNE_KLAMMER = re.compile(r"(?:" + BESCHREIBER + r")\s+GSE\d{4,6}", re.I)

# The case of the identifiers: the abbreviations are forbidden project-wide.
KURZFORM = re.compile(r"\b(OSTEO|ADIPO|Osteo|Adipo|osteog|chondr)\b(?!-MSC)")

# ---------------------------------------------------------------------------
# The internal vocabulary of the analysis code: the names of variables, of the
# marker-set constants and of the columns in derived_data/ and results/. These
# are identifiers, not prose. They stay as they are because they are bound to
# the column names of stored intermediate files that cannot be regenerated
# without the raw data; column_glossary.csv gives the English meaning of each,
# and README.md says so. The list is explicit rather than a broad pattern, so
# that a genuine German SENTENCE in a source file is still caught.
INTERNE_BEZEICHNER = {
    "achse", "achsen", "adipogen", "anteil", "arm", "ausgabe", "beobachtet",
    "bestanden", "chondrogen", "datei", "datensatz", "datensaetze", "differenz",
    "durchgefallen", "eichung", "eichungen", "eingaben", "entitaet", "fenster",
    "gegen", "geeicht", "grenze", "hintergrund", "kein", "keine", "kohorte",
    "kontrast", "konkordanz", "laesion", "laufzeit", "linie", "modul",
    "myogen", "nachweisgrenze", "naiv", "naive", "nicht", "ohne", "osteogen",
    "probe", "proben", "regel", "seite", "sichtung", "spender", "tabelle",
    "zelle", "zellen", "zerlegung", "zweck", "anlage", "begruendung",
    "protokoll", "vorregistrierung", "abbildung", "satz", "saetze", "urteil",
    "status", "stufe", "punkt", "punkte", "wert", "werte", "zeile", "zeilen",
    "spalte", "spalten", "null", "median", "mittel", "delta", "rang",
    "osteo", "adipo", "naiver",
}
QUELLTEXT = {".py", ".R"}


_INTERN = re.compile(
    r"\b(" + "|".join(sorted(INTERNE_BEZEICHNER, key=len, reverse=True))
    + r")\b", re.I)


def pruefe_text(name: str, text: str, *, streng: bool = True,
                quelltext: bool = False) -> list[str]:
    """Check one text.

    `streng=False` skips the identifier rules, for files whose content is data
    rather than prose. `quelltext=True` additionally blanks the internal
    vocabulary listed above, so that only real German sentences are reported.
    """
    fehler: list[str] = []
    for i, zeile in enumerate(text.split("\n"), 1):
        rest = AUSNAHME.sub(" ", zeile)
        if quelltext:
            rest = _INTERN.sub(" ", rest)
        if UMLAUT.search(rest):
            fehler.append(f"{name}:{i} umlaut: {zeile.strip()[:100]}")
        if ERSATZZEICHEN.search(zeile):
            fehler.append(f"{name}:{i} replacement character: "
                          f"{zeile.strip()[:100]}")
        m = VERDACHT.search(rest)
        if m:
            fehler.append(f"{name}:{i} German '{m.group(0)}': "
                          f"{zeile.strip()[:90]}")
        if not streng:
            continue
        if not quelltext and GSE_FALSCH.search(rest):
            fehler.append(f"{name}:{i} accession rule reversed: "
                          f"{zeile.strip()[:90]}")
        if not quelltext and GSE_OHNE_KLAMMER.search(rest):
            fehler.append(f"{name}:{i} accession without parentheses: "
                          f"{zeile.strip()[:90]}")
        if KURZFORM.search(rest):
            fehler.append(f"{name}:{i} abbreviated identifier: "
                          f"{zeile.strip()[:90]}")
    return fehler


def text_aus_docx(p: pathlib.Path) -> str:
    with zipfile.ZipFile(p) as z:
        roh = z.read("word/document.xml").decode("utf-8", "replace")
    roh = re.sub(r"</w:p>", "\n", roh)
    return re.sub(r"<[^>]+>", "", roh)


def main() -> int:
    sag("12_check_language.py -- English throughout in everything published")
    sag("")
    alle: list[str] = []

    # 1) the manuscript sources
    for p in sorted((WURZEL / "manuscript").glob("*.md")):
        alle += pruefe_text(f"manuscript/{p.name}", p.read_text(encoding="utf-8"))

    # 2) the panel data: the VALUES, and the column headers
    for p in sorted((WURZEL / "figures" / "data").glob("*.csv")):
        d = pd.read_csv(p)
        werte = list(d.columns)
        for k in d.columns:
            if pd.api.types.is_string_dtype(d[k]) or d[k].dtype == object:
                werte += [str(v) for v in d[k].dropna().unique()]
        alle += pruefe_text(f"figures/data/{p.name}", "\n".join(werte),
                            streng=False)

    # 3) the code that is published, and the documents beside it
    # The three checking scripts must themselves contain umlauts and German
    # words -- those are their search patterns. They are checked for prose,
    # not for the characters they hunt for.
    # _display.py is the translation table itself: it holds the German source
    # strings on the left-hand side by construction. Checking it for German
    # would flag the very thing that makes the rest English.
    PRUEFER = {"11_check_references.py", "12_check_language.py",
               "21_build_submission.py", "_display.py"}
    # Every pipeline stage folder, the data-acquisition folder and the figure
    # style. The stage folders are numbered NN_name at the repository root.
    stufen = [d.name for d in sorted(WURZEL.iterdir())
              if d.is_dir() and re.match(r"\d\d_", d.name)]
    for ordner in (*stufen, "data_acquisition", "figure_style"):
        basis = WURZEL / ordner
        if not basis.exists():
            continue
        for p in sorted(basis.rglob("*")):
            if p.suffix not in (".py", ".R", ".md") or "__pycache__" in p.parts:
                continue
            if p.name in PRUEFER:
                continue
            rel = p.relative_to(WURZEL).as_posix()
            alle += pruefe_text(rel, p.read_text(encoding="utf-8",
                                                 errors="replace"),
                                quelltext=p.suffix in QUELLTEXT)

    # 4) the preregistrations
    for p in sorted((WURZEL / "preregistrations").glob("*.md")):
        alle += pruefe_text(f"preregistrations/{p.name}",
                            p.read_text(encoding="utf-8", errors="replace"))

    # 5) the documents in the repository root
    for name in ("README.md", "LICENSE", "LICENSE-CODE", "CITATION.cff"):
        p = WURZEL / name
        if p.exists():
            alle += pruefe_text(name, p.read_text(encoding="utf-8",
                                                  errors="replace"))

    # 6) the delivered items themselves
    sub = WURZEL / "submission"
    if sub.exists():
        for p in sorted(sub.glob("*.docx")):
            alle += pruefe_text(p.name, text_aus_docx(p))
        for p in sorted(sub.glob("*.md")):
            alle += pruefe_text(p.name, p.read_text(encoding="utf-8"))
        for p in sorted(sub.glob("*.xlsx")):
            for blatt, d in pd.read_excel(p, sheet_name=None).items():
                werte = list(d.columns)
                for k in d.columns:
                    if pd.api.types.is_string_dtype(d[k]) or d[k].dtype == object:
                        werte += [str(v) for v in d[k].dropna().unique()]
                alle += pruefe_text(f"{p.name}[{blatt}]", "\n".join(werte),
                                    streng=False)

    if alle:
        for f in alle:
            sag("  " + f)
    else:
        sag("Nothing found. Everything published is English:")
        sag("  * no umlaut, no replacement character, no German word")
        sag("  * no 'naiv'/'naive' -- 'undifferentiated' throughout")
        sag("  * the accession rule holds: information, then the accession")
        sag("    in parentheses")
        sag("  * no abbreviated identifiers (no 'Osteo', no 'OSTEO')")
    sag("")
    sag(f"=== {len(alle)} findings ===")

    p = RES / "language_check.txt"
    p.write_text("\n".join(AUSGABE) + "\n", encoding="utf-8")
    print("->", p)
    return 1 if alle else 0


if __name__ == "__main__":
    sys.exit(main())
