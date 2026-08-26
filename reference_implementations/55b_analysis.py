# -*- coding: utf-8 -*-
"""
55b_analysis.py -- phase M-E: the three tests of the preregistration.

Preregistration: `preregistrations/PRAEREG_M_E.md`, dated before the first
statistic. Primary is **Test 1**.

  Test 1  noise null floor: does `iv_1x1` converge more strongly than a
          lesion-free `pseudo_iv_1x1` of the same algebra and sample count?
  Test 2  SNR equalization: does `dWT` still converge when its SHARED share
          is compressed to the signal-to-noise ratio of `iv`?
  Test 3  cross-control set vs quantity (descriptive, decides nothing).

Convergence rule unchanged from the main part: a gene converges if its sign
is the same in >= 90 % of the datasets in which it is measurable.
Universe: genes measurable in >= 16 of 18 datasets.

Output: derived_data/M_kalibrierung/{test1_*.csv, test2_snr.csv,
         test3_kreuz.csv, 55_log.txt}
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _module import ERGEBNISSE, MODUL, TABELLEN, leiter  # noqa: E402

WURZEL = pathlib.Path(__file__).resolve().parents[1]
ANTRAG = WURZEL.parent
AUS = ERGEBNISSE / "M_kalibrierung"
AUS.mkdir(parents=True, exist_ok=True)
KONTRASTE = AUS / "kontraste"
GENE20D = (ANTRAG / "backups" / "_backup_2026-08-19_vor_paperaufbau" /
           "20_Exploration" / "derived_data")

MIN_DATENSAETZE = 16       # universe U1, fixed in advance
KONSISTENZ = 0.90          # convergence rule of the main part
SEED_E = 20260823
N_FLIP = 2000

LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


# ---------------------------------------------------------------------------
def konvergent(S: np.ndarray) -> np.ndarray:
    """Convergence count per column of a (genes x datasets x draws) matrix
    of signs. 0 counts as not measurable."""
    messbar = (S != 0).sum(axis=1)
    pos = (S > 0).sum(axis=1)
    neg = (S < 0).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        anteil = np.maximum(pos, neg) / np.where(messbar > 0, messbar, np.nan)
    return ((anteil >= KONSISTENZ) & (messbar >= MIN_DATENSAETZE)).sum(axis=0)


def lade_vorzeichen(groesse: str, universum: pd.Index) -> np.ndarray:
    """-> (genes x datasets x draws), int8."""
    dateien = sorted(KONTRASTE.glob(f"{groesse}_*.csv.gz"))
    n_d = len(dateien)
    S = None
    for j, f in enumerate(dateien):
        d = pd.read_csv(f, index_col=0)
        d = d.reindex(universum).fillna(0)
        if S is None:
            S = np.zeros((len(universum), n_d, d.shape[1]), dtype=np.int8)
        S[:, j, :] = d.to_numpy(dtype=np.int8)
    return S


# ---------------------------------------------------------------------------
def test1() -> pd.DataFrame:
    log("=" * 78)
    log("TEST 1 (primary) -- the noise null floor for iv-shaped quantities")
    log("=" * 78)

    zaehl: dict[str, int] = {}
    for f in sorted(KONTRASTE.glob("iv_*.csv.gz")):
        for g in pd.read_csv(f, usecols=[0], index_col=0).index:
            zaehl[g] = zaehl.get(g, 0) + 1
    universum = pd.Index(sorted(g for g, k in zaehl.items()
                                if k >= MIN_DATENSAETZE))
    log("Universe: %d genes in >= %d of 18 datasets"
        % (len(universum), MIN_DATENSAETZE))

    zeilen, ziehungen = [], {}
    for groesse, name in (("iv", "iv_1x1"), ("pseudo", "pseudo_iv_1x1"),
                          ("dwt", "dWT_1x1")):
        S = lade_vorzeichen(groesse, universum)
        k = konvergent(S)
        ziehungen[name] = k
        zeilen.append(dict(
            groesse=name, n_ziehungen=len(k), mittel=float(k.mean()),
            sd=float(k.std(ddof=1)), median=float(np.median(k)),
            p2_5=float(np.percentile(k, 2.5)),
            p97_5=float(np.percentile(k, 97.5)),
            min=int(k.min()), max=int(k.max())))
        log("  %-14s mean %8.1f  SD %6.1f  [2.5 %% %6.1f .. 97.5 %% %7.1f]"
            % (name, k.mean(), k.std(ddof=1),
               np.percentile(k, 2.5), np.percentile(k, 97.5)))
        del S

    T = pd.DataFrame(zeilen)
    T.to_csv(AUS / "test1_rauschboden.csv", index=False)
    pd.DataFrame(ziehungen).to_csv(AUS / "test1_ziehungen.csv", index=False)

    iv_m = T.loc[T.groesse == "iv_1x1", "mittel"].iloc[0]
    ps_hi = T.loc[T.groesse == "pseudo_iv_1x1", "p97_5"].iloc[0]
    ps_m = T.loc[T.groesse == "pseudo_iv_1x1", "mittel"].iloc[0]
    log("")
    log("  Decision rule: mean(iv_1x1) %.1f against the 97.5th percentile("
        "pseudo_iv_1x1) %.1f" % (iv_m, ps_hi))
    if iv_m <= ps_hi:
        log("  -> iv sits AT the noise floor of its own construction."
            "  THREAD B.")
    else:
        log("  -> iv carries structure above the noise floor.  continue with test 2.")
    log("  (ratio iv/pseudo on average: %.2f)" % (iv_m / ps_m))
    return T


# ---------------------------------------------------------------------------
def lade20d() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Gene x point for dWT and iv from the existing tables."""
    teile = [pd.read_csv(f) for f in sorted(GENE20D.glob("20d_gene_*.csv"))]
    G = pd.concat(teile, ignore_index=True)
    dwt = G.pivot_table(index="gen", columns="punkt", values="dWT")
    iv = G.pivot_table(index="gen", columns="punkt", values="iv")
    ok = dwt.notna().sum(axis=1) >= MIN_DATENSAETZE
    return dwt[ok], iv.reindex(dwt[ok].index)


def zaehle(M: pd.DataFrame) -> int:
    S = np.sign(M.to_numpy(dtype=float))
    S = np.nan_to_num(S, nan=0.0)[:, :, None]
    return int(konvergent(S.astype(np.int8))[0])


def flip_null(M: pd.DataFrame, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    A = np.nan_to_num(np.sign(M.to_numpy(dtype=float)), nan=0.0).astype(np.int8)
    aus = np.empty(n)
    for i in range(n):
        f = rng.choice(np.array([-1, 1], dtype=np.int8), size=A.shape[1])
        aus[i] = konvergent((A * f[None, :])[:, :, None])[0]
    return aus


def test2() -> pd.DataFrame:
    log("")
    log("=" * 78)
    log("TEST 2 -- the SNR equalization")
    log("=" * 78)
    dwt, iv = lade20d()
    log("Universe: %d genes in >= %d of 18 datasets" % (len(dwt),
                                                         MIN_DATENSAETZE))

    def snr(M: pd.DataFrame) -> tuple[float, pd.Series, pd.Series]:
        m = M.mean(axis=1)
        s = M.std(axis=1, ddof=1)
        return float((m.abs() / s).median()), m, s

    snr_d, m_d, s_d = snr(dwt)
    snr_i, _, _ = snr(iv)
    k = snr_i / snr_d
    log("  SNR(dWT) %.4f | SNR(iv) %.4f | compression factor k = %.4f"
        % (snr_d, snr_i, k))

    dwt_stern = dwt.sub(m_d, axis=0).add(k * m_d, axis=0)

    zeilen = []
    for name, M in (("dWT (beobachtet)", dwt), ("iv (beobachtet)", iv),
                    ("dWT* (SNR-angeglichen)", dwt_stern)):
        n = zaehle(M)
        null = flip_null(M, N_FLIP, SEED_E)
        zeilen.append(dict(groesse=name, konvergent=n,
                           null_mittel=float(null.mean()),
                           null_sd=float(null.std(ddof=1)),
                           null_p97_5=float(np.percentile(null, 97.5)),
                           ueber_null=bool(n > np.percentile(null, 97.5))))
        log("  %-24s convergent %5d | flip null %7.1f +- %5.1f "
            "(97.5 %% %6.1f) | %s"
            % (name, n, null.mean(), null.std(ddof=1),
               np.percentile(null, 97.5),
               "above the null" if n > np.percentile(null, 97.5) else "on the null"))

    T = pd.DataFrame(zeilen)
    T.to_csv(AUS / "test2_snr.csv", index=False)

    n_stern = T.loc[T.groesse == "dWT* (SNR-angeglichen)", "konvergent"].iloc[0]
    hi = T.loc[T.groesse == "dWT* (SNR-angeglichen)", "null_p97_5"].iloc[0]
    n_iv = T.loc[T.groesse == "iv (beobachtet)", "konvergent"].iloc[0]
    log("")
    log("  Decision rule: dWT* %d against null-97.5 %% %.1f and against "
        "2 x iv = %d" % (n_stern, hi, 2 * n_iv))
    if n_stern > hi and n_stern >= 2 * n_iv:
        log("  -> the asymmetry survives the equalization.  THREAD A.")
    else:
        log("  -> the asymmetry cannot be separated from size and noise."
            "  THREAD B.")
    return T


# ---------------------------------------------------------------------------
def test3() -> pd.DataFrame:
    log("")
    log("=" * 78)
    log("TEST 3 -- cross-control set vs quantity (descriptive)")
    log("=" * 78)
    dwt, iv = lade20d()
    prog = set(MODUL.ensembl)
    laes = set(pd.read_csv(ERGEBNISSE / "M_patienten" /
                           "laesionssatz_173.csv").ensembl)
    zeilen = []
    for gname, gset in (("Programmsatz (173)", prog),
                        ("Laesionssatz (173)", laes)):
        idx = dwt.index.intersection(sorted(gset))
        for qname, M in (("dWT", dwt), ("iv", iv)):
            n = zaehle(M.loc[idx])
            null = flip_null(M.loc[idx], 500, SEED_E)
            zeilen.append(dict(genset=gname, groesse=qname, n_gene=len(idx),
                               konvergent=n, null_mittel=float(null.mean()),
                               null_p97_5=float(np.percentile(null, 97.5))))
            log("  %-20s %-4s  n %3d | convergent %3d | null %5.1f "
                "(97.5 %% %5.1f)"
                % (gname, qname, len(idx), n, null.mean(),
                   np.percentile(null, 97.5)))
    T = pd.DataFrame(zeilen)
    T.to_csv(AUS / "test3_kreuz.csv", index=False)
    return T


# ---------------------------------------------------------------------------
def main() -> None:
    log("Phase M-E -- PRAEREG_M_E.md, seed %d" % SEED_E)
    test1()
    test2()
    test3()
    (AUS / "55_log.txt").write_text("\n".join(LOG), encoding="utf-8")
    log("")
    log("-> %s" % AUS)


if __name__ == "__main__":
    main()
