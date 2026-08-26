# -*- coding: utf-8 -*-
"""
20_in_vivo_donor_test.py -- donor-stratified trend test of the in vivo axis

Purpose   The in vivo module trend has to be tested with the specimen as the
          unit of analysis, not the sample: a difference between specimen
          means must not be able to produce a trend. This script computes
          that version from the per-sample contrast values that already
          exist, without loading the atlas again.

Unit      The specimen (`adj_sample` prefix Pcw<stage>[_s<i>]). The atlas
          carries 9 developmental stages across 16 specimens.

Null      Zone ranks are permuted WITHIN the specimen. That is the stricter
          test compared with a free permutation of samples, and it is the
          one reported in the manuscript.

Inputs    derived_data/followup/ws4_modulwert_je_probe.csv
          derived_data/followup/ws4_positivkontrolle_je_probe.csv
          derived_data/followup/ws4_p3_panel_vs_modul.csv
Outputs   results/invivo_spendertest.csv
Runtime   a few seconds
"""
import os, pathlib, numpy as np, pandas as pd
from scipy.stats import spearmanr

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
ERG = WURZEL / "derived_data" / "followup"
AUS = WURZEL / "results"; AUS.mkdir(parents=True, exist_ok=True)

ZONEN = ["MesCond", "ChondroProg", "RestingChon", "ProlifChon",
         "PrehyperChon", "HyperChon"]
RANG = {z: i + 1 for i, z in enumerate(ZONEN)}
NZIEH, SEED = 20000, 20260823


def spezimen(probe: pd.Series) -> pd.Series:
    return probe.str.extract(r"^(Pcw[0-9.]+(?:_s[0-9])?)", expand=False)


def test(df: pd.DataFrame, label: str, rang=None,
         spezimen_fn=None) -> dict:
    """Donor-stratified trend test.

    `rang` and `spezimen_fn` are parameters so that
    code/33_postnatal_growth_plate_test.py can use the same implementation
    on the postnatal growth plate -- the project rule is that no metric is
    implemented twice. Called without arguments the behaviour is identical
    to the run stored under results/.
    """
    d = df.copy()
    d["rang"] = d.zone.map(RANG if rang is None else rang)
    fn = spezimen if spezimen_fn is None else spezimen_fn
    d["spz"] = fn(d.probe)
    obs = spearmanr(d.rang, d.kontrast).statistic
    rng = np.random.default_rng(SEED)
    gruppen = [g.index.to_numpy() for _, g in d.groupby("spz")]
    rangwerte = d["rang"].to_numpy()
    null = np.empty(NZIEH)
    idx = d.index.to_numpy()
    pos = {v: i for i, v in enumerate(idx)}
    grp_pos = [np.array([pos[i] for i in g]) for g in gruppen]
    kontrast = d["kontrast"].to_numpy()
    for b in range(NZIEH):
        r = rangwerte.copy()
        for gp in grp_pos:
            r[gp] = rng.permutation(r[gp])
        null[b] = spearmanr(r, kontrast).statistic
    mu, sd = float(null.mean()), float(null.std(ddof=1))
    z = (obs - mu) / sd
    # MDE80: the smallest statistic that would exceed the 95th percentile of
    # the null in 80 % of repetitions -- the same rule as everywhere else
    mde80 = float(np.quantile(null, 0.95) + 0.8416 * sd)
    p = float((np.sum(null >= obs) + 1) / (NZIEH + 1))
    return {"groesse": label, "n_proben": len(d), "n_spezimen": d.spz.nunique(),
            "rho": obs, "null_mittel": mu, "null_sd": sd, "z": z,
            "p_perm_oben": p, "mde80_rho": mde80,
            "ueber_mde80": bool(obs > mde80), "nzieh": NZIEH, "seed": SEED}


def main() -> None:
    mw = pd.read_csv(ERG / "ws4_modulwert_je_probe.csv")
    pk = pd.read_csv(ERG / "ws4_positivkontrolle_je_probe.csv")
    pk = pk[pk.vergleich == "chondrogen_vs_naiv"]
    p3 = pd.read_csv(ERG / "ws4_p3_panel_vs_modul.csv")
    zeilen = [test(pk, "Positivkontrolle chondrogen-naiv"),
              test(mw, "Modul (173 Gene)"),
              test(p3, "PA309 gegen Modul")]
    out = pd.DataFrame(zeilen)
    out.to_csv(AUS / "invivo_spendertest.csv", index=False)
    print(out.to_string(index=False))
    print(f"\nnumpy {np.__version__} | pandas {pd.__version__}")


if __name__ == "__main__":
    main()
