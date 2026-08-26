> Translated from the German original of 2026-08-18. The content, the dates
> and every number are unchanged.

# Preregistration S1 — an induction-matched null for `FN1`

Written 2026-08-18, **before the first computation**. A stopping criterion of
the project. Once fixed, nothing in this file is changed; deviations are added
as a dated addendum at the end.

## 1. The question

In the predecessor project the `LAMA5` knockout and the `FN1` patient lines lie
in opposite quadrants:

| data set | cargo−machinery decoupling | z_corrected | blunting |
|---|---|---|---|
| LAMA5-KO osteogenic | +0.408 | +9.86 | −0.463 |
| LAMA5-KO chondrogenic | +0.107 | +3.65 | −0.247 |
| FN1 C123R | −0.384 | −5.28 | **−1.022** |
| FN1 C231W | −0.441 | −7.05 | **−1.172** |
| SERPINA3-KD | −0.210 | −1.30 | −0.427 |
| MIR181A1HG-KD | +0.000 | −0.56 | −0.663 |

The preregistered criterion K3 has fallen: the FN1 lines are roughly four times
more strongly blunted than our own system. **The question: does the sign
reversal follow mechanically from that severity, or does it survive having the
severity removed?**

## 2. What the question is not

A purely **additive** severity effect cannot produce a sign reversal in a
contrast between two gene sets — it cancels out. The reversal requires the
blunting to act **differentially** on cargo and machinery. The plausible route
to that is induction dependence: strongly induced genes are blunted more
strongly (a ceiling or regression effect), and cargo and machinery differ in
their induction strength. That is exactly what T1 and T2 address. The existing
null matches induction in only **four** quartile classes — coarse enough that
at fourfold blunting the residual differences within a class could carry the
effect.

## 3. The three tests

The metric, the sets (`S_FRACHT0`, `S_MASCHINE`, `S_DISTAL`), the VIF
correction and the pooling are unchanged. Only the null, or the target
quantity, is changed.

- **T1 — finer induction matching.** Twenty vigintile classes of `dWT` instead
  of four quartile classes; in addition, a nearest-neighbour matching on `dWT`
  as a variant. If z shrinks towards 0, the reversal was a residual induction
  artefact.
- **T2 — induction-residualised interaction (the decisive test).** A loess fit
  `iv ~ dWT` over the pool, with the contrast computed on the residuals. That
  removes **every** monotone dependence of the interaction on induction
  strength — and a ceiling or severity effect is exactly that.
- **T3 — severity-normalised scale.** The contrast divided by the spread of
  `iv` in the pool, for all six data sets. This clarifies the **comparability
  of the magnitudes** between the quadrants, not the sign (division by a
  positive number cannot turn a sign).

MDE80 is reported alongside every null. All tests run on **all six** data sets,
not only on FN1 — LAMA5 is the positive control, SERPINA3 and MIR181A1HG are
the negative controls **within the same procedure**.

## 4. Decision rule

**PASSED** — only if all three hold:
1. T2: in **both** FN1 lines the decoupling stays negative with
   z_corrected < −2;
2. T1 (vigintiles): in both FN1 lines z_corrected < −2;
3. T2: LAMA5 osteogenic stays positive with z_corrected > +2.

→ The typology carries. Continue with the positive control.

**FALLEN** — if (1) or (2) is missed, that is, if in one FN1 line the
residualised or finely matched z rises above −2 or turns sign.

→ The two-quadrant statement is not tenable. **The project ends here.** The
result is recorded as a methods note; no reformulation, no recomputation with
different sets.

**PROCEDURE INVALID** — if condition (3) is missed **while** the FN1 lines
collapse as well. The residualisation has then removed the signal in general
and says nothing about FN1. This is not a stop but a correction of the
procedure; the case is documented explicitly and is **not** read as a pass.

**GREY ZONE** — one FN1 line passes, the other does not. The typology then
counts as **not established**: a methods report, no publication on that basis.

## 5. What would refute the result

The typology is refuted if the sign reversal between LAMA5 and FN1 can be
represented as a function of perturbation severity — concretely, if after T2
both systems lie on the same side of the null or FN1 loses significance.
Conversely: if the reversal survives the removal of the induction dependence,
the severity objection is answered and the capacity-against-cargo axis is not
replaceable by a severity grade.

## 6. Known limits that this test does NOT remove

- `GSE251698` is patient material (iPSC, two patients against **one** control),
  not an isogenic design. The comparison with the isogenic LAMA5 system stays
  asymmetric in that respect.
- The comparison is arm-matched only against **LAMA5 chondrogenic** (+0.107).
  The osteogenic core finding (+0.408) stands beside it across arms and is not
  used as a reference quantity.
- Only two data sets per quadrant. The test settles the severity objection, not
  the question whether these are two **types** or two **individual cases**.

---

## Addendum 2026-08-18 (before the computation, after the rule was released)

- **The strictest version was chosen** and is binding: T2 decides, both FN1
  lines must hold, and the grey zone counts as *not established*.
- **Implementation of T1:** the nearest-neighbour matching is realised as class
  formation with **50** classes rather than as caliper matching — at pool sizes
  of 1 347 to 12 039 genes and set fractions around 10 % that is the same
  resolution at considerably lower computational cost. The ladder reported is
  **4 / 20 / 50** classes; 4 is the previous project standard and serves as a
  reproduction control.
- Class counts above 50 are not computed: there the draw approaches
  determinism, the spread of the null shrinks artificially and z becomes
  anticonservative.
