# *Rhodopseudomonas palustris* CGA009 ME-model — pure-Python port

A standalone **Python** implementation of the *R. palustris* CGA009 genome-scale
**Metabolism & Expression (ME)** model. It runs with only `numpy` + `scipy` — **no GAMS, no
CPLEX, no cobra, no license** — and **exactly reproduces the published Saha-lab growth rates**
for all four substrates.

| Substrate | Published µ (Saha Table 2) | This port | ATP maintenance |
|---|---|---|---|
| p-coumarate | 1.21 | **1.21** ✓ | 85.4 |
| acetate | 0.77 | **0.77** ✓ | 54.0 |
| butyrate | 0.86 | **0.86** ✓ | 56.7 |
| succinate | 0.74 | **0.74** ✓ | 45.7 |

"Reproduces exactly" = the ME LP is **feasible at the published µ and infeasible 0.02 above**,
for every substrate (see `reproduce_saha_results.ipynb`).

## Attribution
The **model itself** is the work of the Saha lab:
> Chowdhury N.B., Alsiyabi A., Saha R. *Characterizing the Interplay of Rubisco and
> Nitrogenase Enzymes in Anaerobic-Photoheterotrophically Grown Rhodopseudomonas palustris
> CGA009 through a Genome-Scale Metabolic and Expression Model.* **Microbiology Spectrum** (2022).
> Original code (GAMS + a fixed-µ Pyomo export): https://github.com/ssbio/palustris_ME_model

This repository is **only a solver-independent Python re-implementation** of that model,
built from the original GAMS data files (redistributed here under `data/` for reproducibility).
All modeling credit belongs to the original authors.

## License
The Python re-implementation and analysis code here are released under the [MIT License](LICENSE)
© 2026 Jackson Comes. The Saha-lab model files (`data/`, `validation/gams_reference/`) are **not**
covered by that grant. As of July 2026 the upstream repository
([ssbio/palustris_ME_model](https://github.com/ssbio/palustris_ME_model)) carries **no license**, so
its authors have not formally granted redistribution rights; those files are included here in good
faith for scientific reproducibility, with citation. If you use this port, please cite the Saha lab's
original model (above); if you plan to redistribute or build on their files beyond forking this repo,
contact the original authors. See the note at the bottom of [`LICENSE`](LICENSE).

## Why this exists
The published model is written in **GAMS** and its growth-maximization is a **nonlinear program**.
The authors' workflow requires running GAMS (or a commercial solver) once per growth rate, which
makes flux-variability analysis and parameter sweeps tedious. This port removes that dependency:
the ME problem is only nonlinear because of growth-rate (µ) coupling, so **fixing µ makes it a
plain linear program** that any Python LP solver handles. A bisection over µ then recovers the
maximum growth rate.

## Install & run
```bash
pip install -r requirements.txt          # numpy, scipy  (optional: jupyter for the notebook)
python palustris_me.py                    # loads the model and prints the Table-2 reproduction
```
```python
from palustris_me import MEModel
me = MEModel()                            # load once (~0.5 s)
up, atp, saha_mu = me.set_medium("coumarate")
status, obj, fluxes = me.solve(1.14, atp_maint=atp)     # solve at a fixed growth rate
mu_max = me.bisect_max_mu(atp_maint=atp)                 # ME-predicted max growth
print(me.flux(fluxes, "R_rxn00018_c0_2"))                # e.g. RuBisCO/CO2-fixation flux
```
FVA / knockouts / medium sweeps: just call `solve()` in a loop at a fixed µ — every call is a
fast LP, no GAMS.

## What the model is (and how the LP is built)
At a fixed growth rate µ the model is exactly the GAMS LP `palustris_ME.gms`:

```
minimize   Σ v_j                                  # pFBA (primalobj)
s.t.       S(µ) · v = 0                            # mass balance with µ-coupling
           v[bio2] = µ                             # growth pinned to µ
           Σ_{ATP-producing} v = atp_maint         # ATP maintenance (atps_const2; substrate-specific)
           Σ_{ATP-producing} v = v10042+v01517+v00148B1+v00148B2   # atps_const
           v[R_rxn10042_c0] = v[R_rxn37614_c0]     # ATP synthase = photosystem (atps_eq_ps2)
           v[R_rxn00288_c0] = 0                     # succTOfum
           lb ≤ v ≤ ub  + fixed/stopped reactions
```
`S(µ)` is the stoichiometric matrix with **~8,500 growth-rate-dependent coupling coefficients**
that encode enzyme dilution and the transcription/translation/ribosome machinery cost. They are
read symbolically from the `.gms` and re-evaluated at each µ.

## How this port was built from the GAMS model (the "edits")
Nothing about the biology was changed. The port is a faithful translation; the only work was
extracting the model from GAMS syntax into matrices a Python LP solver can consume:

1. **Base stoichiometry** — parsed from `data/sij.txt` (`'met'.'rxn' coef`).
2. **Growth-rate coupling (the crucial part)** — the ~8,500 `S('…','…') = <expr in µ>;`
   assignments in `data/palustris_ME.gms` are parsed *with their full right-hand-side
   expression* (handling multi-line statements and all six algebraic forms, e.g.
   `-mu/kcat`, `-mu/(a+b+…)`, `-N*(mu+(kt*r0))/(kt*cribo2/X)`), the E. coli scalars
   (`kcat, kt, r0, cmrna2, cribo2`) are substituted, and each expression is compiled once and
   re-evaluated at every µ. These overwrite the placeholder values in `sij.txt`.
   *(An earlier version that captured only the `-mu/…` enzyme-dilution forms and missed the
   transcription/translation coupling over-predicted growth — that bug is fixed and guarded by
   the validation.)*
3. **Custom constraints** — the four non-mass-balance equations (`atps_const`, `atps_const2`,
   `atps_eq_ps2`, `succTOfum`) are added as explicit rows.
4. **Bounds & stopped reactions** — from `data/lower_bound.txt` / `data/upper_bound.txt`, plus
   the `v.up/v.lo = 0` fixed reactions in the `.gms`.
5. **Growth handling** — `v[bio2]` is pinned to µ (the GAMS `v.lo/v.up('bio2')=mu`); max growth
   is found by bisection (feasible → raise µ; infeasible → lower µ), matching the authors'
   "increase µ until infeasible" procedure.
6. **Substrate conditions** — the four Table-2 media (carbon-source uptake + substrate-specific
   ATP maintenance) are in `SUBSTRATES` in `palustris_me.py`.

## Solver & numerics
Solved with **HiGHS** via `scipy.optimize.linprog(method="highs")`. ME-models are ill-scaled
(coupling coefficients ~1e-6 alongside O(1) stoichiometry), which can make weaker double-precision
solvers falsely report infeasibility near the growth boundary. HiGHS handles this cleanly here
(sharp feasible→infeasible transition at each substrate's published µ). If a future variant
misbehaves at the boundary, escalate precision (`glpk_exact`, SoPlex iterative refinement, or the
quad-precision qMINOS the authors' solveME stack uses).

## Caveat inherited from the model
Only ~64 enzymes carry *R. palustris*-specific efficiencies; the remainder (including the "average"
enzyme) use the **E. coli** `kcat = 234000` and E. coli translation/transcription constants. The ME
framework is therefore rigorous in *form* (and gives a working nitrogenase / N₂-fixing, H₂-evolving
regime), but its per-enzyme kinetic parameters are largely borrowed. This is a property of the
original model, unchanged here.

## Files
```
palustris_me.py                 the model (load, solve, bisect_max_mu, set_medium)
reproduce_saha_results.ipynb    reproduces all four Table-2 growth rates + a flux inspection
validate.py                     command-line version of the reproduction test
requirements.txt                numpy, scipy
data/                           original Saha-lab GAMS files (sij, reactions, metabolites,
                                bounds, palustris_ME.gms) — redistributed for reproducibility
tests/                          pytest suite (test_me_audit.py) + pytest.ini
validation/                     audit + GAMS cross-check
    AUDIT_REPORT.md, CHANGES.md
    gams_reference/             authors' GAMS-Convert instances + compare_to_gams.py
```
Run the tests with a bare `pytest` at the repo root.
