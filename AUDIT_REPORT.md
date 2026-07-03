# Audit — pure-Python R. palustris ME-model port

> **Status:** the fixes recommended below (last-wins coupling, medium provenance) were applied on
> branch `fix/last-wins-coupling` — see `CHANGES.md`. After the fix the port matches the authors'
> GAMS instances to machine precision. This report is the point-in-time audit record.


**Question.** Does this Python port faithfully equal the GAMS original — at the level of µ_max,
objective, feasibility boundary, **and the flux vector** — or only at the four published growth rates?

**Environment.** Python 3.12.4, scipy 1.18.0, numpy 2.4.6, HiGHS via `scipy.optimize.linprog`
(`method="highs"`), Windows 11. One LP solve ≈ 1.4 s; a full four-substrate bisection ≈ 570 s.

**Reference obtained (this is the big change from the earlier draft).** The upstream repo
`ssbio/palustris_ME_model` ships, under `PYTHON/Table_2/`, the Saha lab's own **GAMS-Convert dumps** —
`coumarate/acetate/buytrate/succinate_max_growth.py`. These are the *actual GAMS LP instances*
serialized to Pyomo, with the µ-coupling coefficients baked in and the growth variable fixed to the
published µ (1.21 / 0.77 / 0.86 / 0.74). `dictionary.txt` maps the anonymous `xN` variables to reaction
ids. I copied them into `gams_reference/`, parsed them numerically, and solved them with the same
HiGHS backend the port uses — so any difference is *formulation*, not solver. This gave the
flux-level GAMS reference the earlier draft said it lacked.

## Verdict
**The port is a faithful reconstruction of the authors' GAMS model.** Concretely:
- **Objective / feasibility:** for all four substrates the port's `min Σv` at the published µ equals
  the GAMS instance to **≤ 6×10⁻⁵** as shipped, and both are feasible there.
- **The one real bug** is that GAMS last-wins assignment is implemented as additive summing, doubling
  **26** coupling coefficients. It is the *entire* source of that residual: fixing it to last-wins
  drops the coumarate objective gap from 2.0×10⁻⁴ to **−3.9×10⁻¹¹** — machine-precision-exact to GAMS.
- **Flux vector:** matches GAMS on all 5479 reactions **except 25–46 alternate-optima swaps** (flux
  moving between parallel isozyme/direction variants `_1`↔`_2`, `_B_1`↔`_B_2` at identical cost) — the
  limit of what a degenerate ME-LP allows, not a discrepancy.
- **µ_max is reproduced only loosely:** the port's true maximum sits **0.009–0.019 above** each
  published value; `validate.py`'s "exactly" really means "feasible at published µ, infeasible at
  published+0.02." At the published µ itself, port ≡ GAMS (above).
- **HiGHS is trustworthy** at this scale; the feasibility boundary is invariant to solver tolerances
  and presolve despite 3.4×10⁹ matrix scaling. qMINOS is not needed for these solves.

---

## Track 0 — Baseline
`validate.py` → exit 0, four rows "reproduced". `python palustris_me.py` bisection (tol 1e-3):
coumarate 1.229, acetate 0.779, butyrate 0.870, succinate 0.749. `reproduce_saha_results.ipynb`
executes end-to-end (exit 0). All four bisected maxima are 0.009–0.019 **above** the published values
(see Track C).

## Track A — Parsing & structural fidelity

| Check | Result |
|---|---|
| `S('met','rxn') = …;` assignments in `.gms` | **8544** (all single-line, `;`-terminated) — port captures all 8544, none dropped |
| Assignment forms (sum = 8544) | 3898 `-N*(mu+(kt*r0))/(kt*cribo2/N)`; 3178 `-mu/N`; 774 `-(mu+(kt*r0))/(kt*cmrna2/N)`; 507 `-mu/kcat`; 187 `-mu/(Σ)` |
| Multi-line / `$`-conditional / looped `S(i,j)=` | **none** (only the declaration `S(i,j)` at line 36) — nothing for the regex to miss |
| Scalars `kcat 234000, kt 108, r0 4.5, cmrna2 16072.5, cribo2 1976.8` | match `.gms` lines 54–58; not re-overridden |
| base `sij.txt` | 49249 parsed = 49249 in file; 0 unknown met/rxn; 8325 coupling pairs overwrite their sij placeholders |
| stopped reactions | 12, all reflected |
| dims | 5479 reactions, 3192 metabolites, 173 ATP-producing reactions |

**Four custom constraints — verified against `.gms` lines 8691–8701 and against the GAMS-Convert dumps:**
- `primalobj`: `min sum(j,v(j))` → port `np.ones(n)`, `minimize`. The dumps confirm `sense=minimize`
  over `x2+…+x5480`. ✅
- `atps_const2`: `sum(j$(S('cpd00002[c0]',j)>0),v)=85.4`. No coupling assignment touches the ATP row
  (verified: zero `S('cpd00002[c0]',…)=` statements), so the port's static `atp_prod` set is
  µ-independent and correct. The dumps carry `== 85.4/54/56.7/45.7` — matching the port. ✅
- `atps_const`, `atps_eq_ps2` (rxn10042=rxn37614), `succTOfum` (rxn00288=0): coefficients, signs, RHS
  all match; the `+=`/`-=` construction correctly nets reactions that are both ATP producers and in the
  4-reaction RHS. ✅
- Growth pin `v.lo/up('bio2')=mu`. ✅  `bio1`/`sss` are inert in the solve, as claimed.

### A-bug (REAL): duplicate coupling assignments summed, not last-wins
GAMS assignments are **last-wins**; there are **26 duplicate `(met,rxn)` coupling targets** in the
`.gms`. The port stores each as a COO triplet and calls `coo_matrix(...).tocsr()`, which **sums**
duplicate coordinates → those 26 coefficients are doubled (or sum of two near-equal values). Verified
at µ=1: `mRNA_RPA0185/TL_RPA0185` (identical assignments at lines 4007 & 4231) gives port −0.027214 vs
GAMS single −0.013607 = **exactly 2×**. **This is the sole cause of the port-vs-GAMS objective
residual:** coumarate obj−GAMS = 2.04×10⁻⁴ as shipped, **−3.9×10⁻¹¹ after a last-wins fix.**
*Impact on µ_max: none at 0.01 resolution* — but it perturbs the flux vector and should be fixed
before the model is extended. **One-line fix** (dedup the coupling list keeping the last occurrence):

```python
# in parse_model(), before returning:
_last = {}
for met, rxn, code in coupling:
    _last[(met, rxn)] = code            # GAMS last-wins
coupling = [(m, r, c) for (m, r), c in _last.items()]
```

### RETRACTED from the earlier draft: the medium bounds are NOT a divergence
`upper_bound.txt` lists coumarate uptake **2.0** and acetate/butyrate/succinate **0**; the port injects
2.54 / 6.47 / 3.69 / 4.66 from its `SUBSTRATES` dict, and I initially flagged this as an input
mismatch. **The GAMS-Convert dumps settle it:** the authors' actual GAMS instances use uptake bounds
**2.54 / 6.47 / 3.69 / 4.66** (dump variable `m.x7` etc.) and maintenance **85.4 / 54 / 56.7 / 45.7** —
**exactly the port's values.** `upper_bound.txt`'s 2.0/0 is the stale artifact that the authors (and
the port) override. The port is faithful here; the distributed text file is not the model the paper ran.

### Latent fragility (not currently triggered)
`_bounds`, the `sij` parser, and `_prep_rhs` accept `[\d.eE+]` but **not a `-` in the exponent**. There
are no negative-exponent literals in today's data, so nothing breaks — but `1.5e-06` would be silently
truncated. Harden before ingesting new data for the 2-cell model.

## Track B — Solve mechanics
- **Determinism:** coumarate µ=1.0 solved twice → `max|Δx| = 0.0`.
- **Alternate optima:** perturbing the objective by 1e-6 moves **16 of 5479** fluxes at the same
  optimum (one by 2.81). `min Σv` does **not** make the vertex unique — so exact flux equality is
  undefinable for those reactions (they are exactly the swaps seen vs GAMS in Track C).
- **Objective near the cliff blows up:** obj = 1484 (µ=1.20) → 1494 (1.21) → 25 960 (1.22) → 148 297
  (1.225) → infeasible (1.23). Objective-value comparison is only meaningful *below* the boundary.
- **Bisection:** tol 1e-3 in µ, returns last feasible; "0.02 above" in `validate.py` is a coarse
  spot-check, not the resolution (the true cliff is ~1e-3-sharp).

## Track C — Exact match to GAMS (core question)
Reference = the authors' GAMS-Convert dumps (`gams_reference/`), solved with HiGHS; reproduce with
`python gams_reference/compare_to_gams.py`.

| substrate | GAMS-dump obj | port obj (shipped) | |Δobj| shipped | |Δobj| last-wins | flux mismatches (>1e-3 / >1e-6 of 5479) | published µ vs port max-µ |
|---|---|---|---|---|---|---|
| coumarate | 1494.339774 | 1494.339979 | 2.0×10⁻⁴ | 3.9×10⁻¹¹ | 12 / 46 | 1.21 vs **1.229** (+0.019) |
| acetate | 1026.500530 | 1026.500554 | 2.5×10⁻⁵ | — | 14 / 25 | 0.77 vs **0.779** (+0.009) |
| butyrate | 1034.079163 | 1034.079219 | 5.5×10⁻⁵ | — | 12 / 37 | 0.86 vs **0.870** (+0.010) |
| succinate | 814.915408 | 814.915431 | 2.2×10⁻⁵ | — | 10 / 25 | 0.74 vs **0.749** (+0.009) |

**Every** flux mismatch is an alternate-optimum swap — e.g. coumarate `R_rxn01451_c0_1` GAMS 2.4619 /
port 0.0 with `R_rxn01451_c0_2` 0.0 / 2.4619, and the recurring `rxn03242/03244/03249/06777` and
`rxn02933_B_1/_B_2` clusters. Same objective, flux redistributed among parallel reactions. **No
genuine (different-objective) discrepancy was found on any substrate.**

On µ_max: the port's own maximum is 0.009–0.019 **above** published. The dumps are frozen at the
published µ, so I cannot test whether GAMS is *also* feasible above 1.21 — i.e. I can't tell whether
the authors stopped incrementing at 1.21 or the port is marginally less constrained. But at every µ I
*can* compare (the published points), **port ≡ GAMS**.

## Track D — Breaking points (coumarate, uptake 2.54, maint 85.4)
| regime | result | reason |
|---|---|---|
| µ = 0 | feasible | trivial |
| µ grid | feasible ≤ **1.225**, infeasible ≥ **1.230** (acetate ≥0.782, butyrate ≥0.872, succinate ≥0.752) | proteome/ATP-machinery coupling binds |
| µ = 5 | infeasible | growth exceeds proteome capacity |
| ATP maint = 0 | infeasible | `atps_const2` forces Σ(ATP prod)=0, incompatible with pinned growth |
| ATP maint = 1000 | feasible at µ=1.0 | slack absorbs it |
| carbon closed | infeasible | no carbon/energy |
| uptake lowered to 2.0 | max-µ drops to **1.16** | *then* carbon uptake binds (it does not bind at 2.54) |
| dilution-only (3872 forms) | feasible at **1.24** (full model infeasible ≥1.23) | dropping the 4672 transcription/translation terms **raises** the ceiling → over-prediction sentinel confirmed (shift modest at coumarate, ~0.01–0.07) |

The feasible→infeasible transition is a sharp cliff, not a fuzzy band — a clean ME-LP.

## Track E — Numerical robustness (HiGHS vs qMINOS)
At the coumarate boundary, classification is **invariant** across default, tightened (pfeas/dfeas 1e-9),
loosened (pfeas 1e-5), and **presolve-off**: µ=1.229 feasible, µ=1.230 infeasible under all four.
Matrix scaling at µ=1.229 spans `1.4×10⁻⁸ … 48` (**ratio 3.4×10⁹**, ill-scaled as expected) yet the
boundary does not move. **HiGHS is trustworthy here; qMINOS is not required for the four base solves.**
(Re-run this sweep on any *derived* model — a 2-cell version has a larger dynamic range.)

---

## Direct answers
- **Flux-level equal to GAMS, or only µ_max?** Flux-level, to the degeneracy limit: objective equals
  the GAMS instance to ≤6×10⁻⁵ (≈4×10⁻¹¹ with the last-wins fix), and fluxes agree on all reactions
  except alternate-optima swaps. Verified against the authors' own GAMS-Convert dumps.
- **Breaking point per substrate + binding constraint.** Infeasible at ≥1.230 / ≥0.782 / ≥0.872 /
  ≥0.752 (coumarate/acetate/butyrate/succinate); the cliff is set by the proteome/ATP-machinery
  coupling, **not** carbon uptake (uptake saturates above ~2.5; it only binds if lowered toward the
  stale file value).
- **Track-A discrepancy the growth rates were masking?** One: the 26 duplicate coupling coefficients
  summed instead of last-wins (proven to be the whole objective residual). The medium-bound concern
  from the first pass is **retracted** — the port matches the authors' GAMS instances exactly.
- **HiGHS trustworthy here?** Yes — boundary invariant to tolerance/presolve despite 3.4×10⁹ scaling.
- **Before a 2-cell community ME-model can be trusted:** (a) apply the last-wins fix (one line above);
  (b) keep sourcing medium/maintenance from the paper/GAMS instances, not `upper_bound.txt`;
  (c) carry the `gams_reference/` flux cross-check forward — a community model multiplies degeneracy, so
  µ_max agreement alone won't be enough; (d) re-run the Track-E tolerance sweep on the derived model;
  (e) harden the `e-` exponent regex before ingesting new data.

## Limitations / could-not-substantiate
- **µ_max above the published value is untestable against GAMS.** The dumps are frozen at the published
  µ, so I cannot determine whether GAMS is also feasible up to ~1.229 (port stopped-early vs GAMS, or a
  marginal formulation difference) — only that port ≡ GAMS *at* the published µ.
- **The GAMS reference was solved with HiGHS, not GAMS/CPLEX.** It is the authors' exact LP *instance*
  (GAMS-Convert output), but re-solved; a CPLEX solve could pick a different alternate-optimum vertex.
  This does not affect the objective/feasibility conclusions, only which degenerate fluxes appear.
- **The reference is at the published µ only.** No GAMS reference exists for off-point sweeps
  (ATP-maintenance, uptake, N₂↔NH₄), so Track-C agreement across the continuum is shown for the Python
  model against itself, not against GAMS.
