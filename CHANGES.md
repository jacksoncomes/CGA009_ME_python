# CHANGES — `fix/last-wins-coupling`

Applies the fixes identified in `AUDIT_REPORT.md`. No biology changed; these bring the port into
exact agreement with the authors' GAMS model instances and make the medium provenance explicit.

## 1. Last-wins coupling coefficients (the real bug)
`palustris_me.py :: parse_model()` previously appended every `S(met,rxn)=…` assignment to a list;
duplicate `(met,rxn)` coordinates were then **summed** by `coo_matrix(...).tocsr()`, doubling the
26 coefficients that the `.gms` assigns twice. GAMS assignment semantics are **last-wins**.

Fix: collapse duplicates in `parse_model()` keeping the **last** occurrence (a dict keyed on
`(met,rxn)`), so each matrix cell gets the final assignment, never a sum.

- Coupling entries: **8544 → 8518** (26 duplicates collapsed; all 8544 source assignments are still
  read — this is a collapse, not a drop).
- `MEModel.n_coupling`: 8544 → 8518.

**Effect — objective gap vs the authors' GAMS-Convert instance (min Σv at published µ):**

| substrate | before fix | after fix |
|---|---|---|
| coumarate | 2.04×10⁻⁴ | **3.9×10⁻¹¹** |
| acetate | 2.46×10⁻⁵ | **1.8×10⁻¹⁰** |
| butyrate | 5.54×10⁻⁵ | **3.2×10⁻¹⁰** |
| succinate | 2.22×10⁻⁵ | **7.3×10⁻¹²** |

The bug was the *entire* source of the residual: the port now reproduces each GAMS instance to
machine precision. Reproduce with `python gams_reference/compare_to_gams.py`.

## 2. New test: µ-dependence of the coupling coefficients
`test_me_audit.py :: test_mu_dependence_matches_six_algebraic_forms`. Every prior check verified the
coefficients only at the published µ — a single snapshot. This test independently predicts every
coupling coefficient at two off-point µ values (1.11 and 0.30) from the six algebraic forms + the
substituted E. coli scalars (kcat 234000, kt 108, r0 4.5, cmrna2 16072.5, cribo2 1976.8) and asserts
the port matches to ~1e-9. It is the guard the future community model relies on, where µ is swept over
a real range rather than evaluated at one point. (No GAMS needed — pure internal consistency.)

Also added `test_flux_differences_are_true_alternate_optima`: plugs the GAMS reference flux vector into
the port's **own** LP and shows it is feasible (mass balance + all four custom constraints + bounds)
and equally optimal — proving the residual per-reaction flux differences are alternate optima of the
identical LP, not a discrepancy.

## 3. Quarantined the stale medium file
`data/upper_bound.txt` lists coumarate uptake **2.0** and the other three carbons **0** — a stale
snapshot. The authoritative uptakes (**2.54 / 6.47 / 3.69 / 4.66**) and ATP maintenance
(**85.4 / 54 / 56.7 / 45.7**) are the authors' GAMS-Convert instances (`gams_reference/`), which the
port already uses.

- The file is **not renamed**: it is authoritative for *all non-carbon* upper bounds and is loaded in
  full by `parse_model()`. Renaming it would break the loader and mislabel ~5400 correct bounds. Instead
  a loud `*`-comment header (ignored by the `_bounds()` parser) marks the carbon-EX rows as stale and
  points to the GAMS instances as ground truth.
- **Why the port yields 2.54 despite the file's 2.0:** it is a *runtime override*, not a different
  source file. `MEModel.set_medium()` sets `vub[carbon EX]` from the `SUBSTRATES` dict (now annotated
  with its provenance) *after* the file is loaded. The 2.0 row is never used for coumarate. (Uptake 2.0
  would give µ_max 1.16, not the published 1.21.)

## Verification (all on this branch)
- **`test_me_audit.py`: 26 passed** (21 original checks, several split/tightened, plus the two new
  tests above). Tolerances tightened to the now-exact values, not loosened: the GAMS-objective check is
  `< 1e-6` (was `< 1e-3`); the duplicate-coefficient check asserts the single GAMS value, not 2×.
- **µ_max offsets unchanged** — the fix perturbs 26 coefficients too little to move the boundary
  (bisection, tol 1e-4):

  | substrate | published µ | post-fix max-µ | offset | pre-fix max-µ |
  |---|---|---|---|---|
  | coumarate | 1.21 | 1.2291 | +0.0191 | 1.229 |
  | acetate | 0.77 | 0.7788 | +0.0088 | 0.779 |
  | butyrate | 0.86 | 0.8699 | +0.0099 | 0.870 |
  | succinate | 0.74 | 0.7487 | +0.0087 | 0.749 |

- `validate.py` still exits 0.

## Carry-forward for the 2-cell community model (not started here)
- Keep the `gams_reference/` cross-check in the loop: a community model multiplies the alternate-optima
  degeneracy, so µ_max agreement alone won't be sufficient evidence — assert objective + GAMS-solution
  feasibility as done here.
- Source medium/maintenance from the paper/GAMS instances, never from `upper_bound.txt`.
- The µ-dependence test matters more once µ is swept over a range; extend it to any new coupling forms.
- Latent (not triggered today): the `_bounds`/`sij`/`_prep_rhs` regexes don't accept a `-` in a
  scientific-notation exponent; harden before ingesting new data.
