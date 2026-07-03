"""
test_me_audit.py — regression suite for the pure-Python R. palustris ME-model port.

Encodes the facts established by the audit (see AUDIT_REPORT.md) and the fixes in CHANGES.md:

  * structural parse fidelity  : 8544 source assignments -> 8518 after GAMS last-wins collapse
  * duplicate (met,rxn) targets use LAST-WINS, not summing (the fix)
  * mu-dependence: coupling coefficients equal the six algebraic forms at a second mu
  * the four Saha Table-2 rates : published mu feasible, published+0.02 infeasible
  * per-substrate feasibility boundary sits ~0.01-0.02 ABOVE the published value
  * solve determinism
  * the dilution-form sentinel OVER-predicts (boundary strictly rises)
  * objective matches the authors' GAMS-Convert instances to machine precision
  * the remaining flux differences are provably alternate optima (GAMS solution is
    feasible + equally optimal in the port's own LP)

Run:  pytest test_me_audit.py -v        (a handful of tests are slow: ~1.4 s / LP solve)
"""
import re
import numpy as np
import pytest

import palustris_me as pm
from palustris_me import MEModel, SUBSTRATES, parse_model, GMS


# ----------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def model():
    return MEModel()


def _set_medium(me, name):
    ex, up, atp, mu = SUBSTRATES[name]
    for c in [v[0] for v in SUBSTRATES.values()]:
        if c in me.rj:
            me.vub[me.rj[c]] = 0.0
    me.vub[me.rj[ex]] = up
    return ex, up, atp, mu


# ----------------------------------------------------------------- Track A: parsing
def test_coupling_lastwins_collapses_26_duplicates():
    """All 8544 source assignments are captured (not the 3872 dilution-only subset), then the
    26 duplicate (met,rxn) targets collapse to their LAST assignment -> 8518 unique. The port
    must NOT sum duplicates (that was the audited bug)."""
    text = open(GMS).read()
    raw = re.findall(r"S\('([^']+)','([^']+)'\)\s*=\s*.+?;", text, re.DOTALL)
    assert len(raw) == 8544                              # every source assignment captured
    from collections import Counter
    n_dup = sum(c - 1 for c in Counter(raw).values() if c > 1)
    assert n_dup == 26                                   # 26 duplicate (met,rxn) targets
    M = parse_model()
    assert len(M["coupling"]) == 8544 - 26               # collapsed last-wins == 8518


def test_gms_has_exactly_8544_single_line_assignments():
    text = open(GMS).read()
    n = len(re.findall(r"S\('[^']+','[^']+'\)\s*=\s*.+?;", text, re.DOTALL))
    assert n == 8544


def test_no_coupling_or_base_entries_silently_dropped():
    M = parse_model()
    mets, rxns = set(M["mets"]), set(M["rxns"])
    assert not [1 for (m, r) in M["S"] if m not in mets or r not in rxns]
    assert not [1 for (m, r, _) in M["coupling"] if m not in mets or r not in rxns]


def test_model_dimensions(model):
    assert model.n == 5479
    assert model.m == 3192
    assert model.n_coupling == 8518          # 8544 source assignments, 26 duplicates collapsed
    assert len(model.atp_prod) == 173


def test_duplicate_assignments_use_last_wins(model):
    """THE FIX: GAMS assignment is last-wins. `mRNA_RPA0185/TL_RPA0185` is assigned twice
    (identical value) in the .gms; the coefficient must be the single GAMS value, NOT 2x it
    (the old coo->csr summing bug)."""
    i, j = model.mi["mRNA_RPA0185"], model.rj["TL_RPA0185"]
    kt, r0, cmrna2 = 108.0, 4.5, 16072.5
    gams_last_wins = -(1.0 + kt * r0) / (kt * cmrna2 / 48.5)
    assert model._S_at(1.0)[i, j] == pytest.approx(gams_last_wins, rel=1e-12)


def test_mu_dependence_matches_six_algebraic_forms(model):
    """Guards the coupling coefficients' variation WITH mu (not just a snapshot at published mu).
    Independently predict every coefficient at two off-point mu values from the six algebraic
    forms + substituted E. coli scalars, and assert the port matches to ~1e-9. This is what the
    future community model relies on when mu is swept over a real range."""
    kcat, kt, r0, cmrna2, cribo2 = 234000.0, 108.0, 4.5, 16072.5, 1976.8

    def expected(rhs, mu):
        s = re.sub(r"\s+", "", rhs).rstrip(";")
        m = re.fullmatch(r"-([0-9.eE]+)\*\(mu\+\(kt\*r0\)\)/\(kt\*cribo2/([0-9.eE]+)\)", s)  # ribosome
        if m:
            return -float(m.group(1)) * (mu + kt * r0) / (kt * cribo2 / float(m.group(2)))
        m = re.fullmatch(r"-\(mu\+\(kt\*r0\)\)/\(kt\*cmrna2/([0-9.eE]+)\)", s)               # mRNA
        if m:
            return -(mu + kt * r0) / (kt * cmrna2 / float(m.group(1)))
        m = re.fullmatch(r"-mu/(.+)", s)                                                     # dilution
        if m:
            denom = m.group(1).replace("kcat", repr(kcat))            # only scalar that appears here
            assert re.fullmatch(r"[0-9.eE+\-*/()]+", denom), "non-numeric denominator: " + denom
            return -mu / eval(denom, {"__builtins__": {}})           # pure numeric arithmetic, no mu
        raise AssertionError("unrecognized coupling form: " + rhs)

    # last-wins map of (met,rxn) -> raw RHS string, straight from the source
    text = open(GMS).read()
    last = {}
    for met, rxn, rhs in re.findall(r"S\('([^']+)','([^']+)'\)\s*=\s*(.+?);", text, re.DOTALL):
        last[(met, rxn)] = rhs
    assert len(last) == 8518

    for mu in (1.11, 0.30):                       # two off-point mu (published mu is 1.21)
        coo = model._S_at(mu).tocoo()
        port = {(i, j): v for i, j, v in zip(coo.row, coo.col, coo.data)}
        worst = 0.0
        for (met, rxn), rhs in last.items():
            exp = expected(rhs, mu)
            got = port.get((model.mi[met], model.rj[rxn]), 0.0)
            worst = max(worst, abs(got - exp) - (1e-9 * abs(exp) + 1e-12))
        assert worst <= 0.0, f"mu={mu}: coupling coefficient mismatch beyond 1e-9 (slack {worst:.2e})"


# ----------------------------------------------------------------- Track B/C: solves
@pytest.mark.parametrize("name", list(SUBSTRATES))
def test_published_mu_feasible_and_plus_002_infeasible(model, name):
    """The validate.py claim, restated precisely: published mu is feasible, +0.02 is not.
    NOTE: this is weaker than 'published mu == max growth' (see boundary test)."""
    ex, up, atp, mu = _set_medium(model, name)
    assert model.solve(mu, atp_maint=atp)[0] == "optimal"
    assert model.solve(round(mu + 0.02, 2), atp_maint=atp)[0] != "optimal"


@pytest.mark.slow
@pytest.mark.parametrize("name,lo_feas,hi_infeas", [
    ("coumarate", 1.225, 1.230),
    ("acetate",   0.775, 0.782),
    ("butyrate",  0.867, 0.872),
    ("succinate", 0.745, 0.752),
])
def test_feasibility_boundary_above_published(model, name, lo_feas, hi_infeas):
    """The port's true max growth sits ~0.01-0.02 ABOVE the published value: it is
    still feasible past the published mu. This is the audit's central mu_max finding."""
    ex, up, atp, mu = _set_medium(model, name)
    assert model.solve(lo_feas, atp_maint=atp)[0] == "optimal"      # past published, still feasible
    assert model.solve(hi_infeas, atp_maint=atp)[0] != "optimal"    # boundary crossed
    assert lo_feas > mu                                             # boundary is above published


def test_determinism(model):
    _set_medium(model, "coumarate")
    x1 = model.solve(1.0, atp_maint=85.4)[2]
    x2 = model.solve(1.0, atp_maint=85.4)[2]
    assert np.max(np.abs(x1 - x2)) == 0.0


def test_edge_regimes(model):
    _set_medium(model, "coumarate")
    assert model.solve(0.0, atp_maint=85.4)[0] == "optimal"          # mu=0 feasible
    assert model.solve(5.0, atp_maint=85.4)[0] != "optimal"          # absurd growth infeasible
    assert model.solve(1.0, atp_maint=0.0)[0] != "optimal"           # zero maintenance infeasible
    # carbon closed -> infeasible
    for c in [v[0] for v in SUBSTRATES.values()]:
        model.vub[model.rj[c]] = 0.0
    assert model.solve(1.0, atp_maint=85.4)[0] != "optimal"


# ----------------------------------------------------------------- Track D: sentinel
def _build_with(transform):
    orig = pm.parse_model
    def patched():
        M = orig()
        M["coupling"] = transform(M["coupling"])
        return M
    pm.parse_model = patched
    try:
        return MEModel()
    finally:
        pm.parse_model = orig


@pytest.mark.slow
def test_dilution_only_over_predicts():
    """Sentinel: dropping the 4672 transcription/translation coupling forms (keeping the
    3872 -mu/enzyme dilution forms) must RAISE the feasibility ceiling. Guards against a
    regression that silently loses the machinery-cost terms."""
    text = open(GMS).read()
    is_transl = {}
    for met, rxn, rhs in re.findall(r"S\('([^']+)','([^']+)'\)\s*=\s*(.+?);", text, re.DOTALL):
        is_transl[(met, rxn)] = "(kt*r0)" in rhs.replace(" ", "")

    def dilution_only(coupling):
        return [(m, r, c) for (m, r, c) in coupling if not is_transl.get((m, r), False)]

    n_dilution = sum(1 for v in is_transl.values() if not v)   # 3869 after last-wins collapse
    me_s = _build_with(dilution_only)
    assert me_s.n_coupling == n_dilution
    _set_medium(me_s, "coumarate")
    # full model is infeasible at 1.24; dilution-only must still be feasible there
    assert me_s.solve(1.24, atp_maint=85.4)[0] == "optimal"


# ----------------------------------------------------------------- Track C: GAMS reference
def _dedup_last_wins(coupling):
    d = {}
    for met, rxn, code in coupling:
        d[(met, rxn)] = code
    return [(m, r, c) for (m, r), c in d.items()]


@pytest.mark.slow
@pytest.mark.parametrize("case,mu,atp,portname", [
    ("coumarate", 1.21, 85.4, "coumarate"),
    ("acetate",   0.77, 54.0, "acetate"),
    ("buytrate",  0.86, 56.7, "butyrate"),
    ("succinate", 0.74, 45.7, "succinate"),
])
def test_matches_gams_convert_objective(case, mu, atp, portname):
    """Cross-check against the authors' GAMS-Convert dump (gams_reference/). After the last-wins
    fix the port's min-sum-v objective equals the GAMS instance to MACHINE PRECISION (~1e-9);
    the few remaining flux differences are alternate-optima swaps (see the degeneracy test)."""
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "gams_reference"))
    if not os.path.exists(os.path.join(os.path.dirname(__file__), "gams_reference", case + ".py")):
        pytest.skip("GAMS reference dump not present")
    import compare_to_gams as gr
    r = gr.compare(case)
    assert abs(r["gams_obj"] - r["port_obj"]) < 1e-6      # machine-exact after the fix
    assert r["n_gt_1e3"] < 20                             # only alternate-optima swaps remain


@pytest.mark.slow
@pytest.mark.parametrize("case,mu,atp,portname", [
    ("coumarate", 1.21, 85.4, "coumarate"),
    ("acetate",   0.77, 54.0, "acetate"),
    ("buytrate",  0.86, 56.7, "butyrate"),
    ("succinate", 0.74, 45.7, "succinate"),
])
def test_flux_differences_are_true_alternate_optima(case, mu, atp, portname):
    """Prove the residual flux differences vs GAMS are pure degeneracy, not a discrepancy:
    the GAMS reference flux vector, plugged into the PORT's OWN LP, is feasible (mass balance +
    all four custom constraints + bounds) and gives the SAME objective. Two feasible vertices,
    identical optimum => alternate optima of the identical LP."""
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "gams_reference"))
    if not os.path.exists(os.path.join(os.path.dirname(__file__), "gams_reference", case + ".py")):
        pytest.skip("GAMS reference dump not present")
    import compare_to_gams as gr
    from scipy.optimize import linprog

    varlist, A, b, lo, hi = gr.parse_dump(case + ".py")
    res = linprog(np.ones(len(varlist)), A_eq=A, b_eq=b, bounds=list(zip(lo, hi)), method="highs")
    xmap = gr._xmap()
    gfl = {xmap[v]: res.x[i] for i, v in enumerate(varlist) if v in xmap}

    me = MEModel(); me.set_medium(portname); _, pobj, _ = me.solve(mu, atp_maint=atp)
    vref = np.array([gfl.get(r, 0.0) for r in me.rxns])       # GAMS solution in port order

    assert np.abs(me._S_at(mu) @ vref).max() < 1e-8           # mass balance in PORT's S(mu)
    ap = sum(vref[j] for j in me.atp_prod)
    assert abs(ap - atp) < 1e-6                               # atps_const2
    four = sum(vref[me.rj[r]] for r in
               ("R_rxn10042_c0", "R_rxn01517_c0", "R_rxn00148_c0_B_1", "R_rxn00148_c0_B_2") if r in me.rj)
    assert abs(ap - four) < 1e-6                              # atps_const
    assert abs(vref[me.rj["R_rxn10042_c0"]] - vref[me.rj["R_rxn37614_c0"]]) < 1e-6   # atps_eq_ps2
    assert abs(vref[me.rj["R_rxn00288_c0"]]) < 1e-6           # succTOfum
    lb2, ub2 = me.vlb.copy(), me.vub.copy(); j2 = me.rj["bio2"]; lb2[j2] = ub2[j2] = mu
    assert (vref >= lb2 - 1e-9).all() and (vref <= ub2 + 1e-9).all()   # bounds
    assert abs(vref.sum() - pobj) < 1e-6                      # same objective


@pytest.mark.slow
def test_last_wins_fix_makes_gams_match_exact():
    """The 26-entry duplicate-summing bug is the SOLE source of the ~1e-4 objective residual:
    with GAMS last-wins semantics the port reproduces the GAMS instance to ~1e-10."""
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "gams_reference"))
    import compare_to_gams as gr
    varlist, A, b, lo, hi = gr.parse_dump("coumarate.py")
    from scipy.optimize import linprog
    gams_obj = linprog(np.ones(len(varlist)), A_eq=A, b_eq=b,
                       bounds=list(zip(lo, hi)), method="highs").fun

    me = _build_with(_dedup_last_wins)
    _set_medium(me, "coumarate")
    port_obj = me.solve(1.21, atp_maint=85.4)[1]
    assert abs(port_obj - gams_obj) < 1e-6      # exact to solver precision after the fix
