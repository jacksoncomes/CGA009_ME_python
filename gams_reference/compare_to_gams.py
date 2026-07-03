"""
compare_to_gams.py — flux-level cross-check of the Python port against the Saha lab's own
GAMS-Convert dumps (the authors' GAMS model instance, frozen at the published growth rate).

The four `*_max_growth.py` files in the upstream repo (ssbio/palustris_ME_model, PYTHON/Table_2/)
were produced by GAMS Convert: they are the exact GAMS LP for each substrate, serialized to Pyomo,
with the growth-rate coupling coefficients baked in at the published mu (bio2 fixed to 1.21 / 0.77 /
0.86 / 0.74). `dictionary.txt` maps the anonymous xN variables back to reaction ids.

We parse each dump numerically and solve it with the SAME HiGHS backend the port uses, so any
objective/flux difference is attributable to the *formulation*, not the solver. Result (see
AUDIT_REPORT.md): objective matches to <=6e-5 as shipped, to ~4e-11 once the duplicate-summing bug
is fixed; fluxes match on all 5479 reactions except ~25-46 alternate-optima swaps between parallel
isozyme/direction variants.

Usage:  python compare_to_gams.py [coumarate acetate buytrate succinate]
"""
import os, re, sys
import numpy as np
from scipy.sparse import coo_matrix
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import palustris_me as pm

# published mu / ATP maintenance per substrate, and the port's substrate key
CASES = {
    "coumarate": ("coumarate.py", 1.21, 85.4, "coumarate"),
    "acetate":   ("acetate.py",   0.77, 54.0, "acetate"),
    "buytrate":  ("buytrate.py",  0.86, 56.7, "butyrate"),   # upstream spelling
    "succinate": ("succinate.py", 0.74, 45.7, "succinate"),
}


def _xmap():
    d = {}
    for ln in open(os.path.join(HERE, "dictionary.txt")):
        m = re.match(r"\s*(x\d+)\s+v\('([^']+)'\)", ln)
        if m:
            d[m.group(1)] = m.group(2)
    return d


def parse_dump(pyfile):
    """Return (varlist, A_eq, b_eq, lb, ub) for the GAMS-Convert LP (objective is min sum v)."""
    txt = open(os.path.join(HERE, pyfile)).read()
    lb, ub = {}, {}
    for m in re.finditer(r"m\.(x\d+)\s*=\s*Var\(within=Reals,bounds=\(([^,]+),([^)]+)\)", txt):
        lb[m.group(1)] = float(m.group(2)); ub[m.group(1)] = float(m.group(3))
    varlist = sorted(lb, key=lambda s: int(s[1:]))
    vidx = {v: i for i, v in enumerate(varlist)}
    rows, cols, data, rhs = [], [], [], []
    r = 0
    for m in re.finditer(r"Constraint\(expr=(.*?)\)\s*\n", txt, re.DOTALL):
        expr = m.group(1).replace("\n", " ")
        lhs, rval = expr.rsplit("==", 1)
        rhs.append(float(rval.strip()))
        for term in (" " + lhs.strip()).replace(" - ", " + -").split(" + "):
            term = term.strip()
            if not term:
                continue
            if "*" in term:
                c, v = term.split("*"); coef = float(c.replace(" ", ""))
            elif term.startswith("-"):
                coef, v = -1.0, term[1:]
            else:
                coef, v = 1.0, term
            rows.append(r); cols.append(vidx[v.replace("m.", "").strip()]); data.append(coef)
        r += 1
    A = coo_matrix((data, (rows, cols)), shape=(r, len(varlist))).tocsr()
    lo = np.array([lb[v] for v in varlist]); hi = np.array([ub[v] for v in varlist])
    return varlist, A, np.array(rhs), lo, hi


def compare(case):
    pyf, mu, atp, portname = CASES[case]
    varlist, A, b, lo, hi = parse_dump(pyf)
    res = linprog(np.ones(len(varlist)), A_eq=A, b_eq=b, bounds=list(zip(lo, hi)), method="highs")
    xmap = _xmap()
    gfl = {xmap[v]: res.x[i] for i, v in enumerate(varlist) if v in xmap}

    me = pm.MEModel(); me.set_medium(portname)
    st, pobj, x = me.solve(mu, atp_maint=atp)
    pfl = {rr: x[j] for rr, j in me.rj.items()}

    common = [rr for rr in gfl if rr in pfl]
    dev = np.array([abs(gfl[rr] - pfl[rr]) for rr in common])
    return dict(gams_obj=res.fun, port_obj=pobj, obj_absdiff=abs(res.fun - pobj),
                n_common=len(common), max_dflux=float(dev.max()),
                n_gt_1e6=int((dev > 1e-6).sum()), n_gt_1e3=int((dev > 1e-3).sum()))


if __name__ == "__main__":
    for case in (sys.argv[1:] or list(CASES)):
        r = compare(case)
        print(f"{case:10s} obj GAMS={r['gams_obj']:.6f} port={r['port_obj']:.6f} "
              f"|d|={r['obj_absdiff']:.2e}  flux: {r['n_common']} shared, "
              f"max|d|={r['max_dflux']:.2e}, #>1e-6={r['n_gt_1e6']}, #>1e-3={r['n_gt_1e3']}")
