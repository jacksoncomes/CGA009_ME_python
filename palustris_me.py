"""
palustris_me.py
===============
The *Rhodopseudomonas palustris* CGA009 genome-scale Metabolism & Expression (ME) model,
running independently in pure Python -- no GAMS, no cobra, no license.

It reconstructs the exact ME linear program of Chowdhury, Alsiyabi & Saha (Microbiol.
Spectrum 2022; GitHub ssbio/palustris_ME_model) from the original GAMS data files, and
solves it at any fixed growth rate mu with HiGHS (scipy.optimize.linprog). A bisection
returns the maximum feasible mu -- the ME-predicted growth rate.

The ME LP (exactly as in the GAMS source `palustris_ME.gms`):

    minimize   sum_j v_j                       # pFBA objective (primalobj)
    s.t.       S(mu) . v = 0                    # mass balance with growth-rate coupling
               v[bio2] = mu                     # growth pinned to mu
               sum_{j: ATP produced} v_j = M    # ATP maintenance (atps_const2), substrate-specific
               sum_{j: ATP produced} v_j = v10042 + v01517 + v00148B1 + v00148B2   # atps_const
               v[R_rxn10042_c0] = v[R_rxn37614_c0]   # ATP synthase = photosystem (atps_eq_ps2)
               v[R_rxn00288_c0] = 0             # succTOfum
               lb <= v <= ub, + fixed/stopped reactions

`S(mu)` differs from an ordinary stoichiometric matrix by ~8,500 growth-rate-dependent
"coupling" coefficients that encode enzyme dilution and the transcription/translation/
ribosome machinery cost. These are read symbolically from the .gms and re-evaluated at
each mu, which is what makes the whole thing a plain LP at fixed mu (and hence solvable
without qMINOS/GAMS).

Validated to reproduce the Saha-lab Table-2 growth rates exactly (see reproduce_saha_results.ipynb).

Author of this Python port: (your name).  Original model: Saha lab, ssbio/palustris_ME_model.
"""
from __future__ import annotations
import os
import re
import numpy as np
from scipy.sparse import lil_matrix, coo_matrix, vstack
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
GMS = os.path.join(DATA, "palustris_ME.gms")

# E. coli-borrowed scalars used by the coupling equations (palustris_ME.gms lines 54-58)
SCALARS = {"kcat": 234000.0, "kt": 108.0, "r0": 4.5, "cmrna2": 16072.50, "cribo2": 1976.8}
ATP_C = "cpd00002[c0]"

# Saha-lab Table-2 conditions: name -> (uptake exchange, uptake bound, ATP maintenance, published mu)
SUBSTRATES = {
    "coumarate": ("EX_cpd00604_e0_B", 2.54, 85.4, 1.21),
    "acetate":   ("EX_cpd00029_e0_B", 6.47, 54.0, 0.77),
    "butyrate":  ("EX_cpd00211_e0_B", 3.69, 56.7, 0.86),
    "succinate": ("EX_cpd00036_e0_B", 4.66, 45.7, 0.74),
}
_CARBON_EX = [v[0] for v in SUBSTRATES.values()]


# --------------------------------------------------------------------------- parsing
def _quoted_list(path):
    out = []
    for ln in open(path):
        s = ln.strip()
        if not s or s == "/" or s.startswith("*"):
            continue
        m = re.findall(r"'([^']+)'", s)
        if m:
            out.append(m[0])
    return out


def _bounds(path):
    d = {}
    for ln in open(path):
        m = re.match(r"'([^']+)'\s+(-?\.?[\d.eE+]+)", ln.strip())
        if m:
            d[m.group(1)] = float(m.group(2))
    return d


def _prep_rhs(rhs):
    """Turn a GAMS coupling RHS into a python expression string of `mu` only, substituting
    the E. coli scalars. Handles every form found in the .gms:
        -mu/DENOM,  -mu/kcat,  -mu/(a+b+...),
        -N*(mu+(kt*r0))/(kt*cribo2/X),  -(mu+(kt*r0))/(kt*cmrna2/X)."""
    rhs = re.sub(r"\s+", " ", rhs).strip().rstrip(";").strip()
    for k, v in sorted(SCALARS.items(), key=lambda kv: -len(kv[0])):   # longest name first
        rhs = re.sub(r"\b%s\b" % k, repr(v), rhs)
    if not re.fullmatch(r"[0-9.eE+\-*/() mu]+", rhs):
        raise ValueError("unparseable coupling RHS: " + rhs[:60])
    return rhs


def parse_model():
    rxns = _quoted_list(os.path.join(DATA, "reactions.txt"))
    mets = _quoted_list(os.path.join(DATA, "metabolites.txt"))
    lb = _bounds(os.path.join(DATA, "lower_bound.txt"))
    ub = _bounds(os.path.join(DATA, "upper_bound.txt"))
    # base stoichiometry from sij.txt: 'met'.'rxn' coef
    S = {}
    for ln in open(os.path.join(DATA, "sij.txt")):
        m = re.match(r"'([^']+)'\.'([^']+)'\s+(-?\.?[\d.eE+]+)", ln.strip())
        if m:
            S[(m.group(1), m.group(2))] = float(m.group(3))
    # ALL growth-rate-dependent S overwrites from the .gms (statements may span lines)
    text = open(GMS).read()
    coupling = []
    for met, rxn, rhs in re.findall(r"S\('([^']+)','([^']+)'\)\s*=\s*(.+?);", text, re.DOTALL):
        coupling.append((met, rxn, compile(_prep_rhs(rhs), "<rhs>", "eval")))
    stopped = set(re.findall(r"v\.(?:up|lo)\('([^']+)'\)\s*=\s*0\s*;", text))
    return dict(rxns=rxns, mets=mets, lb=lb, ub=ub, S=S, coupling=coupling, stopped=stopped)


# --------------------------------------------------------------------------- the model
class MEModel:
    """Pure-Python R. palustris ME-model. Load once, then solve() at any mu."""

    def __init__(self):
        M = parse_model()
        self.rxns, self.mets = M["rxns"], M["mets"]
        self.rj = {r: j for j, r in enumerate(self.rxns)}
        self.mi = {m: i for i, m in enumerate(self.mets)}
        self.n, self.m = len(self.rxns), len(self.mets)
        # bounds (GAMS default lb 0 / ub Vmax=1000; files override; stopped reactions -> 0)
        self.vlb = np.zeros(self.n)
        self.vub = np.full(self.n, 1000.0)
        for r, v in M["lb"].items():
            if r in self.rj:
                self.vlb[self.rj[r]] = v
        for r, v in M["ub"].items():
            if r in self.rj:
                self.vub[self.rj[r]] = v
        for r in M["stopped"]:
            if r in self.rj:
                self.vlb[self.rj[r]] = self.vub[self.rj[r]] = 0.0
        # split S into static (mu-independent) and coupling (mu-dependent, overwrites placeholders)
        cpairs = {(p, rx) for p, rx, _ in M["coupling"]}
        srow, scol, sval = [], [], []
        for (met, rx), c in M["S"].items():
            if (met, rx) in cpairs or met not in self.mi or rx not in self.rj:
                continue
            srow.append(self.mi[met]); scol.append(self.rj[rx]); sval.append(c)
        self._sval = np.array(sval, float)
        cr, cc, code = [], [], []
        for p, rx, cd in M["coupling"]:
            if p in self.mi and rx in self.rj:
                cr.append(self.mi[p]); cc.append(self.rj[rx]); code.append(cd)
        self._ccode = code
        self.n_coupling = len(code)
        self._rows = np.array(srow + cr)
        self._cols = np.array(scol + cc)
        # ATP-producing reactions (S[ATP,j] > 0) for the maintenance constraints
        ai = self.mi.get(ATP_C)
        self.atp_prod = [j for i, j, c in zip(srow, scol, sval) if i == ai and c > 0]

    # ------------------------------------------------------------------ medium
    def set_medium(self, substrate):
        """Configure a Table-2 substrate. Returns (uptake, atp_maint, published_mu)."""
        ex, up, atp, mu = SUBSTRATES[substrate]
        for c in _CARBON_EX:
            if c in self.rj:
                self.vub[self.rj[c]] = 0.0
        self.vub[self.rj[ex]] = up
        return up, atp, mu

    # ------------------------------------------------------------------ solve
    def _S_at(self, mu):
        env = {"mu": mu, "__builtins__": {}}
        cval = np.fromiter((eval(cd, env) for cd in self._ccode), float, self.n_coupling)
        data = np.concatenate([self._sval, cval])
        return coo_matrix((data, (self._rows, self._cols)), shape=(self.m, self.n)).tocsr()

    def solve(self, mu, atp_maint=85.4):
        """Solve the fixed-mu ME LP with HiGHS. Returns (status, objective, fluxes | None).
        `objective` is min sum(v); `atp_maint` is the substrate-specific ATP maintenance
        (coumarate 85.4, acetate 54, butyrate 56.7, succinate 45.7)."""
        A_eq = self._S_at(mu)
        extra = lil_matrix((4, self.n)); be = np.zeros(4)
        for j in self.atp_prod:            # atps_const2: sum ATP production = atp_maint
            extra[0, j] = 1.0
        be[0] = atp_maint
        for j in self.atp_prod:            # atps_const: sum ATP prod - (4 synth rxns) = 0
            extra[1, j] += 1.0
        for r in ("R_rxn10042_c0", "R_rxn01517_c0", "R_rxn00148_c0_B_1", "R_rxn00148_c0_B_2"):
            if r in self.rj:
                extra[1, self.rj[r]] -= 1.0
        extra[2, self.rj["R_rxn10042_c0"]] = 1.0     # atps_eq_ps2
        extra[2, self.rj["R_rxn37614_c0"]] = -1.0
        extra[3, self.rj["R_rxn00288_c0"]] = 1.0     # succTOfum
        A = vstack([A_eq, extra.tocsr()]).tocsr()
        b = np.concatenate([np.zeros(self.m), be])
        lb, ub = self.vlb.copy(), self.vub.copy()
        j2 = self.rj["bio2"]; lb[j2] = ub[j2] = mu    # pin growth to mu
        res = linprog(np.ones(self.n), A_eq=A, b_eq=b, bounds=list(zip(lb, ub)), method="highs")
        if res.success:
            return "optimal", float(res.fun), res.x
        return ("infeasible" if "infeasible" in res.message.lower() else res.message), None, None

    def flux(self, solution, reaction):
        """Convenience: flux of a reaction id from a solve() flux vector."""
        return float(solution[self.rj[reaction]]) if reaction in self.rj else float("nan")

    def bisect_max_mu(self, atp_maint=85.4, lo=0.0, hi=2.5, tol=1e-3):
        """Largest mu for which the ME LP is feasible = ME max growth rate."""
        assert self.solve(lo, atp_maint)[0] == "optimal", "lower bound infeasible"
        best = lo
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            if self.solve(mid, atp_maint)[0] == "optimal":
                best = lo = mid
            else:
                hi = mid
        return best


if __name__ == "__main__":
    import time
    t0 = time.time()
    me = MEModel()
    print("ME-model loaded: %d reactions, %d metabolites, %d mu-coupling entries (%.1fs)\n"
          % (me.n, me.m, me.n_coupling, time.time() - t0))
    print("%-11s  %-9s  %-9s  %s" % ("substrate", "Saha mu", "my max mu", "match"))
    print("-" * 44)
    for name in SUBSTRATES:
        up, atp, saha = me.set_medium(name)
        mx = me.bisect_max_mu(atp_maint=atp)
        print("%-11s  %-9.2f  %-9.3f  %s" % (name, saha, mx, "OK" if abs(mx - saha) <= 0.02 else "XX"))
    print("\ntotal %.1fs" % (time.time() - t0))
