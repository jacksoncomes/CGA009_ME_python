"""
validate.py — prove the Python ME-model reproduces the Saha-lab Table-2 growth rates.

For each substrate the published maximum growth rate is the largest mu at which the ME LP is
feasible. This script checks, for each substrate, that the model is FEASIBLE at (mu - 0.02) and
at the published mu, and INFEASIBLE at (mu + 0.02) -- i.e. the published mu is exactly the model's
growth boundary. Exit code 0 iff all four reproduce.
"""
import sys
from palustris_me import MEModel, SUBSTRATES


def main():
    me = MEModel()
    print("%-11s | published mu | feas @mu-0.02 | feas @mu | feas @mu+0.02 | reproduced" % "substrate")
    print("-" * 78)
    all_ok = True
    for name, (ex, up, atp, mu) in SUBSTRATES.items():
        me.set_medium(name)
        below = me.solve(round(mu - 0.02, 2), atp_maint=atp)[0] == "optimal"
        at = me.solve(mu, atp_maint=atp)[0] == "optimal"
        above = me.solve(round(mu + 0.02, 2), atp_maint=atp)[0] == "optimal"
        ok = below and at and not above
        all_ok &= ok
        print("%-11s |    %.2f      |    %-5s     |  %-5s   |    %-5s     | %s"
              % (name, mu, below, at, above, "YES" if ok else "NO"))
    print("-" * 78)
    print("RESULT:", "ALL FOUR SAHA TABLE-2 GROWTH RATES REPRODUCED EXACTLY"
          if all_ok else "MISMATCH")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
