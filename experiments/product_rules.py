"""Residual widths of the multiplication rules on the boxes of results/<set>.json.

    python experiments/product_rules.py            # writes results/product_rules.json and prints a summary

All the rules share the representation and the linear part of a product (see zcaa/zcaa.py)
and differ only in the residual interval:

    faf2    'standard' multiplication of reduced affine forms (Skalna and Hladik 2017, Eq. 5),
            the rule of the fAF2 model of Ibex
    zcaa    monomial-wise envelopes + four-vertex residual (Section 3 of the paper)
    skalna  Chebyshev minimum-error multiplication (Skalna and Hladik 2017, Eq. 10)
    exact   exact range of the nonlinear part (four vertices x joint range)

AF2 and first-order Taylor models (TM) are taken from the stored records, whose W_OPT is
computed with the linear part of ZCAA. In a single product the four rules have the same
linear part; in compositions (Case 4, Pramanik) the linear parts may differ, and the optimal
width of each rule is then recomputed with IbexOpt for its own linear part (<rule>_opt_width),
with the precision and time limit of run_experiments.py.
"""
import json
import os

import numpy as np

from common import RESULTS, optimal_residual_width
from run_experiments import CONFIG
from zcaa import set_box, zcaa_variables, set_product
from benchmarks.functions import ALL

RULES = ["faf2", "zcaa", "skalna", "exact"]
BY_NAME = {f.name: f for f in ALL}


def ibex_settings(set_name):
    """(precision, timeout) used for each function of the set in run_experiments.py."""
    out = {}
    for item in CONFIG[set_name]["functions"]:
        f, prec, timeout = item[0], item[-2], item[-1]
        out[f.name] = (prec, timeout)
    return out


def evaluate_rules(records, settings):
    out = []
    for r in records:
        f = BY_NAME[r["function"]]
        set_box(r["box"])
        rec = dict(r)
        slopes = {}
        for rule in RULES:
            set_product(rule)
            form = f(*zcaa_variables(3))
            rec[f"{rule}_width"] = form.residual_width()
            slopes[rule] = form.slopes
        set_product("zcaa")
        ref = np.array(slopes["zcaa"])
        computed = {}
        for rule in RULES:
            s = np.array(slopes[rule])
            if np.max(np.abs(s - ref)) <= 1e-9 * (1 + np.max(np.abs(ref))):
                rec[f"{rule}_opt_width"] = r["opt_width"]
                continue
            key = tuple(np.round(s, 12))
            if key not in computed:
                computed[key] = optimal_residual_width(f, [float(t) for t in s], r["box"], *settings[f.name])
            rec[f"{rule}_opt_width"] = computed[key]
        out.append(rec)
    return out


def summary(records):
    print(f"{'Function':16s} {'Size':6s} " + " ".join(f"{m:>8s}" for m in ["TM", "AF2"] + RULES)
          + "   exact<zcaa  zcaa<skalna  exact/skalna")
    keys = list(dict.fromkeys((r["function"], r["size"]) for r in records))
    for f, size in keys:
        sel = [r for r in records if r["function"] == f and r["size"] == size]
        cols = []
        for m in ["tm", "af2"] + RULES:
            w = np.array([r[f"{m}_width"] for r in sel])
            wo = np.array([r.get(f"{m}_opt_width", r["opt_width"]) for r in sel])
            cols.append(np.mean((w - wo) / wo * 100))
        we = np.array([r["exact_width"] for r in sel])
        wz = np.array([r["zcaa_width"] for r in sel])
        ws = np.array([r["skalna_width"] for r in sel])
        tol = 1e-9 * np.maximum(wz, ws)
        print(f"{f:16s} {size:6s} " + " ".join(f"{c:8.2f}" for c in cols)
              + f"   {np.sum(we < wz - tol):5d}/{len(sel):<4d} {np.sum(wz < ws - tol):5d}/{len(sel):<4d}"
              + f"   {np.mean(we / ws):.3f}")


if __name__ == "__main__":
    allrecs = {}
    for name in ("controlled", "mixed", "instances"):
        path = os.path.join(RESULTS, f"{name}.json")
        if not os.path.exists(path):
            print(f"{path} not found, skipped")
            continue
        with open(path) as fh:
            recs = evaluate_rules(json.load(fh), ibex_settings(name))
        allrecs[name] = recs
        print(f"=== {name}: mean gap (%) w.r.t. the optimal width of each method ===")
        summary(recs)
        print()
    with open(os.path.join(RESULTS, "product_rules.json"), "w") as fh:
        json.dump(allrecs, fh, indent=0)
