"""Re-run the pre-registered robustness test (Nd = 3, 4; Gamma and 2x2 k-mesh;
U/t = 1.15) and save the full table to results/robustness_test.json.

Usage (from the repository root):
    python scripts/run_robustness_test.py
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from holey_graphene_qubit import robustness_test

def _clean(v):
    if isinstance(v, (np.floating, float)):
        v = float(v); return None if not np.isfinite(v) else v
    if isinstance(v, (np.integer,)): return int(v)
    if isinstance(v, np.bool_): return bool(v)
    if isinstance(v, np.ndarray): return None  # drop arrays (masks, moments)
    if isinstance(v, (list, tuple)): return [_clean(x) for x in v]
    if isinstance(v, dict): return {k: _clean(x) for k, x in v.items()}
    return v

if __name__ == "__main__":
    t0 = time.time()
    out = robustness_test(nk=2, Nds=(3, 4), uot_band=(1.15,))
    rows = []
    for (Nd, knk), d in out["table"].items():
        for rep, rec in d.items():
            r = {k: _clean(v) for k, v in rec.items() if not isinstance(v, np.ndarray)}
            r.update(Nd=Nd, k_mesh=knk)
            rows.append(r)
    res = {"settings": {"Nds": [3, 4], "nk": 2, "U_over_t": [1.15]},
           "rows": rows, "robust": out["robust"],
           "flips": _clean(out["flips"]), "dispersive": _clean(out["dispersive"]),
           "wall_time_s": round(time.time() - t0, 1)}
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs(os.path.join(root, "results"), exist_ok=True)
    with open(os.path.join(root, "results", "robustness_test.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("saved results/robustness_test.json  wall %.1f s" % res["wall_time_s"])
