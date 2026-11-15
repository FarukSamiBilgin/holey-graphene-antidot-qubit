"""
holey_graphene_qubit.py
=======================

Backbone for the project we scoped:

    "Kasitli kusurla tanimlanmis delikli-grafen qubitinin koherans-belirleyici
     ozellikleri ile parazit kenar/midgap momentini, AYNI geometri taramasinda
     karsilastirip, bu iki rakip buyuklugun ortustugu mu yoksa odunlestigi mi
     rejimini haritalamak."

This file is intentionally ONE parameter-driven pipeline so we write it once and
open the scope incrementally (triangle -> hexagon -> circle ; triangular ->
square lattice ; single defect -> double-dot) by changing arguments only.

Pipeline layers (matching what we agreed):
  L1  geometry generator (honeycomb + antidot carving + designed defect)
  L1b Lieb sublattice-imbalance pre-filter           [analytic, instant]
  L2  sparse tight-binding, states near Fermi level   [shift-invert eigsh]
  L2m competing metrics, per geometry:
        - w_edge   : near-zero weight on antidot-edge atoms   (PARASITIC, min)
        - w_defect : near-zero weight on the designed defect   (QUBIT, max)
        - iso      : energy isolation of the defect state in the gap (max)
        - Jproxy   : (double-dot only) splitting of the two defect states

`evaluate_geometry(params) -> dict` is a PURE function: joblib/multiprocessing
can map it over a parameter grid for the parallel scan. No GPU needed; the heavy
call is scipy.sparse.linalg.eigsh with shift-invert.

Requires only numpy + scipy (both preinstalled on Colab; pick a standard CPU
runtime, NOT a GPU runtime).

COLAB USAGE  (do NOT run this file as a script in a notebook cell - the
Jupyter kernel injects its own args and argparse will SystemExit). Instead:

    # cell 1: write the file to disk (paste the whole module), OR upload it
    # cell 2: import and call functions directly
    from holey_graphene_qubit import periodic_validation, evaluate_geometry
    periodic_validation()                       # the two mandatory checks
    periodic_validation(shape="triangle", edge="zigzag")   # parasitic side

Or from a shell cell:  !python holey_graphene_qubit.py --periodic

Run as a script (terminal) for the same validation; --demo is the legacy
finite-flake path kept only for reference.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field, asdict
from typing import Callable

import numpy as np
from scipy.spatial import cKDTree
from scipy import sparse
from scipy.sparse.linalg import eigsh


# ----------------------------------------------------------------------------
# Units
# ----------------------------------------------------------------------------
# All lengths in units of the C-C bond length a_cc = 1.
# Nearest-neighbour hopping magnitude (eV); energies come out in eV.
T_HOP = 2.7


# ----------------------------------------------------------------------------
# L1 - honeycomb sheet
# ----------------------------------------------------------------------------
def make_honeycomb(Lx: float, Ly: float):
    """Generate a finite honeycomb (graphene) patch in the box [0,Lx] x [0,Ly].

    Returns
    -------
    pos : (N,2) float array of atom positions
    sub : (N,)  int array, sublattice label 0 (A) or 1 (B)
    """
    a1 = np.array([1.5,  np.sqrt(3) / 2.0])   # primitive vectors (a_cc = 1)
    a2 = np.array([1.5, -np.sqrt(3) / 2.0])
    basis = [np.array([0.0, 0.0]),            # A
             np.array([1.0, 0.0])]            # B

    # integer range generous enough to cover the box after the shear
    nmax = int(np.ceil((Lx + Ly) / 1.5)) + 4
    mmax = nmax
    ns = np.arange(-nmax, nmax + 1)
    ms = np.arange(-mmax, mmax + 1)
    NN, MM = np.meshgrid(ns, ms, indexing="ij")
    NN = NN.ravel()
    MM = MM.ravel()

    pos_list = []
    sub_list = []
    for b_idx, b in enumerate(basis):
        P = (NN[:, None] * a1[None, :] +
             MM[:, None] * a2[None, :] + b[None, :])
        pos_list.append(P)
        sub_list.append(np.full(P.shape[0], b_idx, dtype=np.int8))

    pos = np.vstack(pos_list)
    sub = np.concatenate(sub_list)

    inside = (pos[:, 0] >= 0) & (pos[:, 0] <= Lx) & \
             (pos[:, 1] >= 0) & (pos[:, 1] <= Ly)
    return pos[inside], sub[inside]


# ----------------------------------------------------------------------------
# L1 - antidot shape predicates  (parameter-driven: this is the "write once"
#       part; which shapes we actually RUN is our methodological choice)
# ----------------------------------------------------------------------------
def _inside_circle(p, center, R):
    return np.linalg.norm(p - center, axis=1) <= R


def _inside_regular_polygon(p, center, R, n_sides, phase):
    """Point-in-regular-polygon (apothem-form). R = circumradius."""
    d = p - center
    ang = np.arctan2(d[:, 1], d[:, 0])
    r = np.linalg.norm(d, axis=1)
    # distance to nearest polygon edge along radial direction
    sector = 2.0 * np.pi / n_sides
    # fold angle into one sector measured from a vertex direction (phase)
    a = np.mod(ang - phase, sector)
    a = np.minimum(a, sector - a)
    # radius of polygon boundary at this angle
    r_bound = R * np.cos(np.pi / n_sides) / np.maximum(np.cos(a), 1e-12)
    return r <= r_bound


def _chi_of(edge):
    """Map `edge` to a chirality fraction chi in [0,1].
    'zigzag'->0.0, 'armchair'->1.0, or pass a float in [0,1] directly for an
    INTERMEDIATE chiral edge (the transition region the scan must probe)."""
    if isinstance(edge, str):
        return {"zigzag": 0.0, "armchair": 1.0}[edge]
    chi = float(edge)
    if not 0.0 <= chi <= 1.0:
        raise ValueError("edge chirality fraction must be in [0,1]")
    return chi


def shape_predicate(shape: str, edge):
    """Return f(pos, center, R)->bool mask of atoms inside one antidot.

    EQUAL-AREA NORMALISATION: R is a reference-circle radius; every shape is
    scaled to that SAME area so shapes differ by SHAPE, not accidental size
    (this fixed the small-R collapse where different shapes returned
    identical N/|dN|).

    CONTINUOUS CHIRALITY: `edge` is either 'zigzag'/'armchair' or a float
    chi in [0,1]. The hole is rotated by phase = chi * (pi/6), interpolating
    the boundary orientation from zigzag-type (chi=0) to armchair-type
    (chi=1). This is an orientation PROXY for edge chirality - not an exact
    crystallographic edge-vector construction - but it is exactly what is
    needed to scan the transition region, where any parasitic-vs-qubit
    trade-off (absent at the two well-known extremes) must reveal itself.
    """
    if shape == "circle":
        return lambda p, c, R: _inside_circle(p, c, R)

    if shape in ("triangle", "hexagon"):
        n = 3 if shape == "triangle" else 6
        phase = _chi_of(edge) * (np.pi / 6.0)
        scale = np.sqrt(2.0 * np.pi / (n * np.sin(2.0 * np.pi / n)))
        return lambda p, c, R: _inside_regular_polygon(
            p, c, R * scale, n, phase)

    raise ValueError(f"unknown shape: {shape}")


# ----------------------------------------------------------------------------
# L1 - antidot-lattice hole centers (the superlattice)
# ----------------------------------------------------------------------------
def hole_centers(lattice: str, L: float, Lx: float, Ly: float):
    """Centers of the periodic hole array. L = superlattice constant (a_cc)."""
    centers = []
    if lattice == "triangular":
        v1 = np.array([L, 0.0])
        v2 = np.array([L / 2.0, L * np.sqrt(3) / 2.0])
    elif lattice == "square":
        v1 = np.array([L, 0.0])
        v2 = np.array([0.0, L])
    else:
        raise ValueError(f"unknown lattice: {lattice}")

    imax = int(np.ceil(Lx / L)) + 2
    jmax = int(np.ceil(Ly / (L * 0.5))) + 2
    for i in range(-2, imax):
        for j in range(-2, jmax):
            c = i * v1 + j * v2
            if -L <= c[0] <= Lx + L and -L <= c[1] <= Ly + L:
                centers.append(c)
    return np.array(centers)


# ----------------------------------------------------------------------------
# L1 - assemble geometry: carve antidots, leave designed defect(s) uncarved
# ----------------------------------------------------------------------------
def _prune_dangling(pos: np.ndarray, sub: np.ndarray, min_coord: int = 2):
    """Iteratively drop atoms with < min_coord neighbours.

    Singly-coordinated carbon is unphysical and injects spurious exact zero
    modes that also wreck shift-invert conditioning. Standard for antidot TB.
    """
    while True:
        tree = cKDTree(pos)
        pairs = tree.query_pairs(r=1.1, output_type="ndarray")
        deg = np.zeros(len(pos), dtype=np.int64)
        if len(pairs):
            np.add.at(deg, pairs[:, 0], 1)
            np.add.at(deg, pairs[:, 1], 1)
        keep = deg >= min_coord
        if keep.all():
            return pos, sub
        pos, sub = pos[keep], sub[keep]
        if len(pos) == 0:
            return pos, sub


@dataclass
class Geometry:
    pos: np.ndarray
    sub: np.ndarray
    edge_atom_mask: np.ndarray     # atoms bordering a carved hole (parasitic)
    defect_atom_mask: np.ndarray   # atoms in the designed-defect region (qubit)
    dN: int                        # global sublattice imbalance N_A - N_B
    meta: dict = field(default_factory=dict)


def build_geometry(params: dict) -> Geometry:
    p = params
    Lx, Ly = p["Lx"], p["Ly"]
    pos, sub = make_honeycomb(Lx, Ly)

    centers = hole_centers(p["lattice"], p["L"], Lx, Ly)
    pred = shape_predicate(p["shape"], p["edge"])
    R = p["R"]

    # choose the designed defect: omit the hole nearest the box center
    # (single) or the two holes straddling the center (double-dot).
    box_c = np.array([Lx / 2.0, Ly / 2.0])
    order = np.argsort(np.linalg.norm(centers - box_c, axis=1))
    if p["defect"] == "single":
        omit = {order[0]}
    elif p["defect"] == "double":
        omit = {order[0], order[1]}
    else:
        raise ValueError("defect must be 'single' or 'double'")

    active_centers = centers[[k for k in range(len(centers))
                              if k not in omit]]
    defect_centers = centers[list(omit)]

    # carve: remove atoms inside any ACTIVE hole
    remove = np.zeros(len(pos), dtype=bool)
    for c in active_centers:
        remove |= pred(pos, c, R)
    pos_k, sub_k = pos[~remove], sub[~remove]

    # remove unphysical dangling atoms created by carving
    pos_k, sub_k = _prune_dangling(pos_k, sub_k)
    if len(pos_k) == 0:
        return Geometry(pos_k, sub_k,
                        np.zeros(0, bool), np.zeros(0, bool), 0,
                        {"N": 0})

    # edge atoms: rim sites just OUTSIDE a carved hole (pruning-safe:
    # defined purely by distance to active centers, no index bookkeeping)
    edge_atom_mask = np.zeros(len(pos_k), dtype=bool)
    for c in active_centers:
        d = np.linalg.norm(pos_k - c, axis=1)
        edge_atom_mask |= (d >= R) & (d <= R + 1.3)

    # defect region: kept sites within R of an omitted center -> exactly the
    # sites that host the designed defect state
    defect_atom_mask = np.zeros(len(pos_k), dtype=bool)
    for c in defect_centers:
        defect_atom_mask |= (np.linalg.norm(pos_k - c, axis=1) <= R)

    dN = int(np.sum(sub_k == 0) - np.sum(sub_k == 1))

    return Geometry(
        pos=pos_k, sub=sub_k,
        edge_atom_mask=edge_atom_mask,
        defect_atom_mask=defect_atom_mask,
        dN=dN,
        meta={"N": int(len(pos_k)),
              "n_holes_active": int(len(active_centers)),
              "n_holes_defect": int(len(defect_centers))},
    )


# ----------------------------------------------------------------------------
# L2 - sparse tight-binding Hamiltonian + near-Fermi spectrum
# ----------------------------------------------------------------------------
def tight_binding_H(pos: np.ndarray) -> sparse.csr_matrix:
    tree = cKDTree(pos)
    pairs = tree.query_pairs(r=1.1, output_type="ndarray")
    N = len(pos)
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    vals = np.full(len(rows), -T_HOP, dtype=np.float64)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))


def _det_v0(N: int, seed: int = 0) -> np.ndarray:
    """Fixed deterministic ARPACK start vector.

    eigsh with shift-invert starts from a RANDOM vector unless v0 is given.
    With many exactly-degenerate zero modes (zigzag holes -> a flat band at
    E=0) the returned degenerate subspace - and any quantity read off the
    near-zero states - depends on that random seed, so the same input gave
    different output across runs. A fixed v0 makes every call bit-identical.
    """
    return np.random.default_rng(seed).standard_normal(N)


def _robust_gap(vals: np.ndarray, ztol: float = 1e-6):
    """Band gap that is INVARIANT to how the degenerate E=0 manifold splits.

    The exact zero modes (flat-band edge states) are physical for zigzag
    holes, not an artifact - but mixing them into occ.max()/vir.min() made
    the gap seed-dependent. Define the gap strictly between the dispersive
    bands (|E| > ztol) and report the zero-mode count separately.
    """
    nz = int(np.sum(np.abs(vals) <= ztol))
    pos = vals[vals > ztol]
    neg = vals[vals < -ztol]
    gap = float(pos.min() - neg.max()) if pos.size and neg.size \
        else float("nan")
    return gap, nz


def near_fermi_spectrum(H: sparse.csr_matrix, k: int = 12):
    """Eigenpairs closest to E = 0 via shift-invert.

    Bipartite TB has EXACT zero modes, so sigma=0 makes (H - sigma I) singular
    and the LU/shift-invert factorisation fails. A tiny offset (1e-3 eV, far
    inside any antidot gap of ~0.1-1 eV) restores conditioning while still
    returning the near-zero defect/edge midgap states. Deterministic v0.
    """
    k = min(k, H.shape[0] - 2)
    vals, vecs = eigsh(H, k=k, sigma=1e-3, which="LM",
                       v0=_det_v0(H.shape[0]))
    order = np.argsort(np.abs(vals))
    return vals[order], vecs[:, order]


# ----------------------------------------------------------------------------
# L2m - competing metrics  (the scientific core)
# ----------------------------------------------------------------------------
def evaluate_geometry(params: dict) -> dict:
    """PURE function: geometry params -> metrics. Safe to map in parallel."""
    g = build_geometry(params)
    out = {**params, **g.meta, "dN": g.dN}

    # L1b Lieb pre-filter: # of guaranteed zero modes >= |dN|.
    # Large |dN| => guaranteed net moment => parasitic. Cheap reject signal.
    out["lieb_imbalance"] = abs(g.dN)

    if g.meta["N"] < 8:
        out["status"] = "too_small"
        return out

    try:
        vals, vecs = near_fermi_spectrum(tight_binding_H(g.pos),
                                         k=params.get("k", 12))
    except Exception as e:                      # pragma: no cover
        out["status"] = f"eig_fail:{type(e).__name__}"
        return out

    prob = vecs ** 2                            # |psi|^2 columns (real H)
    prob /= prob.sum(axis=0, keepdims=True)

    # weight of each near-zero state on parasitic edge vs designed defect
    w_edge_all = prob[g.edge_atom_mask, :].sum(axis=0)
    w_def_all = prob[g.defect_atom_mask, :].sum(axis=0)
    ipr_all = (prob ** 2).sum(axis=0)           # inverse participation ratio

    # designed-defect state = near-zero state most localized on the defect
    j = int(np.argmax(w_def_all))
    E_defect = float(vals[j])

    # gap proxy: smallest |E| among states that are NOT defect/edge-localised
    bulk = np.array([i for i in range(len(vals))
                     if (w_edge_all[i] < 0.25 and w_def_all[i] < 0.25)])
    gap_proxy = float(np.min(np.abs(vals[bulk]))) if bulk.size else float("nan")

    out.update({
        "status": "ok",
        "w_edge":   float(w_edge_all[j]),   # PARASITIC  (minimise)
        "w_defect": float(w_def_all[j]),    # QUBIT      (maximise)
        "ipr_defect": float(ipr_all[j]),    # localisation of qubit state
        "E_defect": E_defect,
        "gap_proxy": gap_proxy,             # bulk gap edge (approx, finite flake)
        "isolation": gap_proxy - abs(E_defect) if np.isfinite(gap_proxy)
                     else float("nan"),
        "spectrum_near0": [float(v) for v in vals[:8]],
    })

    if params["defect"] == "double":
        # two lowest |E| defect-like states -> tunnel/exchange proxy
        dlike = np.argsort(-w_def_all)[:2]
        out["Jproxy"] = float(abs(vals[dlike[0]] - vals[dlike[1]]))

    return out


# ----------------------------------------------------------------------------
# Demo + mini scan
# ----------------------------------------------------------------------------
def _base_params() -> dict:
    return dict(
        Lx=46.0, Ly=40.0,      # flake box (a_cc) - keep modest for first runs
        lattice="triangular",
        shape="triangle",
        edge="zigzag",
        L=12.0,                # superlattice constant
        R=3.2,                 # antidot size
        defect="single",
        k=14,
    )


def _demo():
    print("=" * 70)
    print("DEMO  triangle antidot / triangular lattice / single defect")
    print("=" * 70)
    for edge in ("zigzag", "armchair"):
        p = _base_params()
        p["edge"] = edge
        r = evaluate_geometry(p)
        if r.get("status") != "ok":
            print(f"edge={edge:8s}  status={r.get('status')}  "
                  f"N={r.get('N')}  |dN|={r.get('lieb_imbalance')}")
            continue
        print(f"edge={edge:8s}  N={r['N']:4d}  |dN|={r['lieb_imbalance']:2d}  "
              f"w_edge(parasitic)={r['w_edge']:.3f}  "
              f"w_def(qubit)={r['w_defect']:.3f}  "
              f"IPR={r['ipr_defect']:.3f}  "
              f"|E_def|={abs(r['E_defect']):.3f}  "
              f"gap~{r['gap_proxy']:.3f}")

    print("\nMINI SCAN  edge=zigzag, sweep antidot size R "
          "(parasitic vs qubit trade-off probe)")
    print("-" * 70)
    for R in (2.2, 2.8, 3.4, 4.0, 4.6):
        p = _base_params()
        p["R"] = R
        r = evaluate_geometry(p)
        if r.get("status") != "ok":
            print(f"R={R:.1f}  status={r.get('status')}")
            continue
        print(f"R={R:.1f}  |dN|={r['lieb_imbalance']:2d}  "
              f"w_edge={r['w_edge']:.3f}  w_def={r['w_defect']:.3f}  "
              f"IPR={r['ipr_defect']:.3f}  iso={r['isolation']:.3f}")


# ============================================================================
# PERIODIC COMMENSURATE SUPERCELL  (current scientific path)
# ----------------------------------------------------------------------------
# The finite-flake path above showed the rectangular flake's own boundary
# swamps the defect state. Pedersen 2008/2013 avoid this with a periodic
# antidot supercell. Here:
#   - graphene Bravais vectors ga1, ga2 (a_cc = 1), 2-atom basis
#   - antidot cell      = L x L graphene cells, one hole at the cell centre
#   - defect supercell  = Nd x Nd antidot cells with the CENTRAL hole OMITTED
#                         (exactly Pedersen's designed-defect construction)
# Everything is integer-commensurate by construction, so PBC are exact.
# ============================================================================

GA1 = np.array([1.5,  np.sqrt(3) / 2.0])   # graphene Bravais vectors
GA2 = np.array([1.5, -np.sqrt(3) / 2.0])
BASIS = (np.array([0.0, 0.0]), np.array([1.0, 0.0]))   # A, B


def _supercell_atoms(L: int, Nd: int):
    """Atoms tiling one Nd*L x Nd*L graphene-cell supercell."""
    M = Nd * L
    ns, ms = np.meshgrid(np.arange(M), np.arange(M), indexing="ij")
    ns, ms = ns.ravel(), ms.ravel()
    pos, sub = [], []
    for b_idx, b in enumerate(BASIS):
        pos.append(ns[:, None] * GA1 + ms[:, None] * GA2 + b)
        sub.append(np.full(ns.size, b_idx, np.int8))
    S1 = M * GA1
    S2 = M * GA2
    return np.vstack(pos), np.concatenate(sub), S1, S2


def _frac(disp, S1, S2):
    """Cartesian -> fractional coords in the (S1,S2) supercell basis."""
    Bmat = np.linalg.inv(np.column_stack([S1, S2]))
    return disp @ Bmat.T


def _min_image(disp, S1, S2):
    f = _frac(disp, S1, S2)
    f -= np.round(f)
    return f @ np.column_stack([S1, S2]).T


def build_periodic_geometry(shape, edge, L, R, Nd, defect):
    pos, sub, S1, S2 = _supercell_atoms(L, Nd)

    # antidot-cell centres inside the supercell
    centers = np.array([(i + 0.5) * (L * GA1) + (j + 0.5) * (L * GA2)
                        for i in range(Nd) for j in range(Nd)])
    cidx = list(range(len(centers)))
    cell_c = 0.5 * (S1 + S2)
    order = np.argsort(np.linalg.norm(centers - cell_c, axis=1))
    omit = {order[0]} if defect == "single" else {order[0], order[1]}
    active = centers[[k for k in cidx if k not in omit]]
    defc = centers[list(omit)]

    pred = shape_predicate(shape, edge)
    remove = np.zeros(len(pos), bool)
    for c in active:                       # min-image so holes wrap correctly
        d = pos - c
        dm = _min_image(d, S1, S2)
        remove |= pred(dm + c, c, R)
    pos, sub = pos[~remove], sub[~remove]

    # PBC dangling-atom prune
    for _ in range(12):
        pairs, _off = _pbc_pairs(pos, S1, S2)
        deg = np.zeros(len(pos), np.int64)
        if len(pairs):
            np.add.at(deg, pairs[:, 0], 1)
            np.add.at(deg, pairs[:, 1], 1)
        keep = deg >= 2
        if keep.all():
            break
        pos, sub = pos[keep], sub[keep]

    # edge atoms: rim sites just OUTSIDE a carved (active) hole
    edge_mask = np.zeros(len(pos), bool)
    for c in active:
        dm = np.linalg.norm(_min_image(pos - c, S1, S2), axis=1)
        edge_mask |= (dm >= R) & (dm <= R + 1.3)

    # defect region = Voronoi cell of the OMITTED hole(s): an atom belongs to
    # the designed defect if its nearest antidot centre (min-image) is an
    # omitted one. This captures the whole missing-hole cell where the defect
    # state actually lives - radius R undersized it (the earlier 'still
    # buried' readout) since the state spreads over the ~L-scale cell.
    all_c = np.vstack([active, defc]) if len(active) else defc
    is_def_center = np.array([False] * len(active) + [True] * len(defc))
    d_all = np.stack([np.linalg.norm(_min_image(pos - c, S1, S2), axis=1)
                      for c in all_c], axis=1)
    nearest = np.argmin(d_all, axis=1)
    def_mask = is_def_center[nearest] & ~edge_mask

    dN = int(np.sum(sub == 0) - np.sum(sub == 1))
    return dict(pos=pos, sub=sub, S1=S1, S2=S2,
                edge_mask=edge_mask, def_mask=def_mask, dN=dN)


def _pbc_pairs(pos, S1, S2):
    """Neighbour pairs under PBC; returns (pairs, offsets) with offsets the
    integer (i,j) image of the second atom in (S1,S2) units."""
    imgs = []
    offs = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            imgs.append(pos + i * S1 + j * S2)
            offs.append((i, j))
    home = cKDTree(pos)
    pr_list = []
    off_list = []
    for shift, (i, j) in zip(imgs, offs):
        tree = cKDTree(shift)
        res = home.query_ball_tree(tree, r=1.1)
        for a, bs in enumerate(res):
            for b in bs:
                if (i, j) == (0, 0) and a >= b:
                    continue
                pr_list.append((a, b))
                off_list.append((i, j))
    if not pr_list:
        return np.zeros((0, 2), int), np.zeros((0, 2), int)
    return np.array(pr_list), np.array(off_list)


def _recip(S1, S2):
    """2D reciprocal vectors: S_i . B_j = 2 pi delta_ij."""
    Smat = np.column_stack([S1, S2])              # columns S1,S2
    B = 2.0 * np.pi * np.linalg.inv(Smat).T       # rows B1,B2
    return B[0], B[1]


def _kmesh(S1, S2, nk):
    """Gamma-centred nk x nk k-point mesh in the supercell BZ.
    nk=1 -> just Gamma (recovers the original behaviour exactly)."""
    if nk <= 1:
        return [np.zeros(2)]
    B1, B2 = _recip(S1, S2)
    fr = (np.arange(nk) + 0.5) / nk - 0.5         # symmetric, Gamma-centred
    return [a * B1 + b * B2 for a in fr for b in fr]


def _bloch_H(pos, S1, S2, k):
    pairs, offs = _pbc_pairs(pos, S1, S2)
    N = len(pos)
    if len(pairs) == 0:
        return sparse.csr_matrix((N, N), dtype=complex)
    Rvec = offs[:, 0:1] * S1 + offs[:, 1:2] * S2
    phase = np.exp(1j * (Rvec @ k))
    a, b = pairs[:, 0], pairs[:, 1]
    data = np.concatenate([-T_HOP * phase, -T_HOP * np.conj(phase)])
    rows = np.concatenate([a, b])
    cols = np.concatenate([b, a])
    return sparse.csr_matrix((data, (rows, cols)), shape=(N, N))


def _near0(H, k=16):
    k = min(k, H.shape[0] - 2)
    vals, vecs = eigsh(H, k=k, sigma=1e-3, which="LM",
                       v0=_det_v0(H.shape[0]))
    o = np.argsort(np.abs(vals))
    return vals[o], vecs[:, o]


# ============================================================================
# L3 - MEAN-FIELD HUBBARD  (the parasitic-moment axis becomes physical here)
# ----------------------------------------------------------------------------
# At pure TB level the zigzag zero modes are just degenerate states - NOT a
# magnetic moment. The moment (hence the whole "parasitic edge moment vs
# qubit" competition our project is about) only appears once on-site
# Coulomb U is included. Collinear pi-orbital mean-field Hubbard:
#
#     H_sigma = H_TB + U * diag( <n_{i,-sigma}> - 1/2 )
#
# solved self-consistently at half filling. Result: local moment
# m_i = <n_{i,up}> - <n_{i,dn}>.  Reported resolved onto the parasitic edge
# mask and the designed-defect mask - exactly the two competing quantities.
# MF Hubbard is U-sensitive, so every call sweeps a U/t band, never one
# fragile value.
# ============================================================================

U_OVER_T_BAND = (1.0, 1.15, 1.3)


def _two_color(H0):
    """2-colour the (bipartite) adjacency of H0 -> sublattice sign +/-1."""
    N = H0.shape[0]
    A = (np.abs(H0) > 1e-9)
    color = np.zeros(N, dtype=np.int8)
    seen = np.zeros(N, bool)
    for s in range(N):
        if seen[s]:
            continue
        color[s] = 1
        seen[s] = True
        stack = [s]
        while stack:
            u = stack.pop()
            for v in np.where(A[u])[0]:
                if not seen[v]:
                    color[v] = -color[u]
                    seen[v] = True
                    stack.append(v)
    return color.astype(float)


def mean_field_hubbard(pos, S1, S2, U, max_iter=400, tol=1e-5,
                       mix=0.3, m0=0.5, seed=0):
    """Self-consistent collinear MF Hubbard at Gamma, half filling.

    Gamma suffices here: edge/defect moments are spatially localized
    (flat-band magnetism), so the k=0 cell captures them; a k-mesh
    refinement is a deliberate later step.
    """
    H0 = np.real(_bloch_H(pos, S1, S2, k=np.zeros(2)).toarray())
    N = H0.shape[0]
    Ne = N                                  # one pi electron per site

    rng = np.random.default_rng(seed)
    sub_sign = _two_color(H0)               # staggered AFM seed
    n_up = 0.5 + 0.5 * m0 * sub_sign + 0.01 * rng.standard_normal(N)
    n_dn = 0.5 - 0.5 * m0 * sub_sign + 0.01 * rng.standard_normal(N)

    converged = False
    it = 0
    for it in range(max_iter):
        Hu = H0 + np.diag(U * (n_dn - 0.5))
        Hd = H0 + np.diag(U * (n_up - 0.5))
        eu, vu = np.linalg.eigh(Hu)
        ed, vd = np.linalg.eigh(Hd)

        # common chemical potential: fill Ne lowest levels across BOTH spins
        allE = np.concatenate([eu, ed])
        idx = np.argsort(allE)[:Ne]
        occ_u = idx[idx < N]                 # columns of vu
        occ_d = idx[idx >= N] - N            # columns of vd

        new_up = (np.abs(vu[:, occ_u]) ** 2).sum(axis=1) \
            if occ_u.size else np.zeros(N)
        new_dn = (np.abs(vd[:, occ_d]) ** 2).sum(axis=1) \
            if occ_d.size else np.zeros(N)

        delta = max(np.max(np.abs(new_up - n_up)),
                    np.max(np.abs(new_dn - n_dn)))
        n_up = (1 - mix) * n_up + mix * new_up
        n_dn = (1 - mix) * n_dn + mix * new_dn
        if delta < tol:
            converged = True
            break

    m = n_up - n_dn
    return {"m": m, "iters": it + 1, "converged": converged,
            "Mabs": float(np.sum(np.abs(m))),
            "Mnet": float(np.sum(m))}


def _moment_metrics(m, edge_mask, def_mask):
    """Correct parasitic/qubit moment metrics from the MF moment field m_i.

    Sum|m_i| over a whole cell was wrong: it adds staggered AFM fluctuations
    so even a moment-free balanced cell scored ~2. The physical quantities:
      - mu_edge : MEAN |m_i| per rim atom  (Vanevic's parasitic edge moment;
                  size-independent, ~0 for armchair, ~0.2-0.4 for zigzag)
      - mu_def  : MEAN |m_i| per defect-region atom (moment intruding on the
                  qubit)
      - Mnet    : |sum m_i|  (global net spin; Lieb lower bound = |dN|)
    """
    ne = int(edge_mask.sum())
    nd = int(def_mask.sum())
    mu_edge = float(np.sum(np.abs(m[edge_mask])) / ne) if ne else 0.0
    mu_def = float(np.sum(np.abs(m[def_mask])) / nd) if nd else 0.0
    return {"mu_edge": mu_edge, "mu_def": mu_def,
            "Mnet": float(abs(np.sum(m))), "n_edge": ne, "n_def": nd}


def hubbard_validation(L=8, R=3.2, Nd=3):
    """MANDATORY check for the new layer, now on the DEFECTED supercell where
    the edge/defect masks are meaningful (pristine had no masks; sum|m_i|
    over it was the wrong metric). Decisive, U/t-robust ordering required:
      balanced hexagon/armchair -> rim is armchair -> mu_edge ~ 0
      zigzag triangle           -> zigzag rim      -> mu_edge large
    Else the MF/metric is wrong and we stop before the scan.
    """
    print("=" * 70)
    print(f"HUBBARD VALIDATION  defected supercell  L={L} R={R} Nd={Nd}  "
          f"U/t in {U_OVER_T_BAND}")
    print("=" * 70)
    cases = [("hexagon", "armchair", "balanced rim -> expect mu_edge ~ 0"),
             ("triangle", "zigzag", "zigzag rim  -> expect mu_edge large")]
    results = {}
    for shape, edge, note in cases:
        g = build_periodic_geometry(shape, edge, L, R, Nd=Nd,
                                    defect="single")
        print(f"\n{shape}/{edge:8s}  N={len(g['pos'])}  |dN|={abs(g['dN'])}"
              f"  n_edge={int(g['edge_mask'].sum())}"
              f"  n_def={int(g['def_mask'].sum())}   ({note})")
        row = {}
        for uot in U_OVER_T_BAND:
            r = mean_field_hubbard(g["pos"], g["S1"], g["S2"],
                                   U=uot * T_HOP)
            mm = _moment_metrics(r["m"], g["edge_mask"], g["def_mask"])
            row[uot] = mm
            flag = "" if r["converged"] else "  *NOT CONVERGED*"
            print(f"   U/t={uot:4.2f}  iters={r['iters']:3d}  "
                  f"mu_edge={mm['mu_edge']:.4f}  "
                  f"mu_def={mm['mu_def']:.4f}  "
                  f"|M_net|={mm['Mnet']:.2f}{flag}")
        results[(shape, edge)] = row

    bal = max(results[("hexagon", "armchair")][u]["mu_edge"]
              for u in U_OVER_T_BAND)
    zig = min(results[("triangle", "zigzag")][u]["mu_edge"]
              for u in U_OVER_T_BAND)
    ok = (bal < 0.05) and (zig > 0.12) and (zig > 4 * max(bal, 1e-6))
    print("-" * 70)
    print(f"balanced max mu_edge = {bal:.4f}   "
          f"zigzag min mu_edge = {zig:.4f}")
    if ok:
        print("PASS -> Hubbard + corrected metric give the right physics "
              "across\nthe U/t band. The parasitic-moment axis is now real "
              "and the\ngeometry scan can finally begin.")
    else:
        print("FAIL -> ordering still not robust; diagnose before any scan "
              "(incremental discipline).")
    return {"balanced_mu_edge": bal, "zigzag_mu_edge": zig, "pass": bool(ok)}


# ============================================================================
# L4 - CORE SCAN  (the actual scientific deliverable)
# ----------------------------------------------------------------------------
# For every geometry, measure the TWO competing axes IN THE SAME RUN:
#
#   QUBIT  (maximise):  q = isolation * w_def
#       isolation = gap(dispersive) - |E_defect|   (defect level sits how
#                   deep inside the gap, away from band edges)
#       w_def     = weight of that level on the designed-defect region
#
#   PARASITIC (minimise):  p = mu_edge   (MF-Hubbard mean |moment| per rim
#                   atom; |dN| and |M_net| reported alongside as the
#                   Lieb-forced floor)
#
# Then extract the Pareto front in (q, p). The scientific question is
# answered by the SHAPE of that front: if some geometry is good on BOTH
# (high q, low p) the axes OVERLAP -> the literature's "armchair = free"
# story holds; if every high-q geometry also has high p, they TRADE OFF
# -> a genuine result that questions that story (the Vanevic competition,
# now quantified). This is exactly the gap vs Pedersen 2008/2013.
# ============================================================================

def evaluate_periodic(shape, edge, L=8, R=3.2, Nd=3, defect="single",
                      uot_band=U_OVER_T_BAND, nk=1):
    """PURE per-geometry evaluation -> combined qubit + parasitic metrics.

    nk=1 -> qubit axis at Gamma only (original behaviour, exact).
    nk>1 -> qubit axis stress-tested over an nk x nk k-mesh: isolation is
    taken as the WORST (min) over k, and the defect level's dispersion
    (max-min of E across k) is reported. A genuine localized qubit state is
    a near-flat band -> tiny dispersion; large dispersion means the Gamma
    'isolation' was an artifact (the level is really a band).
    """
    g = build_periodic_geometry(shape, edge, L, R, Nd=Nd, defect=defect)
    rec = dict(shape=shape, edge=edge, L=L, R=R, Nd=Nd, defect=defect, nk=nk,
               N=len(g["pos"]), dN=abs(g["dN"]),
               n_edge=int(g["edge_mask"].sum()),
               n_def=int(g["def_mask"].sum()))
    if rec["N"] < 16 or rec["n_def"] == 0:
        rec["status"] = "degenerate_geometry"
        return rec

    # ---- QUBIT axis: TB near-zero spectrum over a k-mesh ----------------
    ks = _kmesh(g["S1"], g["S2"], nk)
    iso_k, edef_k, wdef_k, gap_k = [], [], [], []
    for k in ks:
        vals, vecs = _near0(_bloch_H(g["pos"], g["S1"], g["S2"], k=k), k=24)
        prob = np.abs(vecs) ** 2
        prob /= prob.sum(0, keepdims=True)
        wdef = prob[g["def_mask"], :].sum(0)
        j = int(np.argmax(wdef))
        gp, nz = _robust_gap(vals)
        iso_k.append((gp - abs(vals[j])) if np.isfinite(gp)
                     else float("nan"))
        edef_k.append(float(vals[j]))
        wdef_k.append(float(wdef[j]))
        gap_k.append(gp)
        if k is ks[0]:
            ipr0 = float((prob[:, j] ** 2).sum())
            nz0 = nz
    iso_arr = np.array(iso_k, float)
    isolation = float(np.nanmin(iso_arr)) if np.isfinite(iso_arr).any() \
        else float("nan")                       # WORST-case over the BZ
    edef_disp = float(np.nanmax(edef_k) - np.nanmin(edef_k))  # flatness test
    rec.update(E_defect=float(edef_k[0]), w_def=float(wdef_k[0]),
               ipr_def=ipr0, gap=float(gap_k[0]), zero_modes=nz0,
               isolation=isolation, Edef_disp=edef_disp, nk_used=len(ks))

    # ---- PARASITIC axis: MF Hubbard over the U/t band (localized moment,
    #      Gamma is the standard for flat-band magnetism) ----------------
    mu_edge_band, mu_def_band, mnet_band = [], [], []
    for uot in uot_band:
        r = mean_field_hubbard(g["pos"], g["S1"], g["S2"], U=uot * T_HOP)
        mm = _moment_metrics(r["m"], g["edge_mask"], g["def_mask"])
        mu_edge_band.append(mm["mu_edge"])
        mu_def_band.append(mm["mu_def"])
        mnet_band.append(mm["Mnet"])
    # report worst-case (max) over the band -> robust, not a fragile U pick
    rec.update(mu_edge=float(np.max(mu_edge_band)),
               mu_def=float(np.max(mu_def_band)),
               Mnet=float(np.max(mnet_band)),
               mu_edge_band=[float(x) for x in mu_edge_band])

    # ---- combined scalars for the Pareto analysis -----------------------
    # A gapless geometry (zigzag flat-band: isolation = nan) is NOT invalid
    # data to drop - it is genuinely the WORST possible qubit. Coerce to
    # q = 0 so it stays in the comparison; silently excluding it biased the
    # verdict toward armchair (the earlier artifact).
    q = (rec["isolation"] * rec["w_def"]) \
        if np.isfinite(rec["isolation"]) else 0.0
    rec["q_qubit"] = float(q) if np.isfinite(q) else 0.0
    rec["p_parasitic"] = rec["mu_edge"]
    rec["status"] = "ok"
    return rec


def pareto_front(rows):
    """Indices of non-dominated geometries: maximise q_qubit, minimise
    p_parasitic. A point is dominated if another is >= in q AND <= in p
    (strictly better in at least one)."""
    ok = [r for r in rows if r.get("status") == "ok"
          and np.isfinite(r.get("q_qubit", np.nan))
          and np.isfinite(r.get("p_parasitic", np.nan))]
    front = []
    for i, a in enumerate(ok):
        dominated = False
        for b in ok:
            if b is a:
                continue
            if (b["q_qubit"] >= a["q_qubit"] and
                    b["p_parasitic"] <= a["p_parasitic"] and
                    (b["q_qubit"] > a["q_qubit"] or
                     b["p_parasitic"] < a["p_parasitic"])):
                dominated = True
                break
        if not dominated:
            front.append(a)
    return ok, front


def _edge_str(e):
    return e if isinstance(e, str) else f"chi={float(e):.2f}"


def core_scan(grid=None, uot_band=U_OVER_T_BAND, Nd=3, defect="single",
              n_jobs=-1):
    """Run the core scan over a geometry grid and read off the verdict.

    grid: list of (shape, edge, R) where edge is 'zigzag'/'armchair' OR a
    chirality fraction chi in [0,1]. Default now SPANS THE TRANSITION
    (hexagon, R sweep x chi sweep) - the two extremes were already known to
    separate (Vanevic); a trade-off, if it exists, lives in between.

    Parallel across CPU cores via joblib if available (each geometry is an
    independent pure call); falls back to serial. No GPU. The full default
    grid is ~tens of geometries x 3 U/t MF solves - minutes-to-~1h on a
    multicore CPU, which is why this parallelises.
    """
    if grid is None:
        # R is the REAL continuous knob: as it crosses successive lattice
        # shells the carved atom set genuinely changes. Orientation
        # (zigzag/armchair) is the discrete chirality contrast. We
        # OVER-generate then deduplicate by structural fingerprint, so the
        # scan reports the number of GENUINELY DISTINCT geometries instead
        # of silently double-counting rotated copies (the chi-freezing bug).
        Rs = [round(r, 3) for r in np.arange(2.4, 4.41, 0.15)]
        grid = [("hexagon", e, R) for R in Rs
                for e in ("zigzag", "armchair")]

    print("=" * 78)
    print(f"CORE SCAN  Nd={Nd} defect={defect}  U/t in {uot_band}")
    print("=" * 78)

    # ---- dedup pass: keep only structurally distinct geometries ---------
    import hashlib

    def _fingerprint(shape, edge, R):
        g = build_periodic_geometry(shape, edge, 8, R, Nd=Nd, defect=defect)
        if len(g["pos"]) == 0:
            return None, g
        P = g["pos"][np.lexsort((g["pos"][:, 1], g["pos"][:, 0]))]
        h = hashlib.sha1(np.round(P, 4).tobytes()).hexdigest()[:16]
        return (len(g["pos"]), abs(g["dN"]), h), g

    seen = {}
    uniq = []
    for (shape, edge, R) in grid:
        fp, _g = _fingerprint(shape, edge, R)
        if fp is None or fp in seen:
            continue
        seen[fp] = (shape, edge, R)
        uniq.append((shape, edge, R))
    n_dup = len(grid) - len(uniq)
    print(f"{len(grid)} candidates -> {len(uniq)} structurally DISTINCT "
          f"geometries ({n_dup} duplicates removed)")

    def _one(g):
        shape, edge, R = g
        return evaluate_periodic(shape, edge, R=R, Nd=Nd, defect=defect,
                                 uot_band=uot_band)

    rows = None
    try:
        from joblib import Parallel, delayed
        rows = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(_one)(g) for g in uniq)
        print(f"(ran in parallel via joblib, n_jobs={n_jobs})")
    except Exception as e:
        print(f"(serial fallback: {type(e).__name__}; "
              f"install joblib for parallel)")
        rows = [_one(g) for g in grid]

    print(f"{'shape':9s} {'edge':10s} {'R':>4s} {'N':>5s} {'|dN|':>4s} "
          f"{'gap':>6s} {'|Edef|':>7s} {'isol':>6s} {'w_def':>6s} "
          f"{'mu_edge':>8s} {'q':>7s} {'p':>7s}")
    print("-" * 78)
    for r in rows:
        if r.get("status") != "ok":
            print(f"{r.get('shape',''):9s} {_edge_str(r.get('edge','')):10s} "
                  f"-> {r.get('status')}")
            continue
        print(f"{r['shape']:9s} {_edge_str(r['edge']):10s} {r['R']:4.1f} "
              f"{r['N']:5d} {r['dN']:4d} {r['gap']:6.3f} "
              f"{abs(r['E_defect']):7.3f} {r['isolation']:6.3f} "
              f"{r['w_def']:6.3f} {r['mu_edge']:8.4f} "
              f"{r['q_qubit']:7.4f} {r['p_parasitic']:7.4f}")

    ok, front = pareto_front(rows)
    print("-" * 78)
    if not ok:
        print("No valid geometries - widen the grid.")
        return {"rows": rows, "front": []}

    print("Pareto-optimal geometries (non-dominated in q up / p down):")
    for r in sorted(front, key=lambda x: -x["q_qubit"]):
        print(f"  {r['shape']}/{_edge_str(r['edge'])} R={r['R']}  "
              f"q={r['q_qubit']:.4f}  p={r['p_parasitic']:.4f}")

    # verdict: do the axes overlap or trade off?
    best_q = max(ok, key=lambda r: r["q_qubit"])
    low_p = min(ok, key=lambda r: r["p_parasitic"])
    overlap = (id(best_q) == id(low_p))
    print("-" * 78)
    if overlap:
        print("VERDICT (preliminary): the highest-qubit geometry is ALSO "
              "the lowest-\nparasitic one -> axes OVERLAP on this grid. "
              "Check the front size:\nif the front is a single point the "
              "axes truly align; if it is a\ncurve, a trade-off still exists "
              "among near-optimal geometries.")
    else:
        print("VERDICT (preliminary): the highest-qubit geometry is NOT the "
              "lowest-\nparasitic one -> the axes TRADE OFF (Pareto tension). "
              "This is the\nscientifically valuable outcome - the Vanevic "
              "competition, quantified,\nand the point of departure from "
              "Pedersen 2008/2013.")
    print(f"Pareto front has {len(front)} point(s) of {len(ok)} valid "
          f"geometries.")
    print("NOTE: preliminary - single defect, Gamma point, one Nd. The "
          "verdict is\na direction, not proof; robustness (Nd, k-mesh, "
          "disorder) comes next.")
    return {"rows": rows, "front": front}


# ============================================================================
# L5 - ROBUSTNESS TEST  (the honest go / no-go)
# ----------------------------------------------------------------------------
# The core scan's "trade-off" rests on one correlation seen at Gamma, one Nd:
# larger sublattice imbalance |dN| -> larger qubit isolation q AND larger
# parasitic moment p. This layer asks ONE pre-registered question, decided
# BEFORE looking at the numbers:
#
#   Does the ORDER  (|dN| up  =>  q up  AND  p up)  survive
#       (a) a k-mesh instead of Gamma, and
#       (b) a change of supercell size Nd ?
#
# Verdict rule (fixed now, not tunable to the data):
#   * ROBUST  : the sign of (q_hi - q_lo) and (p_hi - p_lo) is positive in
#               EVERY tested condition (Gamma & k-mesh, each Nd).
#   * ARTIFACT: the sign flips in ANY condition.
# Plus an independent honesty flag: if the defect level disperses
# (Edef_disp > 0.10 * gap) the "qubit" is really a band -> the q axis
# itself was a Gamma artifact for that geometry, reported explicitly.
# ============================================================================

def robustness_test(reps=None, nk=2, Nds=(3, 4), uot_band=(1.15,),
                     n_jobs=-1):
    """reps: list of (shape, edge, R). Default auto-picks two structurally
    distinct geometries - the lowest-|dN| and highest-|dN| ones - from a
    small candidate set, i.e. the two carrying the correlation under test.
    """
    if reps is None:
        import hashlib
        cand = [("hexagon", e, round(r, 3))
                for r in np.arange(2.4, 3.21, 0.2)
                for e in ("zigzag", "armchair")]
        seen, info = {}, []
        for (s, e, R) in cand:
            g = build_periodic_geometry(s, e, 8, R, Nd=3, defect="single")
            if len(g["pos"]) == 0:
                continue
            P = g["pos"][np.lexsort((g["pos"][:, 1], g["pos"][:, 0]))]
            h = hashlib.sha1(np.round(P, 4).tobytes()).hexdigest()[:12]
            key = (len(g["pos"]), abs(g["dN"]), h)
            if key in seen:
                continue
            seen[key] = 1
            info.append((abs(g["dN"]), s, e, R))
        info.sort()
        lo = info[0]
        hi = info[-1]
        reps = [(lo[1], lo[2], lo[3]), (hi[1], hi[2], hi[3])]
        print(f"auto reps: LO |dN|~{lo[0]} {reps[0]}   "
              f"HI |dN|~{hi[0]} {reps[1]}")

    print("=" * 78)
    print(f"ROBUSTNESS TEST  reps={reps}  nk={nk}  Nds={Nds}  "
          f"U/t in {uot_band}")
    print("=" * 78)

    conditions = []
    for Nd in Nds:
        for use_nk in (1, nk):
            conditions.append((Nd, use_nk))

    def _ev(args):
        (s, e, R), Nd, knk = args
        return (Nd, knk, evaluate_periodic(s, e, R=R, Nd=Nd,
                defect="single", uot_band=uot_band, nk=knk))

    jobs = [((s, e, R), Nd, knk) for (s, e, R) in reps
            for (Nd, knk) in conditions]
    try:
        from joblib import Parallel, delayed
        out = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(_ev)(j) for j in jobs)
    except Exception as e:
        print(f"(serial fallback: {type(e).__name__})")
        out = [_ev(j) for j in jobs]

    # index results: results[(Nd,knk)] = [rec_lo, rec_hi] in reps order
    res = {}
    for (s, e, R), idx in zip(reps, range(len(reps))):
        pass
    table = {}
    for (Nd, knk, rec) in out:
        table.setdefault((Nd, knk), {})
        # which rep? match by (shape,edge,R)
        keyrep = (rec["shape"], rec["edge"], rec["R"])
        table[(Nd, knk)][keyrep] = rec

    print(f"{'Nd':>3s} {'k':>4s} {'rep':<22s} {'|dN|':>4s} {'gap':>6s} "
          f"{'isol':>6s} {'Edisp':>6s} {'w_def':>6s} {'mu_edge':>8s} "
          f"{'q':>7s} {'p':>7s} flat?")
    print("-" * 78)
    verdict_ok = True
    flips = []
    disp_flags = []
    for (Nd, knk) in conditions:
        recs = [table[(Nd, knk)][(s, e, R)] for (s, e, R) in reps]
        for rec in recs:
            flat = (np.isfinite(rec.get("gap", np.nan)) and rec["gap"] > 0
                    and rec.get("Edef_disp", 9) <= 0.10 * rec["gap"])
            if not flat:
                disp_flags.append((Nd, knk, rec["shape"], rec["edge"],
                                   rec["R"]))
            print(f"{Nd:3d} {('G' if knk==1 else f'{knk}x{knk}'):>4s} "
                  f"{rec['shape']+'/'+_edge_str(rec['edge'])+' R'+str(rec['R']):<22s} "
                  f"{rec['dN']:4d} {rec.get('gap',float('nan')):6.3f} "
                  f"{rec.get('isolation',float('nan')):6.3f} "
                  f"{rec.get('Edef_disp',float('nan')):6.3f} "
                  f"{rec.get('w_def',float('nan')):6.3f} "
                  f"{rec.get('mu_edge',float('nan')):8.4f} "
                  f"{rec.get('q_qubit',float('nan')):7.4f} "
                  f"{rec.get('p_parasitic',float('nan')):7.4f} "
                  f"{'yes' if flat else 'NO'}")
        lo, hi = recs[0], recs[1]
        dq = hi.get("q_qubit", np.nan) - lo.get("q_qubit", np.nan)
        dp = hi.get("p_parasitic", np.nan) - lo.get("p_parasitic", np.nan)
        ok_c = np.isfinite(dq) and np.isfinite(dp) and dq > 0 and dp > 0
        if not ok_c:
            verdict_ok = False
            flips.append((Nd, knk, round(float(dq), 4),
                          round(float(dp), 4)))
        print(f"      -> dq(hi-lo)={dq:+.4f}  dp(hi-lo)={dp:+.4f}  "
              f"order {'HOLDS' if ok_c else 'BROKEN'}")
        print("-" * 78)

    print("PRE-REGISTERED VERDICT:")
    if verdict_ok and not disp_flags:
        print("  ROBUST -> the |dN|->(q,p) correlation survives k-mesh AND "
              "Nd change,\n  and the defect level stays flat. The trade-off "
              "is a real direction,\n  not a Gamma/finite-cell artifact. A "
              "small but defensible result.")
    elif verdict_ok and disp_flags:
        print("  MIXED -> ordering survives, BUT the defect level disperses "
              "in some\n  conditions (not a flat qubit band there). The q "
              "axis is partly a\n  Gamma artifact; the result is weaker than "
              "it looked.")
        for d in disp_flags:
            print(f"    dispersive: Nd={d[0]} k={d[1]} {d[2]}/{d[3]} R={d[4]}")
    else:
        print("  ARTIFACT -> the ordering flips under k-mesh and/or Nd "
              "change.\n  The core-scan 'trade-off' was a finite-cell / "
              "Gamma artifact. The\n  honest outcome: this method at this "
              "scale does NOT establish the\n  effect. A clean negative - "
              "still a publishable methodological note.")
        for f in flips:
            print(f"    flip at Nd={f[0]} k={f[1]}: dq={f[2]} dp={f[3]}")
    return {"table": table, "robust": bool(verdict_ok and not disp_flags),
            "flips": flips, "dispersive": disp_flags}


def periodic_validation(shape="hexagon", edge="armchair", L=8, R=3.2, Nd=3):
    """The two MANDATORY checks before any metric scan (incremental
    discipline): (1) the periodic machinery removes outer-edge artifacts
    -> pristine cell shows a clean TB gap; (2) omitting one antidot puts a
    discrete state INSIDE that gap, localised on the defect -> the qubit
    metric is resolvable (it was NOT on the finite flake).

    Baseline is hexagon+armchair: a triangular hole is intrinsically
    sublattice-imbalanced (3-fold symmetry vs bipartite balance), so it is a
    poor 'clean, moment-free qubit' baseline. A hexagonal armchair hole can
    be balanced -> clean gap and a single clean defect level. Triangular
    zigzag is brought in LATER, for the parasitic-moment axis, where the
    imbalance is the wanted signal.
    """
    print("=" * 70)
    print(f"PERIODIC VALIDATION  shape={shape} edge={edge} "
          f"L={L} R={R} Nd={Nd}  (TB-only)")
    print("=" * 70)

    # ---- check 1: pristine cell, clean gap (machinery sanity) -------------
    g1 = _pristine_cell(shape, edge, L, R)
    vals, _ = _near0(_bloch_H(g1["pos"], g1["S1"], g1["S2"],
                              k=np.zeros(2)), k=16)
    gap, nz = _robust_gap(vals)
    ok1 = np.isfinite(gap) and gap > 0.05
    print(f"[1] pristine  N={len(g1['pos'])}  |dN|={abs(g1['dN'])}  "
          f"zero_modes={nz}  TB gap (dispersive bands) ~ {gap:.3f} eV   "
          f"{'OK clean gap' if ok1 else 'NO gap (machinery suspect)'}")

    # ---- check 2: defect supercell, state inside the gap ------------------
    g2 = build_periodic_geometry(shape, edge, L, R, Nd=Nd, defect="single")
    vals, vecs = _near0(_bloch_H(g2["pos"], g2["S1"], g2["S2"],
                                 k=np.zeros(2)), k=24)
    prob = np.abs(vecs) ** 2
    prob /= prob.sum(0, keepdims=True)
    wdef = prob[g2["def_mask"], :].sum(0)
    wedg = prob[g2["edge_mask"], :].sum(0)
    j = int(np.argmax(wdef))
    ok2 = wdef[j] > 0.30
    print(f"[2] defected  N={len(g2['pos'])}  |dN|={abs(g2['dN'])}")
    print(f"    spectrum near 0 (eV): "
          f"{np.round(np.sort(vals)[:8], 3).tolist()}")
    print(f"    defect state: E={vals[j]:+.3f}  w_defect={wdef[j]:.3f}  "
          f"w_edge={wedg[j]:.3f}  "
          f"{'OK qubit state resolved' if ok2 else 'still buried'}")
    print("-" * 70)
    if ok1 and ok2:
        print("BOTH CHECKS PASS -> periodic machinery validated; the "
              "next layer\n(mean-field Hubbard for the parasitic-moment "
              "axis) can be added on top.")
    else:
        print("A check failed -> stop here and diagnose before scanning "
              "(this is the incremental discipline working as intended).")
    return {"gap": gap, "ok_pristine": bool(ok1),
            "w_defect": float(wdef[j]), "ok_defect": bool(ok2)}


def _pristine_cell(shape, edge, L, R):
    """One antidot cell with its hole carved (NO omission) for the gap test."""
    pos, sub, S1, S2 = _supercell_atoms(L, Nd=1)
    c = 0.5 * (S1 + S2)
    pred = shape_predicate(shape, edge)
    dm = _min_image(pos - c, S1, S2)
    pos, sub = pos[~pred(dm + c, c, R)], sub[~pred(dm + c, c, R)]
    for _ in range(12):
        pairs, _o = _pbc_pairs(pos, S1, S2)
        deg = np.zeros(len(pos), np.int64)
        if len(pairs):
            np.add.at(deg, pairs[:, 0], 1)
            np.add.at(deg, pairs[:, 1], 1)
        keep = deg >= 2
        if keep.all():
            break
        pos, sub = pos[keep], sub[keep]
    dN = int(np.sum(sub == 0) - np.sum(sub == 1))
    return dict(pos=pos, sub=sub, S1=S1, S2=S2, dN=dN)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shape", default="triangle",
                    choices=["triangle", "hexagon", "circle"])
    ap.add_argument("--lattice", default="triangular",
                    choices=["triangular", "square"])
    ap.add_argument("--edge", default="zigzag",
                    choices=["zigzag", "armchair"])
    ap.add_argument("--defect", default="single",
                    choices=["single", "double"])
    ap.add_argument("--demo", action="store_true",
                    help="run the finite-flake demo (legacy/reference)")
    ap.add_argument("--periodic", action="store_true",
                    help="run the periodic-supercell validation (current path)")
    args, _unknown = ap.parse_known_args()   # ignore stray Jupyter kernel args

    if args.periodic:
        periodic_validation()
    elif args.demo:
        _demo()
    else:
        p = _base_params()
        p.update(shape=args.shape, lattice=args.lattice,
                 edge=args.edge, defect=args.defect)
        from pprint import pprint
        pprint(evaluate_geometry(p))
