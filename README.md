# Holey-graphene (antidot-lattice) defect qubits vs. parasitic edge magnetism: a pre-registered negative result

**Status:** Manuscript in preparation (unpublished computational study). **Outcome: negative (ARTIFACT by the pre-registered criterion).**

Author: Faruk Sami Bilgin, Independent Researcher (ORCID [0009-0001-9064-4410](https://orcid.org/0009-0001-9064-4410))

## Question

Leaving one antidot uncarved in a graphene antidot lattice creates a localized in-gap "designed defect" state proposed as a spin-qubit host (Pedersen et al., PRL 100, 136804, 2008). Separately, local magnetic moments on antidot edges were identified as a competing dephasing source (Vanević et al., 2009). These two observations had not been confronted in the same geometry scan: across a family of antidot geometries, do qubit quality and parasitic edge moment **overlap**, or do they **trade off** — and does the answer survive changes of scale?

## Method

Single-file pipeline `scripts/holey_graphene_qubit.py` (NumPy + SciPy + joblib, CPU only):
- **Geometry:** periodic, integer-commensurate supercells (Pedersen construction): antidot cell of L×L graphene cells; defect supercell of Nd×Nd antidot cells with the central hole omitted; hexagonal/triangular/circular hole shapes with zigzag or armchair edges. a_cc = 1, nearest-neighbour hopping t = 2.7 eV. (An earlier finite-flake path, kept as `--demo`, was abandoned because the flake boundary swamped the defect state.)
- **Qubit axis:** sparse Bloch Hamiltonian, near-Fermi states by shift-invert ARPACK with a fixed start vector (determinism); quality score **q = isolation × w_def**, isolation = gap − |E_defect|.
- **Parasitic axis:** collinear π-orbital **mean-field Hubbard** at half filling; **p = μ_edge** = mean |m_i| on antidot-rim atoms.
- **Duplicate removal:** geometries are de-duplicated by a structural fingerprint (small rotations of a discrete lattice remove the same atoms).
- **Parallel scan:** each geometry is an independent pure function mapped with joblib.

### Pre-registered decision criterion (fixed before the data were seen)
Does the ordering (|dN| ↑ ⇒ q ↑ AND p ↑) survive (a) a 2×2 k-mesh instead of Γ and (b) a change of supercell size Nd = 3 → 4?
- **ROBUST:** the signs of (q_hi − q_lo) and (p_hi − p_lo) stay positive in **every** tested condition.
- **ARTIFACT:** the sign flips in **any** condition.
- Independent flag: if the defect level disperses (E_disp > 0.10 × gap), the "qubit" is a band in that geometry.

## Results

Representatives chosen automatically: LO = hexagon/armchair, R = 2.8 (|dN| ≈ 8); HI = hexagon/zigzag, R = 2.8 (|dN| ≈ 16); U/t = 1.15. Full table: [`results/robustness_test_table.csv`](results/robustness_test_table.csv).

| Nd | k | geometry | \|dN\| | gap (eV) | isolation | E_disp | w_def | μ_edge | q | p | flat? |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | Γ | hex/armchair | 8 | 0.948 | 0.208 | 0.000 | 0.358 | 0.0298 | 0.0746 | 0.0298 | yes |
| 3 | Γ | hex/zigzag | 16 | 1.215 | 0.607 | 0.000 | 0.253 | 0.0962 | 0.1534 | 0.0962 | yes |
| 3 | 2×2 | hex/armchair | 8 | 0.941 | 0.187 | 1.466 | 0.244 | 0.0298 | 0.0457 | 0.0298 | NO |
| 3 | 2×2 | hex/zigzag | 16 | 1.211 | 0.590 | 0.000 | 0.277 | 0.0962 | 0.1633 | 0.0962 | yes |
| 4 | Γ | hex/armchair | 15 | 0.960 | 0.480 | 0.000 | 0.190 | 0.0307 | 0.0912 | 0.0307 | yes |
| 4 | Γ | hex/zigzag | 30 | gapless | – | 0.000 | 0.034 | 0.0967 | 0.0000 | 0.0967 | NO |
| 4 | 2×2 | hex/armchair | 15 | 0.964 | 0.482 | 0.964 | 0.196 | 0.0307 | 0.0946 | 0.0307 | NO |
| 4 | 2×2 | hex/zigzag | 30 | gapless | – | 0.000 | 0.030 | 0.0967 | 0.0000 | 0.0967 | NO |

- **Nd = 3:** ordering holds — dq(hi−lo) = +0.079 (Γ), +0.118 (2×2); dp = +0.066.
- **Nd = 4:** ordering breaks — dq = −0.091 (Γ), −0.095 (2×2); dp = +0.066. The zigzag geometry is gapless (q forced to 0) and the armchair defect level disperses (E_disp = 0.964).
- The dispersion flag fired in four of the eight conditions.

**Pre-registered verdict: ARTIFACT.** The apparent parasitic-moment / qubit-quality trade-off seen at Γ in the small supercell is a by-product of the finite cell / Γ-point choice; this method at this scale does not establish the effect.

### Independent re-run (2026-11-15) — same verdict, different details
The test was re-run for this repository (`python scripts/run_robustness_test.py`, Python 3.13.3, numpy 2.2.5, scipy 1.18.0, joblib 1.5.3; wall time 398.7 s). Output: [`results/robustness_test.json`](results/robustness_test.json), [`results/robustness_test_rerun_2026-11-15.log`](results/robustness_test_rerun_2026-11-15.log).

- The pre-registered verdict is again **ARTIFACT**.
- But several entries differ from the original run: at Nd = 4, Γ the zigzag geometry is now gapped (gap 1.217, q = 0.1396) and the ordering **holds** there (dq = +0.0499); the only flip is at Nd = 4, 2×2 (dq = −0.0931, dp = +0.066). At Nd = 3 the signs agree with the report (dq = +0.0718 Γ, +0.1130 2×2; dp = +0.0663), but w_def, q and the dispersion flags differ (e.g. the Nd = 3, 2×2 dispersion flag moves from armchair to zigzag).
- Interpretation: some near-zero-mode quantities are not numerically reproducible across library versions (shift-invert ARPACK near exact zero modes). The verdict does not depend on this, since it requires only one sign flip, but individual table entries should be read as version-dependent. This reproducibility issue is unresolved.

### What this does and does not say
It does not say the architecture is closed to qubits. It says the Nd = 3 trend (larger sublattice imbalance feeding both the defect state and the parasitic moment, consistent in direction with Vanević's concern) is not scale-robust. Limits: only ~4–6 structurally distinct geometries exist at this scale; single defect only; mean-field Hubbard; only Nd = 3 and 4 (~1700–2000 atoms at Nd = 4, near the desktop limit of this approach).

## How the results were verified
Eight validation layers, each with a pass/fail check on the smallest exact system, and each caught error reported (finite-flake boundary contamination → periodic supercell; seed-dependent ARPACK output → fixed v0, bit-identical over 3 runs; a wrong Σ|m_i| moment metric → per-atom μ_edge and M_net; shape collapse at small R → equal-area normalisation; silently dropped gapless geometries → q = 0; rotated duplicates → structural fingerprint). The final verdict is the mechanical application of the rule stated before the run. Full record (Turkish): [`docs/tr/holey_graphene_qubit_negative_result.md`](docs/tr/holey_graphene_qubit_negative_result.md).

## Reproducing
```bash
pip install -r requirements.txt
python scripts/holey_graphene_qubit.py --periodic     # periodic-machinery validation
python scripts/run_robustness_test.py                 # pre-registered test -> results/robustness_test.json
```
The original test was run on an Intel i7-12700H / 16 GB desktop with Python 3.9.25, numpy 1.26.4, scipy 1.13.1, joblib 1.5.3 (per the report). The main table above is transcribed from the report; the re-run is described above.

## Related work by the author
- Twisted hBN preprint: https://doi.org/10.5281/zenodo.24201814
- Twisted MoS2 preprint: https://www.alphaxiv.org/abs/2611.twisted-mos2-hole-spin-coherence
- Moiré microwave loss: https://doi.org/10.5281/zenodo.24204124
- josim-cpr-bridge: https://github.com/FarukSamiBilgin/josim-cpr-bridge

## Citation
```bibtex
@misc{bilgin_holey_graphene_qubit,
  author = {Bilgin, Faruk Sami},
  title  = {A pre-registered negative result on the parasitic-edge-moment versus qubit-quality trade-off in holey-graphene defect qubits},
  year   = {2026},
  note   = {Manuscript in preparation},
  url    = {https://github.com/FarukSamiBilgin/holey-graphene-antidot-qubit}
}
```

AI assistance: AI coding assistants were used for code development and for adversarial review. All load-bearing results were re-derived by independent methods and checked by automated tests; the author is responsible for every claim.

## License
Code: MIT (`LICENSE`). Text, figures and data: CC BY 4.0 (`LICENSE-CONTENT`).
