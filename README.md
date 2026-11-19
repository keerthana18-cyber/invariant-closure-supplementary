# Invariant closure theory: supplementary audit

Supplementary code for the manuscript "A Singularity-Driven Invariant Closure Theory for
Traveling Wave Solutions of Nonlinear Partial Differential Equations".

## Contents

| Path | Purpose |
|---|---|
| `scripts/search.py` | Closure ring: reduction modulo `Y^2 - Q(z)` and the derivation `D` |
| `scripts/run.py` | Affine normalizations of `Q` used in the supplementary audit |
| `scripts/audit612.py` | Main exact-degree audit, 612 cases for KdV, mKdV, Burgers, Boussinesq and Gardner (Table "audit") |
| `results/results612.jsonl` | One JSON record per case of the main audit: case, status, certificate type, seconds |
| `scripts/box.py` | Family II exact-degree audit for Fisher and KdV-Burgers (72 cases, Table "Supplementary audit") |
| `scripts/boxI.py` | Family I audit for Fisher and KdV-Burgers (8 cases) |
| `scripts/verify_new.py` | Exact symbolic verification of the Fisher front, the KdV-Burgers shock and the Gardner kink |
| `results/box_results.json` | Machine-readable status of every Family II case |
| `figures/make_profiles.py` | Regenerates `profiles.png` (Figure "profiles", 600 dpi) and `profiles.pdf` |

## Running

```bash
pip install -r requirements.txt
python scripts/verify_new.py          # seconds; every residual prints 0
python scripts/audit612.py 600 fresh  # about one minute; without 'fresh' it resumes from results/results612.jsonl
python scripts/boxI.py                # a few minutes
python scripts/box.py                 # about one hour on two cores
python figures/make_profiles.py
```

The audit runs in search class (B) of the manuscript: no Stage 1 (singularity) information is imposed, so every certificate holds without any assumption on the poles of the auxiliary function.

A case is recorded as `NO_SOLUTION_CERTIFIED` when the reduced Groebner basis of the
coefficient system, with the exact-degree and nonzero-dissipation conditions imposed through a
Rabinowitsch variable, equals `{1}`.

## Environment

Python 3.13, SymPy 1.14.0, mpmath 1.3.0.

## Independent verification (`verification/`)

These scripts check the audit without using `scripts/search.py` or SymPy's Groebner bases.

| File | Check |
|---|---|
| `indep.py` | Independent derivation of the coefficient systems (chain rule on symbols, polynomial division by `y^2 - Q`) |
| `singular.py` | Interface to Singular 4.3.2; clears denominators and aborts on any Singular error |
| `verify612.py norm` | All 612 cases in Singular; for every inconsistent case, `lift` returns cofactors `c_i` with `sum c_i f_i = 1`, verified by exact expansion. Results: `verify_norm.jsonl` |
| `verify612.py gen` | The same without normalizing `Q` (all coefficients unknown). Results: `verify_gen.jsonl` |
| `verify_supp.py` | Singular cofactor certificates for the 80-case supplementary audit (Fisher, KdV-Burgers). Results: `verify_supp.jsonl` |
| `gardner_2_2_4_3_certificate.py` | Explicit certificate for the one case where Singular's `lift` returned no cofactors |
| `squares.sing` | Ideal of the square locus for the normalized quadratic and quartic, by elimination |
| `normcheck.py` | The affine normalization commutes with the derivation, for every `m <= 5` and `k = 1, 2` |
| `verify_elliptic.py` | Numerical check (30 digits) of the KdV cnoidal and Gardner elliptic closed forms |
| `found13.py` | The 13 found closures: exact residual with `indep.py`, and numerical integration of the auxiliary equation with evaluation of the original reduced ODE. Results: `found13_check.json` |

### Status labels in the verification results

| File | Label | Meaning |
|---|---|---|
| `verify_norm.jsonl`, `verify_gen.jsonl` | `CERT_INCONSISTENT` | Singular returned cofactors `c_i` with `sum c_i f_i = 1`, confirmed by exact expansion: no closure exists |
| `verify_supp.jsonl` | `CERT_VERIFIED` | As above, for the supplementary audit |
| all | `EXISTS`, `EXISTS_NONSQUARE` | The system is consistent: a closure exists (these are the found closures) |
| `verify_norm.jsonl` | `CERT_FAILED_CHECK` | Only the case Gardner, `(k,m)=(2,2)`, `(d_A,d_B)=(4,3)`: Singular's `lift` returned no cofactors. The case is settled by the explicit certificate in `gardner_2_2_4_3_certificate.py` |
| `verify_gen.jsonl` | `UNDECIDED` | Only the same Gardner case, which exceeded the 30-minute limit without normalization; settled by the same certificate |

Requires Singular (`apt-get install singular`).

```bash
python verification/verify612.py norm fresh   # about 25 minutes; without 'fresh' it resumes
python verification/verify_supp.py            # under one minute
python verification/found13.py                # about one minute
python verification/gardner_2_2_4_3_certificate.py
python verification/normcheck.py
Singular -q verification/squares.sing
```

All scripts locate their inputs and write their outputs relative to their own location,
so they can be run from any working directory.

## License

The code and the result files are released under the MIT License (see `LICENSE`).
