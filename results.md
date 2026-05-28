# Results

Verified on 2026-05-28 by running the scripts from the repository root as `python3 scripts/<name>.py`.

## Verification summary

- `python3 -m py_compile scripts/*.py` passed.
- Every file under `scripts/` now executes cleanly from the repo root.
- Run logs were saved under `outputs/run_logs/`.

## Fixes made during verification

- The plotting scripts were importing `scripts.plot_io`, which fails under direct execution (`python3 scripts/foo.py`). They now fall back to `from plot_io import save_figure` when needed.
- `scripts/bonferroni_failed_readressing.py` was previously treating undefined jump-mode tests as if they had `p=0`. It now propagates `NaN` correctly and ranks finite results only.

## Script-by-script status

### `scripts/plot_io.py`

- Exit status: success.
- Notes: helper module only; no output by design.

### `scripts/5_copula_test.py`

- Exit status: success.
- Output: numeric summary only.
- Strongest tail-dependence signal in the tested set was `(p, q) = (3, 5)` with observed tail dependence `0.154867` versus null mean `0.109991`, permutation `p=0.0160`.

### `scripts/bonferroni_failed_readressing.py`

- Exit status: success.
- Output: numeric summary only.
- Level-mode scan produced several strong signals after BH correction; the strongest printed rows include:
  - `(2, 5), N=500, alpha=0.85`, observed `0.4167`, null `0.1618`, `p=0.0000`, `BH-q=0.0000`
  - `(5, 11), N=500, alpha=0.85`, observed `0.1750`, null `0.0801`, `p=0.0000`, `BH-q=0.0000`
  - `(2, 3), N=5000, alpha=0.90`, observed `0.1411`, null `0.0876`, `p=0.0000`, `BH-q=0.0000`
- Jump / carry-change mode no longer reports bogus zero p-values. The best printed jump-mode row was `(5, 11), N=500, alpha=0.90` with `p=0.2660`, so there is no compelling jump-mode signal in the top rows.

### `scripts/borwein_integral_logistic_6.py`

- Exit status: success.
- Output: numeric summary only.
- `rcrit = 1 + sqrt(6) = 3.44948974278317809819728407471`
- Logistic Lyapunov at `rcrit` came out near zero: `-0.00012013927484664507908`.
- Largest c10 jump across the tower ladder was `L3 -> L4 = 0.039597`.

### `scripts/prime_index_dependence_test.py`

- Exit status: success.
- Output: numeric summary only.
- Prime-indexed subsequence test used `N=50000`, `warmup>100`, `5108` prime indices, `B=300`.
- Observed mean off-diagonal cross-base correlation was `0.116390`.
- Uniform control mean was `0.092338 +/- 0.005663`, Monte Carlo `p=0.0033`.
- Residue-matched control mean was `0.094929 +/- 0.004369`, Monte Carlo `p=0.0033`.
- Gap-linked effects remained visible after residue matching, especially for bases `2`, `5`, and `11` with matched-control `q=0.0042`.

### `scripts/pi_e_copulas.py`

- Exit status: success.
- Output: eight copula figures under `outputs/prime_copulas/`.
- Notable observed Spearman correlations:
  - `(2, 5)`: `0.223641`
  - `(61, 1103)`: `0.325592`
  - `(7, 1103)`: `0.364694`

### `scripts/cf_recurrence.py`

- Exit status: success.
- Output files:
  - `recurrence.png`
  - `outputs/recurrence/sqrt2.png`
  - `outputs/recurrence/e.png`
  - `outputs/recurrence/pi.png`

### `scripts/spinor_cocycles.py`

- Exit status: success.
- Output files saved under `outputs/spinor_cocycles/`.
- Summary values:
  - Mean cocycle collision: `0.7996125336724942`
  - Mean spinor mismatch Newton/Collatz: `0.6168496654521637`
  - Mean spinor mismatch Newton/Pisot: `0.5688242172532343`
  - Mean spinor mismatch Collatz/Pisot: `0.4863814833564132`

### `scripts/spinor_split_plot.py`

- Exit status: success.
- Output files:
  - `spinor_split.png`
  - `outputs/spinor_split/growth.png`
  - `outputs/spinor_split/scatter.png`
  - `outputs/spinor_split/histogram.png`

### `scripts/tower_maps.py`

- Exit status: success.
- Output files:
  - `gt_tower.png`
  - `pt_tower.png`
  - `pisot_tower.png`
  - split panels under `outputs/gt_tower/`, `outputs/pt_tower/`, and `outputs/pisot_tower/`

## Main generated artifacts

Top-level combined figures:

- `gt_tower.png`
- `pt_tower.png`
- `pisot_tower.png`
- `recurrence.png`
- `spinor_split.png`

Split figures and logs:

- `outputs/prime_copulas/`
- `outputs/gt_tower/`
- `outputs/pt_tower/`
- `outputs/pisot_tower/`
- `outputs/recurrence/`
- `outputs/spinor_cocycles/`
- `outputs/spinor_split/`
- `outputs/run_logs/`
