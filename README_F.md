# SMERFS Fortran Interface

Fortran 90 interface and example for the SMERFS C library.

## Overview

This adds two things on top of the existing C kernels in `src/`:

| File | Purpose |
|---|---|
| `src/smerfs_interface.f90` | Fortran 90 module wrapping all public C functions via `iso_c_binding` |
| `src/CMakeLists.txt` | Builds `libsmerfs_c.a` (C kernels) and `libsmerfs_f.a` (Fortran interface) |
| `examples_f/example1_f.F90` | Fortran translation of `examples/example1.py` |
| `examples_f/CMakeLists.txt` | Builds the `example1_f` executable |

## Build environment

The build requires the **nag-stack** Lmod environment module, which sets up
the NAG Fortran compiler (`nagfor`) and Apple clang as the C compiler:

```bash
ml nag-stack
```

## Build instructions

### 1. Build the libraries (`src/`)

```bash
ml nag-stack
mkdir -p build/src
cmake src/ -B build/src -DCMAKE_BUILD_TYPE=Release
cmake --build build/src
```

Artifacts produced:

| File | Contents |
|---|---|
| `build/src/libsmerfs_c.a` | C kernels (cov, linalg, smerfs, ziggurat) |
| `build/src/libsmerfs_f.a` | Fortran interface objects |
| `build/src/modules/smerfs_interface.mod` | NAG `.mod` file for `USE smerfs_interface` |

### 2. Build the Fortran example (`examples_f/`)

```bash
cmake examples_f/ -B build/examples_f \
      -DCMAKE_BUILD_TYPE=Release \
      -DSMERFS_SRC_BUILD=$(pwd)/build/src
cmake --build build/examples_f
```

### 3. Run the example

```bash
cd build/examples_f
./example1_f
```

Outputs:

| File | Contents |
|---|---|
| `realisation.csv` | 128 × 256 grid — one realisation of the GRF |
| `correl.csv` | 1001 rows (`z, C(z)`) — analytic covariance function |

## Fortran interface module

`src/smerfs_interface.f90` exposes seven public wrapper subroutines/functions.
All scalar integers and reals from C are mapped via `iso_c_binding` kinds;
C `double complex` maps to Fortran `complex(c_double_complex)`.

### `update_cov_f`

```fortran
subroutine update_cov_f(m_max, N, M, norm_re, norm_im, llp1_re, llp1_im, &
                         F, H, tau_power, eta_ratio2, cov, cross_cov, rc)
```

Accumulates one partial-fraction-decomposition term into the covariance and
cross-covariance matrices.  Wraps `update_cov` in `cov.c`.

| Argument | Type | Intent | Description |
|---|---|---|---|
| `m_max` | `integer` | in | Final m index (0 .. m_max) |
| `N` | `integer` | in | Number of z points |
| `M` | `integer` | in | Markov order |
| `norm_re/im` | `real(c_double)` | in | Real/imag parts of the normalisation factor |
| `llp1_re/im` | `real(c_double)` | in | Real/imag parts of l(l+1) |
| `F` | `complex(c_double_complex)(*)` | in | (m_max+M+1)×N hypergeometric table at x=(1-z)/2 |
| `H` | `complex(c_double_complex)(*)` | in | (m_max+M+1)×N hypergeometric table at y=(1+z)/2 |
| `tau_power` | `real(c_double)(*)` | in | (2M-1)×N twiddle factors τ^p |
| `eta_ratio2` | `real(c_double)(*)` | in | (N-1) ratios τ_{i+1}/τ_i |
| `cov` | `real(c_double)(*)` | inout | (m_max+1)×N×M×M covariance accumulator |
| `cross_cov` | `real(c_double)(*)` | inout | (m_max+1)×(N-1)×M×M cross-covariance accumulator |
| `rc` | `integer` | out | 0 = success, −1 = out-of-memory |

### `inverse_f`

```fortran
subroutine inverse_f(N, M, matrices, out, rc)
```

Inverse of N symmetric M×M matrices (M = 2 or 3).  Wraps `inverse` in `linalg.c`.

| Argument | Type | Intent | Description |
|---|---|---|---|
| `N` | `integer` | in | Number of matrices |
| `M` | `integer` | in | Matrix dimension (2 or 3) |
| `matrices` | `real(c_double)(N*M*M)` | in | Input matrices |
| `out` | `real(c_double)(N*M*M)` | out | Inverses |
| `rc` | `integer` | out | 0 = success; i+1 = matrix i singular; −1 = M unsupported |

### `cholesky_f`

```fortran
subroutine cholesky_f(N, M, matrices, out, rc)
```

Batched lower Cholesky factorisation of N symmetric M×M matrices (M = 2).
Wraps `cholesky` in `linalg.c`.

| Argument | Type | Intent | Description |
|---|---|---|---|
| `N` | `integer` | in | Number of matrices |
| `M` | `integer` | in | Matrix dimension (currently 2 only) |
| `matrices` | `real(c_double)(N*M*M)` | in | Input matrices |
| `out` | `real(c_double)(N*M*M)` | out | Lower triangular Cholesky factors |
| `rc` | `integer` | out | 0 = success; i+1 = matrix i not pos-def; −1 = M unsupported |

### `state_space_f`

```fortran
subroutine state_space_f(N, M, cross_cov, cov, innov, trans, rc)
```

Constructs Kalman-filter innovation and transition matrices.
Wraps `state_space` in `linalg.c`.

| Argument | Type | Intent | Description |
|---|---|---|---|
| `N` | `integer` | in | Number of state-space steps |
| `M` | `integer` | in | State dimension (2 or 3) |
| `cross_cov` | `real(c_double)((N-1)*M*M)` | in | Cross-covariance matrices |
| `cov` | `real(c_double)(N*M*M)` | in | Covariance matrices |
| `innov` | `real(c_double)(N*M*M)` | out | Innovation (Cholesky) matrices |
| `trans` | `real(c_double)((N-1)*M*M)` | out | Transition matrices |
| `rc` | `integer` | out | 0 = success; i+1 = step i failed; −1 = M unsupported |

### `hyp_llp1_f`

```fortran
subroutine hyp_llp1_f(llp1_re, llp1_im, m, nz, zvals, out, rc)
```

Evaluates the Gauss hypergeometric function ₂F₁(−l, l+1; 1+m; z) over an
array of z values for fixed complex l(l+1) and integer m.
Wraps `hyp_llp1` in `smerfs.c`.

| Argument | Type | Intent | Description |
|---|---|---|---|
| `llp1_re/im` | `real(c_double)` | in | Real/imag parts of l(l+1) |
| `m` | `integer` | in | Integer order (≥ 0) |
| `nz` | `integer` | in | Number of z values |
| `zvals` | `real(c_double)(nz)` | in | Real z values |
| `out` | `complex(c_double_complex)(nz)` | out | ₂F₁ values |
| `rc` | `integer` | out | 0 = success; 1 = failure (non-convergence) |

### `hyp_lmz_f`

```fortran
function hyp_lmz_f(l, m, z, gamma_ratio0, psi0) result(res)
```

Evaluates ₂F₁(−l, l+1; 1+m; z) for large z (0.5 < z < 1) using the
psi/logarithmic series expansion (AMS55 15.3.11).
Wraps `hyp_lmz` in `smerfs.c`.

| Argument | Type | Intent | Description |
|---|---|---|---|
| `l` | `complex(c_double_complex)` | in | Complex degree l |
| `m` | `integer` | in | Integer order (≥ 0) |
| `z` | `real(c_double)` | in | Real argument (0.5 < z < 1) |
| `gamma_ratio0` | `complex(c_double_complex)` | in | Precomputed Γ(m+1)Γ(m)/Γ(m−l)/Γ(m+l+1) |
| `psi0` | `complex(c_double_complex)` | in | Precomputed ψ(1)+ψ(1+m)−ψ(m−l)−ψ(m+l+1) |
| `res` | `complex(c_double_complex)` | — | ₂F₁ value (function result) |

### `zigg_f`

```fortran
subroutine zigg_f(num_needed, num_ints, rand_ints, out, rc)
```

Generates float-precision standard normal variates using the Ziggurat method
(Marsaglia & Tsang 2000).  Wraps `zigg` in `ziggurat.c`.

| Argument | Type | Intent | Description |
|---|---|---|---|
| `num_needed` | `integer` | in | Number of normals requested |
| `num_ints` | `integer` | in | Number of 32-bit random integers supplied |
| `rand_ints` | `integer(c_int32_t)(num_ints)` | in | Raw random integers |
| `out` | `real(c_float)(num_needed)` | out | Standard normal variates |
| `rc` | `integer` | out | 0 = all produced; > 0 = number still outstanding |

## Example: `examples_f/example1_f.F90`

Fortran translation of `examples/example1.py`.  Reproduces the two outputs
of the Python script as CSV files instead of plots.

### What it does

```
nz = 128   ! iso-latitude rings
nphi = 256 ! phi points
coeffs = (1.0, 0.0, 1e-4)   ! C_l = 1/(1 + 1e-4*(l(l+1))^2)
```

1. **Partial fraction decomposition** of the power-spectrum denominator
   (closed-form quadratic formula for `M = 2`).

2. **Covariance tables** — evaluates ₂F₁ via `hyp_llp1_f` and accumulates
   `cov` and `cross_cov` via `update_cov_f` for each partial-fraction pole.

3. **State-space filter** — calls `state_space_f` for each Fourier m-mode
   to obtain innovation matrices (`innov`) and transition matrices (`trans`).

4. **Realisation** — generates complex Gaussian noise via `zigg_f`, walks the
   Kalman state-space model from the equator upward and then downward, and
   applies an inverse real FFT to obtain the spatial field.

5. **Analytic covariance** — evaluates
   C(cos θ) = Σ_l (2l+1) C_l P_l(cos θ) / 4π
   using the three-term Legendre recurrence.

### Memory layout convention

All arrays passed to the C kernels use **reversed Fortran dimension order** so
that Fortran column-major storage byte-matches C row-major layout without
transposing or packing:

```
Python/C  array[m][n][p][q]   declared shape (m_max+1, N, M, M)
Fortran   array(q, p, n, m)   declared as    (M, M, N, 0:m_max)
```

### Compiler note

The `-ieee=full` flag (nagfor) is required to prevent the NAG runtime from
aborting on transient NaN/Inf values that arise in intermediate steps of the
complex hypergeometric series before the series converges.
