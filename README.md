# burgers-eq-fd — finite-difference solutions of Burgers' equation (MATLAB)

Solves the 1-D viscous/inviscid Burgers equation on a periodic domain

    u_t + u u_x = nu u_xx ,   0 <= x < L

with seven finite-difference schemes, compares them with exact solutions, and
measures conservation, shock behavior and order of accuracy.

Needs MATLAB R2016b or newer (implicit expansion, local functions in scripts).
No toolboxes are required.

## Quick start

```matlab
cd 2026/burgers_fd
burgers_setup                                   % adds src/ to the path

[x, u, info] = solve_burgers(@sin, 'N', 200, 'tFinal', 2, 'scheme', 'godunov');
plot(x, u)

runtests('tests')                               % run the test suite
run('examples/example_plateau_shock.m')         % figures go to results/
```

## Schemes

| name             | type                    | order (space/time) | notes |
|------------------|-------------------------|--------------------|-------|
| `upwind_nc`      | explicit, non-conservative | 1 / 1           | wrong shock speed, mass drifts: a cautionary example |
| `godunov`        | explicit, conservative  | 1 / 1              | monotone, entropy-satisfying; the robust shock choice |
| `lax_friedrichs` | explicit, conservative  | 1 / 1              | very diffusive; `nu = 0` only |
| `lax_wendroff`   | explicit, conservative  | 2 / 2              | Richtmyer two-step; oscillates at shocks |
| `maccormack`     | explicit, conservative  | 2 / 2              | predictor-corrector; oscillates at shocks |
| `ftcs`           | explicit, conservative  | 2 / 1              | viscous flow only (unstable when `nu = 0`) |
| `crank_nicolson` | implicit, conservative  | 2 / 2              | Newton iteration with sparse Jacobian; no CFL limit |

`solve_burgers` picks the time step from the current `max|u|` and `nu`:
`dt = cfl / (max|u|/dx + 2 nu/dx^2)` for explicit schemes (default `cfl = 0.5`),
`dt = cfl dx / max|u|` for Crank-Nicolson. It stops with `burgers:unstable` if the
solution blows up.

## Reference solutions (`src/exact`)

* `burgers_cole_hopf(x, t, nu)` — exact periodic viscous solution (Cole-Hopf transform).
* `burgers_characteristics(x, t, u0fun, umax)` — exact inviscid solution before the shock forms (`burgers_breaking_time` gives that time; 1 for `sin(x)`).
* `burgers_exact_plateau(x, t, ub, uh, xa, xb)` — exact rarefaction + shock for a raised plateau, valid until the fan catches the shock.

## Layout

    burgers_setup.m          path setup
    src/solve_burgers.m      driver (name-value options)
    src/get_scheme.m         scheme lookup
    src/schemes/step_*.m     one-step update functions
    src/exact/               exact solutions
    src/convergence_study.m  error and order measurement
    src/ic_*.m               initial conditions
    examples/                four runnable scripts that write figures to results/
    tests/                   function-based tests (runtests('tests'))

## Examples

| script | shows |
|--------|-------|
| `example_viscous_exact`   | three schemes vs the Cole-Hopf solution |
| `example_shock_formation` | `sin(x)` steepening into a shock; Godunov vs MacCormack |
| `example_plateau_shock`   | conservation, shock speed and oscillations on a rarefaction + shock |
| `example_convergence`     | observed orders of accuracy |

## Verification status

MATLAB was not available where this project was written, so **the `.m` files
have not been executed in MATLAB**. What was done instead:

* every scheme, exact solution and test condition was transcribed into NumPy and run there; the thresholds in `tests/` come from those runs (e.g. observed orders ~0.9 for Godunov, ~2.0 for MacCormack / Crank-Nicolson on the viscous problem; plateau shock position within 0.0075 of exact for the conservative schemes, 0.04 off for `upwind_nc`);
* the exact Cole-Hopf solution was checked against the PDE residual, the characteristics solution against its defining identity;
* the `.m` files were checked for balanced blocks, brackets and quotes with a script.

Expect to find a typo or two on first run; `runtests('tests')` is the first thing to try. Please report any failure messages.
