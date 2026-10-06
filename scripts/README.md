# Reconstructing Figure 3

**Source:** L. Blackmore, B. Açıkmeşe, J. M. Carson III, *Lossless convexification
of control constraints for a class of nonlinear optimal control problems*,
Systems & Control Letters 61 (2012), 863–870.
DOI: https://doi.org/10.1016/j.sysconle.2012.04.010

## Result and scope

The main figure is an independently **solved, fixed-time MILP trajectory at
490 s**, not a digitized trajectory masquerading as a simulation. With the
figure-consistent settings below, it differs from the 15 published Figure 3
vertices by approximately **47 m RMS and 103 m maximum**, on an 80 km descent.
PDF graphics have finite drawing precision; these comparisons are not against
the authors' original numerical arrays. The inputs were inferred by examining
the published figures, so this is a reconstruction, not a blind replication.

It is **not** a complete reproduction of the paper's free-time optimum. The
reconstructed model gives a lower discrete cost near **451.9 s**. Nor does the
matching trajectory exactly satisfy the original Euclidean upper thrust bound:
the inferred outer polygon permits up to **10.048386 m/s²**. Both limitations
are prominently reported, rather than hidden by rescaling or clipping the plot.

## Run

From the repository root:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r scripts/requirements-figure.txt
python scripts/regenerate_figure.py
python -m unittest discover -s tests -v
python -m http.server 8000
# Open http://localhost:8000/figure.html
```

`regenerate_figure.py` also works from another working directory. It reads the
included `assets/figure/reference.json`; the PDF is not needed for routine
regeneration. HiGHS ships with SciPy. No commercial license, API key, optimizer
initial guess, or external plotting CDN is required. The website also works
when opened directly from disk and its SVG figures do not require JavaScript.

For optional re-extraction, place the supplied eight-page PDF at
`lossless convexification.pdf`, then run:

```sh
python scripts/extract_paper_reference.py
python scripts/regenerate_figure.py
```

The extractor is specific to this PDF's vector drawing order and checks the
expected paths. Its output records the source PDF's SHA-256. It extracts:

- Figure 2's gravity curves, including an inverse-square fit to the smooth
  curve (approximately 9.7795 m/s² and radius 100.165 km, rounded for the model).
- Figure 3's 15 vertices, **for validation only**.
- Figure 4's initial downrange velocity, approximately 199.987 m/s.
- A cropped original Figure 3 image for the comparison panel.

`solve()` receives only the gravity configuration, interval count, initial
velocity, polygon choice and final time. It never receives the Figure 3
reference. Do not interpolate the Figure 2 blue drawing vertices to construct
the field: a rendered polyline can connect samples across tangent-region
switches. The model uses the tangent construction described in Section 5.

## Stated versus inferred settings

See [`figure_config.json`](figure_config.json) for machine-readable provenance.

| Setting | Written description | Figure-consistent reconstruction |
|---|---|---|
| Initial position | [20,000, 80,000] m | Same |
| Target position/velocity | Zero | Same |
| Initial velocity | Zero | [200, 0] m/s, from Fig. 4 |
| Grid | N = 15 steps | 15 plotted states, 14 intervals of 35 s |
| Drag / mass depletion | Both zero | Same; control normalized to acceleration |
| Slack bounds | 2 and 10 m/s² | Same |
| Smooth gravity constants | Not listed | g0 = 9.78 m/s², R = 100,000 m, inferred from Fig. 2 |
| PWL regions | Tangent lines; numerical partition not listed | Tangents at 0,20,40,60,80,100 km; switches at 10,30,50,70,90 km |
| Norm approximation | Inner 32-gon | Outer 32-gon, inferred from Fig. 4's component/norm plateaus |
| Final time | Reported optimum 490 s | Fixed at 490 s for the main comparison; searched independently as a diagnostic |

The first/last gravity zones are clipped to the plotted domain [0, 100] km;
the outer zones' tangent centers therefore lie at the domain endpoints rather
than at the clipped interval midpoints. These choices explain the gravity
profiles closely but remain inferences, not settings explicitly supplied by
the authors. No parameter is adjusted by minimizing Figure 3 position error.

The generator also solves two alternatives at 490 s:
1. Change **only** the polygon from outer to inner (same facet orientations).
2. Additionally use the text's zero initial velocity and 15 intervals.

## Discrete optimization

Write `u = φ/m` in acceleration units, with constant mass. This normalizes the
fuel proxy by a constant factor. Set `Δt = tf / intervals`. Appendix (A.2) gives:

```text
r[k+1] = r[k] + Δt v[k] + 0.5 Δt² (u[k] + ghat(r[k]))
v[k+1] = v[k] + Δt (u[k] + ghat(r[k]))
J      = Δt sum(sigma[k])
```

Gravity is frozen at the **start** of each interval. These updates are not
exact integrations of the smooth inverse-square field over an interval.
Position plots connect state samples with straight lines, matching Figure 3.
The velocity equality is enforced as an equality: Appendix (A.7) appears to
print the same inequality direction twice. Initial conditions and all control
intervals, including `k=0`, are explicitly included despite inconsistent index
ranges in the printed Problem 5.

For downward gravity `g_y(h) = -g0 / (1+h/R)^2`, each region uses the tangent
`g_y(c) + g_y'(c)(h-c)`. Binary `z[k,i]` selects exactly one region. Instead of
loose big-M constraints, a disaggregated (convex-hull) MILP encoding uses:

```text
sum_i z[k,i] = 1
L[i] z[k,i] <= w[k,i] <= H[i] z[k,i]
h[k] = sum_i w[k,i]
g_y[k] = sum_i (a[i] w[k,i] + b[i] z[k,i])
z[k,i] in {0,1}
```

This selects the same region-wise affine dynamics as a correctly bounded
big-M encoding. At shared boundaries, either adjacent region may be selected;
the tangent model can have small jumps. No monotonic-altitude constraint is
added. Bounds are altitude [0,100] km, downrange [-100,100] km, and velocity
components [-2000,2000] m/s. Non-altitude numerical envelopes are inactive at
the reported solution. Position and velocity decision variables are scaled
by 10,000 m and 100 m/s for numerical conditioning; exported values use SI.

Let `n_j = (cos(2πj/32), sin(2πj/32))`, for `j=0,...,31`:

- **Outer:** `n_j^T u <= sigma`, `2 <= sigma <= 10`. The Euclidean norm can
  reach `sigma / cos(π/32)`, so this is *not* an exact SOCP constraint.
- **Inner:** `n_j^T u <= sigma cos(π/32)`. This enforces `norm(u) <= sigma`,
  but `sigma >= 2` alone does not ensure `norm(u) >= 2` after transcription.

The output reports actual Euclidean norms, slack/norm gaps and both original
bound violations. The continuous-time losslessness theorem is not invoked as
a guarantee for the discretized, polygon-approximated solutions.

## Optimality and validation

Each fixed-time solve requests HiGHS relative MILP gap `1e-8` and a 60 s time
limit. A timeout or other nonoptimal termination raises an error; it is never
published as an optimum. The final-time diagnostic scans 380–620 s at 10 s
spacing and applies golden-section refinement to sampled local-minimum
brackets, to 0.05 s width. This is a best-found bounded search, **not a global
certificate** for the possibly non-unimodal mixed-integer value function.

Automated checks include independent dynamics/terminal residuals, inferred
gravity-region consistency, polygon feasibility, actual norm bounds/slack
gaps, comparison against published vertices, zero-velocity/inner-polygon
alternatives, and agreement of the website's generated JS and JSON data.
The per-solve solver status, gap, objective and checks are in `results.json`.
No tests pretend that the small outer-polygon bound violation is zero.

## Files

- `regenerate_figure.py`: MILP, time search, comparison and artifact generation.
- `figure_config.json`: inferred settings and provenance.
- `extract_paper_reference.py`: optional PDF vector-path extraction.
- `render_figure.py`: scientific SVGs generated from solver results.
- `requirements-figure.txt`: Python dependencies.
- `../tests/test_figure.py`: numerical and artifact tests.
- `../assets/figure/reference.json`: PDF-derived validation/model evidence.
- `../assets/figure/results.json` and `results.js`: identical generated payloads.
- `../assets/figure/*.svg`: generated trajectory, overlay, profiles, alternatives
  and final-time-search figures.

The original cropped figure is attributed to the paper and is used for
scholarly comparison. This reconstruction does not claim author endorsement.
