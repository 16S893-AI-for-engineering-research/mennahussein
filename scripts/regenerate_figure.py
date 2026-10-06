#!/usr/bin/env python3
"""Reconstruct the paper's discrete landing experiment, not its trajectory pixels.

Run from any directory: python scripts/regenerate_figure.py
See scripts/README.md for provenance, assumptions, and numerical limitations.
"""
from pathlib import Path
import argparse
import json
import math

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parents[1]


class InfeasibleProblem(RuntimeError):
    """HiGHS certified infeasibility, not merely a timeout or solver failure."""


def solve(tf, gravity, intervals=14, vx0=200.0, time_limit=60.0, polygon='inner'):
    """Fixed-time MILP; SI output, scaled position/velocity decision variables.

    Uses Appendix A dynamics (correcting the printed velocity-equality typo).
    A disaggregated/hull encoding replaces its loose big-M gravity encoding.
    The outer polygon reproduces the figures, but differs from the text's
    inner polygon; neither is asserted to inherit continuous-time losslessness.
    There are no constraints or objectives involving the Figure 3 reference.
    """
    n = intervals
    if n < 1 or int(n) != n or tf <= 0:
        raise ValueError('tf must be positive and intervals a positive integer')
    h = np.asarray(gravity['altitude_m']) / 1e4
    if 'tangent_altitudes_m' in gravity:
        centers = np.asarray(gravity['tangent_altitudes_m'])
        g0, radius = gravity['surface_gravity_m_s2'], gravity['radius_m']
        values = -g0 / (1 + centers / radius)**2
        slopes = 2 * g0 / radius / (1 + centers / radius)**3 * 1e4
        intercepts = values - slopes * centers / 1e4
    else:
        g = -np.asarray(gravity['magnitude_m_s2'])
        slopes = np.diff(g) / np.diff(h)
        intercepts = g[:-1] - slopes * h[:-1]
    nz = len(slopes)
    count = 0

    def variables(shape):
        nonlocal count
        size = math.prod(shape)
        result = np.arange(count, count + size).reshape(shape)
        count += size
        return result

    r = variables((n + 1, 2))
    v = variables((n + 1, 2))
    u = variables((n, 2))
    sigma = variables((n,))
    gy = variables((n,))
    z = variables((n, nz))
    w = variables((n, nz))  # altitude contribution of each active PWL segment
    lower = np.full(count, -np.inf)
    upper = np.full(count, np.inf)
    lower[r], upper[r] = -10, 10
    lower[r[:, 1]], upper[r[:, 1]] = h[0], h[-1]
    lower[v], upper[v] = -20, 20
    lower[u], upper[u] = -10, 10
    lower[sigma], upper[sigma] = 2, 10
    lower[gy], upper[gy] = -10, 0
    lower[z], upper[z] = 0, 1
    lower[w], upper[w] = 0, h[-1]
    lower[r[0]] = upper[r[0]] = [2, 8]
    lower[v[0]] = upper[v[0]] = [vx0 / 100, 0]
    lower[r[-1]] = upper[r[-1]] = 0
    lower[v[-1]] = upper[v[-1]] = 0
    integrality = np.zeros(count)
    integrality[z] = 1
    rows, cols, vals, lbs, ubs = [], [], [], [], []

    def constraint(terms, lb=0, ub=0):
        row = len(lbs)
        for col, value in terms:
            rows.append(row)
            cols.append(int(col))
            vals.append(value)
        lbs.append(lb)
        ubs.append(ub)

    dt = tf / n
    angles = np.arange(32) * 2 * np.pi / 32
    normals = np.column_stack((np.cos(angles), np.sin(angles)))
    if polygon == 'inner':
        # Inscribed regular polygon: vertices lie on the circle, not outside it.
        normals /= np.cos(np.pi / 32)
    elif polygon != 'outer':
        raise ValueError('polygon must be inner or outer')
    for k in range(n):
        for axis in range(2):
            terms = [(r[k + 1, axis], 1), (r[k, axis], -1),
                     (v[k, axis], -dt / 100), (u[k, axis], -.5 * dt**2 / 1e4)]
            if axis == 1:
                terms.append((gy[k], -.5 * dt**2 / 1e4))
            constraint(terms)
            terms = [(v[k + 1, axis], 1), (v[k, axis], -1), (u[k, axis], -dt / 100)]
            if axis == 1:
                terms.append((gy[k], -dt / 100))
            constraint(terms)
        constraint([(j, 1) for j in z[k]], 1, 1)
        constraint([(r[k, 1], -1)] + [(j, 1) for j in w[k]])
        constraint([(gy[k], 1)] + [(w[k, i], -slopes[i]) for i in range(nz)]
                   + [(z[k, i], -intercepts[i]) for i in range(nz)])
        for i in range(nz):
            constraint([(w[k, i], 1), (z[k, i], -h[i])], 0, np.inf)
            constraint([(w[k, i], 1), (z[k, i], -h[i + 1])], -np.inf, 0)
        for a in normals:
            constraint([(u[k, 0], a[0]), (u[k, 1], a[1]), (sigma[k], -1)], -np.inf, 0)

    A = coo_matrix((vals, (rows, cols)), shape=(len(lbs), count)).tocsc()
    objective = np.zeros(count)
    objective[sigma] = dt
    result = milp(objective, integrality=integrality, bounds=Bounds(lower, upper),
                  constraints=LinearConstraint(A, lbs, ubs),
                  options={'time_limit': time_limit, 'mip_rel_gap': 1e-8})
    if result.status == 2:
        raise InfeasibleProblem(f'Infeasible at tf={tf}')
    if not result.success:
        raise RuntimeError(f'MILP failed at tf={tf}: {result.message}')
    x = result.x
    position, velocity, control = x[r] * 1e4, x[v] * 100, x[u]
    magnitudes = np.linalg.norm(control, axis=1)
    acceleration = control + np.column_stack((np.zeros(n), x[gy]))
    position_residual = np.diff(position, axis=0) - dt * velocity[:-1] - .5 * dt**2 * acceleration
    velocity_residual = np.diff(velocity, axis=0) - dt * acceleration
    return {
        'tf_s': float(tf), 'intervals': n, 'initial_velocity_m_s': [vx0, 0], 'polygon': polygon,
        'time_s': np.linspace(0, tf, n + 1).tolist(),
        'position_m': position.tolist(), 'velocity_m_s': velocity.tolist(),
        'control_m_s2': control.tolist(), 'sigma_m_s2': x[sigma].tolist(),
        'thrust_magnitude_m_s2': magnitudes.tolist(), 'gravity_y_m_s2': x[gy].tolist(),
        'gravity_segment': np.argmax(x[z], axis=1).tolist(),
        'cost_m_s': float(result.fun),
        'checks': {
            'solver': 'SciPy milp / HiGHS', 'status': result.message,
            'mip_relative_gap': float(result.mip_gap),
            'terminal_position_error_m': float(np.max(np.abs(position[-1]))),
            'terminal_velocity_error_m_s': float(np.max(np.abs(velocity[-1]))),
            'max_position_dynamics_residual_m': float(np.max(np.abs(position_residual))),
            'max_velocity_dynamics_residual_m_s': float(np.max(np.abs(velocity_residual))),
            'min_thrust_m_s2': float(magnitudes.min()),
            'max_thrust_m_s2': float(magnitudes.max()),
            'max_sigma_minus_norm_m_s2': float(np.max(x[sigma] - magnitudes)),
            'max_abs_sigma_minus_norm_m_s2': float(np.max(np.abs(x[sigma] - magnitudes))),
            'upper_norm_violation_m_s2': float(max(0, magnitudes.max() - 10)),
            'lower_norm_violation_m_s2': float(max(0, 2 - magnitudes.min())),
            'max_polygon_violation_m_s2': float(max(0, np.max(control @ normals.T - x[sigma, None]))),
            'max_gravity_residual_m_s2': float(np.max(np.abs(x[gy] - (slopes[np.argmax(x[z], axis=1)] * position[:-1, 1] / 1e4 + intercepts[np.argmax(x[z], axis=1)])))),
        },
    }


def compare(solution, reference):
    """Node-to-node error, not a fit objective or point-to-curve distance."""
    actual = np.asarray(solution['position_m'])
    expected = np.asarray(reference['figure3']['position_m'])
    if actual.shape != expected.shape:
        return None
    error = np.linalg.norm(actual - expected, axis=1)
    return {
        'node_error_m': error.tolist(),
        'rms_node_error_m': float(np.sqrt(np.mean(error**2))),
        'max_node_error_m': float(error.max()),
        'max_error_percent_initial_altitude': float(error.max() / 80000 * 100),
        'reference_precision_note': 'PDF coordinates have finite drawing precision (order tens of metres), not original author data.',
    }


def time_search(config):
    """Coarse scan + golden-section refinement of sampled local minima.

    Each fixed-time MILP is certified to tolerance. This outer search is NOT a
    global certificate: the mixed-integer value function need not be unimodal.
    """
    cache = {}

    def evaluate(tf):
        key = float(tf)
        if key not in cache:
            try:
                cache[key] = solve(tf, config['gravity'], config['intervals'],
                                   config['initial_downrange_velocity_m_s'],
                                   polygon=config['polygon'])
            except InfeasibleProblem:
                cache[key] = None
        return cache[key]['cost_m_s'] if cache[key] else np.inf

    grid = np.arange(380.0, 621.0, 10.0)
    costs = [evaluate(t) for t in grid]
    ratio = (math.sqrt(5) - 1) / 2
    for i in range(1, len(grid) - 1):
        if not (np.isfinite(costs[i]) and costs[i] <= min(costs[i - 1], costs[i + 1])):
            continue
        lo, hi = grid[i - 1], grid[i + 1]
        c, d = hi - ratio * (hi - lo), lo + ratio * (hi - lo)
        while hi - lo > .05:
            if evaluate(c) < evaluate(d):
                hi, d = d, c
                c = hi - ratio * (hi - lo)
            else:
                lo, c = c, d
                d = lo + ratio * (hi - lo)
        evaluate((lo + hi) / 2)
    feasible = [s for s in cache.values() if s is not None]
    if not feasible:
        raise RuntimeError('No feasible final time in the search window')
    best = min(feasible, key=lambda s: s['cost_m_s'])
    return {
        'method': '10 s coarse scan on [380, 620] s, then golden-section refinement of sampled local minima to 0.05 s bracket width.',
        'caveat': 'Best found in this bounded search, not a global continuous-time optimum. Does not recover the paper\'s reported 490 s optimum.',
        'evaluations': [{'tf_s': t, 'cost_m_s': s['cost_m_s'] if s else None}
                        for t, s in sorted(cache.items())],
        'best': best,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT / 'scripts/figure_config.json')
    parser.add_argument('--output', type=Path, default=ROOT / 'assets/figure')
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    reference = json.loads((ROOT / 'assets/figure/reference.json').read_text())
    solution = solve(config['reported_final_time_s'], config['gravity'], config['intervals'],
                     config['initial_downrange_velocity_m_s'], polygon=config['polygon'])
    inner = solve(config['reported_final_time_s'], config['gravity'], config['intervals'],
                  config['initial_downrange_velocity_m_s'], polygon='inner')
    # Isolate the polygon change above; separately try the written N and velocity.
    written = solve(config['reported_final_time_s'], config['gravity'], 15, 0, polygon='inner')
    import scipy
    bundle = {
        'schema_version': 1,
        'software': {'numpy': np.__version__, 'scipy': scipy.__version__},
        'config': config, 'solution': solution, 'comparison': compare(solution, reference),
        'inner_polygon_solution': inner, 'written_parameters_solution': written,
        'time_search': time_search(config),
        'reference': reference,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(bundle, indent=2, allow_nan=False)
    (args.output / 'results.json').write_text(payload + '\n')
    # A local JS bundle also works when the website is opened with file://.
    (args.output / 'results.js').write_text('window.FIGURE_RESULTS = ' + payload + ';\n')
    from render_figure import render
    render(bundle, args.output)
    print(json.dumps({'comparison': bundle['comparison'], 'checks': solution['checks'],
                      'best_searched_time_s': bundle['time_search']['best']['tf_s']}, indent=2))


if __name__ == '__main__':
    main()
