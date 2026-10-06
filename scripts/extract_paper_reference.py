#!/usr/bin/env python3
"""Extract model input (Fig. 2) and independent validation (Fig. 3) from the PDF.

Uses vector paths, not raster tracing. Drawing indices and axis limits are
specific to the supplied eight-page PDF. Fig. 2 supports the inferred gravity configuration;
Fig. 3 coordinates are used exclusively AFTER optimization for comparison.
"""
import hashlib
import json
from pathlib import Path

import fitz
import numpy as np
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'assets/figure'


def vertices(drawing):
    items = drawing['items']
    assert all(item[0] == 'l' for item in items), 'Expected a vector polyline'
    return np.array([tuple(items[0][1])] + [tuple(item[2]) for item in items])


def coordinates(points, rect, xlim, ylim):
    return np.column_stack((
        xlim[0] + (points[:, 0] - rect.x0) / rect.width * (xlim[1] - xlim[0]),
        ylim[1] - (points[:, 1] - rect.y0) / rect.height * (ylim[1] - ylim[0]),
    ))


def main():
    pdf = ROOT / 'lossless convexification.pdf'
    digest = hashlib.sha256(pdf.read_bytes()).hexdigest()
    doc = fitz.open(pdf)
    assert len(doc) == 8
    p5, p6 = doc[4].get_drawings(), doc[5].get_drawings()
    assert p5[204]['color'] == (0, 0, 1) and len(p5[204]['items']) == 6
    assert p6[49]['color'] == (0, 0, 1) and len(p6[49]['items']) == 14
    gravity = coordinates(vertices(p5[204]), p5[202]['rect'], [0, 100000], [2, 10])
    trajectory = coordinates(vertices(p6[49]), p6[47]['rect'], [-40000, 60000], [-10000, 90000])
    velocity = coordinates(vertices(p6[291]), p6[290]['rect'], [0, 500], [-500, 500])
    # Fit only the gravity curve, never the trajectory. The config uses rounded
    # physical parameters, and remains a separately documented modeling choice.
    smooth_full = coordinates(vertices(p5[203]), p5[202]['rect'], [0, 100000], [2, 10])
    model = lambda h, g0, radius: g0 / (1 + h / radius)**2
    fit, _ = curve_fit(model, smooth_full[:, 0], smooth_full[:, 1], p0=[9.8, 100000])
    fit_error = np.max(np.abs(model(smooth_full[:, 0], *fit) - smooth_full[:, 1]))
    smooth = smooth_full[::40]
    result = {
        'source': {
            'title': 'Lossless convexification of control constraints for a class of nonlinear optimal control problems',
            'authors': 'Blackmore, Açıkmeşe & Carson (2012)',
            'doi': '10.1016/j.sysconle.2012.04.010',
            'pdf_sha256': digest,
            'method': 'Vector paths calibrated against printed axis limits; finite PDF drawing precision.',
        },
        'gravity': {
            'provenance': 'Blue PWL curve in Figure 2, printed page 867. Supports inference of region boundaries; the solver uses tangent lines in figure_config.json, not linear interpolation between these displayed vertices.',
            'altitude_m': gravity[:, 0].tolist(),
            'magnitude_m_s2': gravity[:, 1].tolist(),
            'smooth_reference': smooth.tolist(),
            'inverse_square_fit': {
                'surface_gravity_m_s2': float(fit[0]), 'radius_m': float(fit[1]),
                'max_curve_residual_m_s2': float(fit_error),
                'note': 'Fit to the red Figure 2 vector path; includes finite PDF drawing precision. Config rounds to 9.78 m/s^2 and 100 km.',
            },
        },
        'figure3': {
            'provenance': 'Blue trajectory in Figure 3, printed page 868. Validation only, never passed to solve().',
            'position_m': trajectory.tolist(),
            'node_count': len(trajectory),
            'estimated_drawing_precision_m': 50,
        },
        'figure4': {
            'initial_downrange_velocity_m_s': float(velocity[0, 1]),
            'last_plotted_time_s': float(velocity[-1, 0]),
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'reference.json').write_text(json.dumps(result, indent=2) + '\n')
    doc[5].get_pixmap(matrix=fitz.Matrix(3, 3), clip=fitz.Rect(88, 87, 259, 255)).save(OUT / 'original-figure-3.png')
    print(f'Extracted {len(gravity)} gravity knots and {len(trajectory)} trajectory reference nodes to {OUT}')
    print('Figure 4 initial velocity:', velocity[0, 1], 'm/s')


if __name__ == '__main__':
    main()
