"""Independent numerical and generated-artifact checks. Run with unittest."""
import json
from pathlib import Path
import sys
import unittest
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit
import xml.etree.ElementTree as ET

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from regenerate_figure import compare, solve  # noqa: E402


class LandingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads((ROOT / 'scripts/figure_config.json').read_text())
        cls.reference = json.loads((ROOT / 'assets/figure/reference.json').read_text())
        cls.bundle = json.loads((ROOT / 'assets/figure/results.json').read_text())
        cls.s = solve(490, cls.config['gravity'], intervals=14, vx0=200, polygon='outer')

    def test_boundary_conditions_and_time_grid(self):
        s = self.s
        np.testing.assert_allclose(s['position_m'][0], [20000, 80000], atol=1e-7)
        np.testing.assert_allclose(s['velocity_m_s'][0], [200, 0], atol=1e-7)
        np.testing.assert_allclose(s['position_m'][-1], [0, 0], atol=1e-7)
        np.testing.assert_allclose(s['velocity_m_s'][-1], [0, 0], atol=1e-7)
        np.testing.assert_allclose(np.diff(s['time_s']), 35, atol=1e-12)
        self.assertEqual(len(s['position_m']), 15)
        self.assertEqual(len(s['control_m_s2']), 14)

    def test_dynamics_independently(self):
        s = self.s
        r, v, u = (np.array(s[key]) for key in ('position_m', 'velocity_m_s', 'control_m_s2'))
        a = u.copy()
        a[:, 1] += np.array(s['gravity_y_m_s2'])
        np.testing.assert_allclose(np.diff(r, axis=0), 35 * v[:-1] + .5 * 35**2 * a, atol=1e-6)
        np.testing.assert_allclose(np.diff(v, axis=0), 35 * a, atol=1e-7)
        self.assertGreater(np.max(r[:, 0]), 31000)  # Correct initial rightward bend.
        self.assertTrue(np.all(r[:, 1] >= -1e-7))

    def test_tangent_gravity_and_active_regions(self):
        model = self.config['gravity']
        edges = np.array(model['altitude_m'])
        centers = np.array(model['tangent_altitudes_m'])
        for r, gy, i in zip(self.s['position_m'], self.s['gravity_y_m_s2'], self.s['gravity_segment']):
            h, c = r[1], centers[i]
            self.assertGreaterEqual(h, edges[i] - 1e-6)
            self.assertLessEqual(h, edges[i + 1] + 1e-6)
            radius, g0 = model['radius_m'], model['surface_gravity_m_s2']
            expected = -g0 / (1 + c / radius)**2 + 2 * g0 / radius / (1 + c / radius)**3 * (h - c)
            self.assertAlmostEqual(gy, expected, places=8)

    def test_outer_polygon_and_honest_norm_violation(self):
        s = self.s
        u, sigma = np.array(s['control_m_s2']), np.array(s['sigma_m_s2'])
        angles = np.arange(32) * 2 * np.pi / 32
        normals = np.column_stack((np.cos(angles), np.sin(angles)))
        self.assertLessEqual(np.max(u @ normals.T - sigma[:, None]), 1e-7)
        self.assertGreaterEqual(np.min(sigma), 2 - 1e-7)
        self.assertLessEqual(np.max(sigma), 10 + 1e-7)
        norms = np.linalg.norm(u, axis=1)
        np.testing.assert_allclose(norms, s['thrust_magnitude_m_s2'], atol=1e-10)
        self.assertGreater(np.max(norms), 10.04)  # Do not silently clip to ten.
        self.assertLessEqual(np.max(norms), 10 / np.cos(np.pi / 32) + 1e-7)
        self.assertAlmostEqual(s['cost_m_s'], 35 * sum(sigma), places=6)
        self.assertLessEqual(s['checks']['mip_relative_gap'], 1e-8)

    def test_against_reference_not_an_input(self):
        result = compare(self.s, self.reference)
        self.assertLess(result['max_node_error_m'], 200)
        self.assertLess(result['rms_node_error_m'], 100)
        self.assertGreater(result['max_node_error_m'], 1)  # Not copied coordinates.
        self.assertEqual(self.reference['figure3']['node_count'], 15)
        self.assertAlmostEqual(self.reference['figure4']['initial_downrange_velocity_m_s'], 200, delta=1)

    def test_inner_polygon_alternative(self):
        s = solve(490, self.config['gravity'], polygon='inner')
        self.assertLessEqual(max(s['thrust_magnitude_m_s2']), 10 + 1e-7)
        self.assertGreater(s['cost_m_s'], self.s['cost_m_s'])
        self.assertGreater(compare(s, self.reference)['max_node_error_m'], 1000)
        # Finite transcription/polygon approximation need not preserve the lower bound.
        violation = max(0, 2 - min(s['thrust_magnitude_m_s2']))
        self.assertAlmostEqual(violation, s['checks']['lower_norm_violation_m_s2'])

    def test_written_parameters_alternative(self):
        s = self.bundle['written_parameters_solution']
        self.assertEqual(s['intervals'], 15)
        np.testing.assert_allclose(s['velocity_m_s'][0], [0, 0], atol=1e-7)
        np.testing.assert_allclose(s['position_m'][-1], [0, 0], atol=1e-7)
        self.assertEqual(s['polygon'], 'inner')
        self.assertLessEqual(max(np.array(s['position_m'])[:, 0]), 20000 + 1e-6)

    def test_time_search_discrepancy_is_reported(self):
        search = self.bundle['time_search']
        best = search['best']
        self.assertTrue(440 < best['tf_s'] < 470)
        self.assertLess(best['cost_m_s'], self.s['cost_m_s'])
        self.assertIn('not a global', search['caveat'])
        self.assertAlmostEqual(best['cost_m_s'], min(e['cost_m_s'] for e in search['evaluations'] if e['cost_m_s'] is not None))

    def test_generated_artifacts_are_current(self):
        np.testing.assert_allclose(self.s['position_m'], self.bundle['solution']['position_m'], atol=1e-3)
        self.assertEqual(self.config, self.bundle['config'])
        js = (ROOT / 'assets/figure/results.js').read_text()
        prefix = 'window.FIGURE_RESULTS = '
        self.assertTrue(js.startswith(prefix))
        self.assertEqual(json.loads(js[len(prefix):].rstrip().removesuffix(';')), self.bundle)
        self.assertEqual(self.reference, self.bundle['reference'])
        for name in ['regenerated-figure-3', 'comparison', 'profiles', 'alternatives', 'time-search']:
            root = ET.parse(ROOT / f'assets/figure/{name}.svg').getroot()
            self.assertEqual(root.tag, '{http://www.w3.org/2000/svg}svg')

    def test_invalid_configuration_rejected(self):
        with self.assertRaises(ValueError):
            solve(0, self.config['gravity'])
        with self.assertRaises(ValueError):
            solve(490, self.config['gravity'], polygon='circle')


class PageTests(unittest.TestCase):
    def test_local_assets_and_accessible_images(self):
        class Parser(HTMLParser):
            def __init__(self):
                super().__init__()
                self.links = []
                self.images = []

            def handle_starttag(self, tag, attrs):
                attrs = dict(attrs)
                for attr in ('src', 'href'):
                    if attr in attrs:
                        self.links.append(attrs[attr])
                if tag == 'img':
                    self.images.append(attrs)

        parser = Parser()
        parser.feed((ROOT / 'figure.html').read_text())
        for link in parser.links:
            path = urlsplit(link)
            if path.scheme or not path.path:
                continue
            relative = unquote(path.path)
            if relative == 'lossless convexification.pdf':
                continue  # Optional local source PDF; not required for regeneration.
            self.assertTrue((ROOT / relative).exists(), relative)
        self.assertGreaterEqual(len(parser.images), 6)
        for image in parser.images:
            self.assertTrue(image.get('alt'))


if __name__ == '__main__':
    unittest.main()
