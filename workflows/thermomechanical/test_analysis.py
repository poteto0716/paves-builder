import unittest

import numpy as np

from analysis import fit_modulus, fit_tg, hyperbola_slope, smooth_hyperbola


class AnalysisTest(unittest.TestCase):
    def test_density_hyperbola_recovers_tg_and_cte(self):
        temperature = np.linspace(200.0, 800.0, 601)
        density = smooth_hyperbola(temperature, 1.1, 420.0, -0.00025, -0.00045, 15.0)
        result = fit_tg(temperature, density, 10.0)
        self.assertAlmostEqual(result['tg_K'], 420.0, delta=0.2)
        self.assertAlmostEqual(result['density_asymptote_slope_glass_g_cm3_K'], -0.00025, delta=2e-7)
        self.assertAlmostEqual(result['density_asymptote_slope_rubber_g_cm3_K'], -0.00070, delta=2e-7)
        # CTEs: the true curve's own slope and density at Tg -/+ 50 K.
        truth = (1.1, 420.0, -0.00025, -0.00045, 15.0)
        self.assertEqual((result['cte_glass_temperature_K'], result['cte_rubber_temperature_K']),
                         (result['tg_K'] - 50, result['tg_K'] + 50))
        for phase in ('glass', 'rubber'):
            at = result[f'cte_{phase}_temperature_K']
            slope = float(hyperbola_slope(at, *truth))
            self.assertAlmostEqual(result[f'density_slope_{phase}_g_cm3_K'], slope, delta=1e-6)
            expected = -slope / (3 * float(smooth_hyperbola(np.asarray([at]), *truth)[0]))
            self.assertAlmostEqual(result[f'alpha_linear_{phase}_1_K'], expected, delta=4e-7)
        self.assertAlmostEqual(result['cte_mixing_fraction'], 0.5 * (1 - 50 / np.hypot(50, 15)), delta=1e-3)
        self.assertTrue(result['cte_temperatures_within_data'])

    def test_broad_transition_keeps_the_sign_of_the_visible_slope(self):
        # A fit whose glass asymptote turned positive (density rising with T) although
        # the curve falls everywhere in the data: the CTEs follow the curve, not the asymptote.
        temperature = np.repeat(np.arange(800.0, 99.0, -20.0), 25)
        density = smooth_hyperbola(temperature, 1.10, 450.0, -1.5e-4, -3e-4, 300.0)
        density = density + np.random.default_rng(50).normal(0, .008, len(temperature))
        result = fit_tg(temperature, density, 10.0)
        self.assertGreater(result['density_asymptote_slope_glass_g_cm3_K'], 0)
        self.assertGreater(result['alpha_linear_glass_1_K'], 0)
        self.assertGreater(result['alpha_linear_rubber_1_K'], result['alpha_linear_glass_1_K'])

    def test_linear_stress_fit_recovers_modulus(self):
        strain = np.linspace(0.0, 0.02, 21)
        result = fit_modulus(strain, 2500.0 * strain + 3.0)
        self.assertAlmostEqual(result['tensile_modulus_GPa'], 2.5, places=12)
        self.assertAlmostEqual(result['intercept_MPa'], 3.0, places=12)


if __name__ == '__main__':
    unittest.main()
