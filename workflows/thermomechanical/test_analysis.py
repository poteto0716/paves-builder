import unittest

import numpy as np

from analysis import fit_modulus, fit_tg, smooth_hyperbola


class AnalysisTest(unittest.TestCase):
    def test_density_hyperbola_recovers_tg_and_cte(self):
        temperature = np.linspace(200.0, 800.0, 601)
        density = smooth_hyperbola(temperature, 1.1, 420.0, -0.00025, -0.00045, 15.0)
        result = fit_tg(temperature, density, 10.0)
        self.assertAlmostEqual(result['tg_K'], 420.0, delta=0.2)
        self.assertAlmostEqual(result['density_slope_glass_g_cm3_K'], -0.00025, delta=2e-7)
        self.assertAlmostEqual(result['density_slope_rubber_g_cm3_K'], -0.00070, delta=2e-7)
        expected_glass = 0.00025 / (3 * result['density_at_tg_g_cm3'])
        self.assertAlmostEqual(result['alpha_linear_glass_1_K'], expected_glass, delta=5e-8)

    def test_linear_stress_fit_recovers_modulus(self):
        strain = np.linspace(0.0, 0.02, 21)
        result = fit_modulus(strain, 2500.0 * strain + 3.0)
        self.assertAlmostEqual(result['tensile_modulus_GPa'], 2.5, places=12)
        self.assertAlmostEqual(result['intercept_MPa'], 3.0, places=12)


if __name__ == '__main__':
    unittest.main()
