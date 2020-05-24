import unittest
import main_new as mn


class Tests(unittest.TestCase):

    # def test_props(self):
    #     self.assertAlmostEqual(800.9, mn.props('dens', 300))
    #     self.assertAlmostEqual(1966.9, mn.props('cp', 300))
    #     self.assertAlmostEqual(0.1112, mn.props('cond', 300))
    #     self.assertAlmostEqual(0.001569, mn.props('visc', 300))

    def test_props_300(self):
        self.assertAlmostEqual(800.9, mn.props('Density', 300))
        self.assertAlmostEqual(1966.9, mn.props('Cp', 300))
        self.assertAlmostEqual(0.1112, mn.props('Therm. Cond.', 300))
        self.assertAlmostEqual(0.0015686, mn.props('Viscosity', 300))


if __name__ == "main":
    unittest.main()
