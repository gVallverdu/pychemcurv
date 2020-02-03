#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append("../")

import unittest
import numpy as np
from pychemcurv import VertexAtom, Hybridization


__author__ = "Germain Salvato-Vallverdu"
__copyright__ = "University of Pau and Pays Adour"
__email__ = "germain.vallverdu@univ-pau.fr"

class VertexAtomTests(unittest.TestCase):
    """ Tests about the pychemcurv.VertexAtom class """

    def setUp(self):
        self.theta_sp3 = np.arccos(-1 / 3)
        self.theta_sp2 = np.pi / 2
        self.l = 1.3

        # sp3 pyramid
        coords = [[0, 0, -self.l * np.cos(self.theta_sp3)]]
        IB = self.l * np.sin(self.theta_sp3)
        for angle in [0, 2 * np.pi / 3, 4 * np.pi / 3]:
            coords.append([IB * np.cos(angle), IB * np.sin(angle), 0])
        coords = np.array(coords, dtype=np.float64)
        self.va_sp3 = VertexAtom(coords[0], coords[1:])

        # sp2 case
        coords = [[0, 0, 0]]
        for angle in [0, 2 * np.pi / 3, 4 * np.pi / 3]:
            coords.append([self.l * np.cos(angle), self.l * np.sin(angle), 0])
        coords = np.array(coords, dtype=np.float64)
        self.va_sp2 = VertexAtom(coords[0], coords[1:])

        # random case:
        coords = [[-2.62985741,  6.99670582, -2.89817324],
                  [-2.32058737,  5.49122664, -3.13957301],
                  [-2.92519373,  6.96241176, -1.65009278],
                  [-1.62640146,  7.93539179, -3.17337668]]
        self.va_rand = VertexAtom(coords[0], coords[1:])


    def test_pyrA(self):
        self.assertAlmostEqual(self.va_sp3.pyrA, 
                               np.degrees(self.theta_sp3) - 90.)
        self.assertAlmostEqual(self.va_sp3.pyrA_r,
                               self.theta_sp3 - np.pi / 2)
        self.assertAlmostEqual(self.va_sp2.pyrA, 0)
        self.assertAlmostEqual(self.va_rand.pyrA, 18.7104053164)
        self.assertAlmostEqual(self.va_rand.pyrA_r, np.radians(18.7104053164))
        
    def test_angular_defect(self):
        self.assertAlmostEqual(self.va_sp3.angular_defect,
                               2 * np.pi - 3 * np.arccos(- 1 / 3))
        self.assertAlmostEqual(self.va_sp2.angular_defect, 0)
        self.assertAlmostEqual(np.degrees(self.va_rand.angular_defect), 29.83127456)

    def test_spherical_curvature(self):
        self.assertAlmostEqual(self.va_sp3.spherical_curvature, 0.5128205128205)
        self.assertTrue(np.isnan(self.va_sp2.spherical_curvature))
        self.assertAlmostEqual(self.va_rand.spherical_curvature, 0.4523719038)

    def test_improper(self):
        self.assertAlmostEqual(self.va_sp3.improper, np.radians(-35.2643896828))
        self.assertAlmostEqual(self.va_sp2.improper, 0.)
        self.assertAlmostEqual(self.va_rand.improper, np.radians(-30.021240733))

    def test_pyr_distance(self):
        self.assertAlmostEqual(self.va_sp2.pyr_distance, 0.)
        self.assertAlmostEqual(self.va_sp3.pyr_distance, 
                               self.l * np.sin(self.theta_sp3 - np.pi / 2))
        self.assertAlmostEqual(self.va_rand.pyr_distance, 0.4515551342307116)

class HybridizationTests(unittest.TestCase):
    """ Tests about the pychemcurv.Hybridization class """

    def setUp(self):
        self.pyrA_sp3 = np.degrees(np.arccos(-1 / 3)) - 90
        self.h = Hybridization(pyrA=np.array([0, 18.2, self.pyrA_sp3]))

        coords = [[-2.62985741,  6.99670582, -2.89817324],
                  [-2.32058737,  5.49122664, -3.13957301],
                  [-2.92519373,  6.96241176, -1.65009278],
                  [-1.62640146,  7.93539179, -3.17337668]]
        self.h1 = Hybridization(vertex=VertexAtom(coords[0], coords[1:]))
        self.h2 = Hybridization(pyrA=0.326558177, radians=True)

    def test_pyrA(self):
        self.assertTrue(np.allclose(self.h.pyrA, [ 0., 18.2, 19.47122063]))
        self.assertTrue(np.allclose(self.h.pyrA_r, 
                                    [ 0., np.radians(18.2), np.arccos(-1 / 3) - np.pi / 2]))
        self.assertAlmostEqual(self.h1.pyrA, self.h2.pyrA)
        self.assertAlmostEqual(self.h1.pyrA_r, self.h2.pyrA_r)

    def test_coeffs(self):
        self.assertTrue(np.allclose(self.h.c_pi, [0, 0.46496976, 1/2]))
        self.assertTrue(np.allclose(self.h.lambda_pi**2, [1, 0.78380312, 3 / 4]))

        self.assertTrue(np.allclose(self.h.n, [2., 2.82749177, 3.]))
        self.assertTrue(np.allclose(self.h.m, [0, 0.27583059, 1/3]))

    def test_coeffs_struct(self):
        h1_dict = self.h1.as_dict()
        h2_dict = self.h2.as_dict()
        print(h1_dict)
        print(h2_dict)
        for key in h1_dict:
            self.assertAlmostEqual(h1_dict[key], h2_dict[key])



if __name__ == "__main__":
    unittest.main()