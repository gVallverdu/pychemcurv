#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append("../")

import unittest
import numpy as np
from pychemcurv import VertexAtom


__author__ = "Germain Salvato-Vallverdu"
__copyright__ = "University of Pau and Pays Adour"
__email__ = "germain.vallverdu@univ-pau.fr"

class VertexAtomTests(unittest.TestCase):

    def setUp(self):
        self.theta_sp3 = np.arccos(-1 / 3)
        self.theta_sp2 = np.pi / 2
        l = 1.3

        # sp3 pyramid
        coords = [[0, 0, -l * np.cos(self.theta_sp3)]]
        IB = l * np.sin(self.theta_sp3)
        for angle in [0, 2 * np.pi / 3, 4 * np.pi / 3]:
            coords.append([IB * np.cos(angle), IB * np.sin(angle), 0])
        coords = np.array(coords, dtype=np.float64)
        self.va_sp3 = VertexAtom(coords[0], coords[1:])

        # sp2 case
        coords = [[0, 0, 0]]
        for angle in [0, 2 * np.pi / 3, 4 * np.pi / 3]:
            coords.append([l * np.cos(angle), l * np.sin(angle), 0])
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
        self.assertAlmostEqual(self.va_sp2.pyrA, 0)
        self.assertAlmostEqual(self.va_rand.pyrA, 18.7104053164)
        
    def test_angular_defect(self):
        ang_sp3 = np.degrees(np.arccos(- 1 / 3))
        self.assertAlmostEqual(self.va_sp3.angular_defect,
                              360 - 3 * ang_sp3)
        self.assertAlmostEqual(self.va_sp2.angular_defect, 0)
        self.assertAlmostEqual(self.va_rand.angular_defect, 29.83127456)

    def test_spherical_curvature(self):
        self.assertAlmostEqual(self.va_sp3.spherical_curvature, 0.5128205128205)
        self.assertTrue(np.isnan(self.va_sp2.spherical_curvature))
        self.assertAlmostEqual(self.va_rand.spherical_curvature, 0.4523719038)

    def test_improper(self):
        self.assertAlmostEqual(self.va_sp3.improper, -35.2643896828)
        self.assertAlmostEqual(self.va_sp2.improper, 0.)
        self.assertAlmostEqual(self.va_rand.improper, -30.021240733)


if __name__ == "__main__":
    unittest.main()