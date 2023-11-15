#!/usr/bin/env python
# coding: utf-8

# import sys
# sys.path.append("../")

from pytest import approx
from pychemcurv import VertexAtom, TrivalentVertex, POAV1, POAV2
import numpy as np

__author__ = "Germain Salvato-Vallverdu"
__copyright__ = "University of Pau and Pays Adour"
__email__ = "germain.vallverdu@univ-pau.fr"


class TestVertexAtom:
    """ Test for class pychemcurv.core.VertexAtom """

    def setup_method(self):

        self.theta_sp3 = np.arccos(-1 / 3)
        self.theta_sp2 = np.pi / 2
        self.CC_bond = 1.3

        # sp3 pyramid
        coords = [[0, 0, -self.CC_bond * np.cos(self.theta_sp3)]]
        IB = self.CC_bond * np.sin(self.theta_sp3)
        for angle in [0, 2 * np.pi / 3, 4 * np.pi / 3]:
            coords.append([IB * np.cos(angle), IB * np.sin(angle), 0])
        coords = np.array(coords, dtype=np.float64)
        self.va_sp3 = VertexAtom(coords[0], coords[1:])

        # squared pyramid
        theta = np.radians(100.0)
        coords = [[0, 0, -self.CC_bond * np.cos(theta)]]
        IB = self.CC_bond * np.sin(theta)
        for i in range(4):
            angle = i * np.pi / 2
            coords.append([IB * np.cos(angle), IB * np.sin(angle), 0])
        coords = np.array(coords, dtype=np.float64)
        self.va_sq = VertexAtom(coords[0], coords[1:])

        # random case:
        coords = [[-2.62985741,  6.99670582, -2.89817324],
                  [-2.32058737,  5.49122664, -3.13957301],
                  [-2.92519373,  6.96241176, -1.65009278],
                  [-1.62640146,  7.93539179, -3.17337668]]
        self.va_rand = VertexAtom(coords[0], coords[1:])

    def test_a_star_a_shape(self):
        assert self.va_sp3.a.shape == (3, )
        assert self.va_sq.a.shape == (3, )
        assert self.va_sp3.star_a.shape == (3, 3)
        assert self.va_sq.star_a.shape == (4, 3)

    def test_a_star_a_values(self):
        assert self.va_rand.a == approx(
            [-2.62985741,  6.99670582, -2.89817324])
        star_a = np.array([[-2.32058737,  5.49122664, -3.13957301],
                           [-2.92519373,  6.96241176, -1.65009278],
                           [-1.62640146,  7.93539179, -3.17337668]])
        assert self.va_rand.star_a.flatten() == approx(star_a.flatten())

    def test_reg_star_a(self):
        reg = np.array([[-2.43106709,  6.02902499, -3.05333841],
                        [-2.86004832,  6.96997636, -1.92539491],
                        [-1.91379553,  7.66654812, -3.09455724]]).flatten()
        assert self.va_rand.reg_star_a.flatten() == approx(reg)

    def test_reg_normal(self):
        assert self.va_rand.reg_normal == approx(
            [-0.81987531,  0.24597365, -0.51701203])

    def test_angles(self):
        angles = self.va_sp3.get_angles()
        for a in angles.values():
            assert a == approx(np.arccos(-1/3))

        angles = self.va_sp3.get_angles(radians=False)
        for a in angles.values():
            assert a == approx(np.degrees(np.arccos(-1/3)))

    def test_distances(self):
        for d in self.va_sp3.distances:
            assert d == approx(self.CC_bond)

        for d in self.va_sq.distances:
            assert d == approx(self.CC_bond)

    def test_normal(self):
        assert self.va_sp3.normal == approx([0., 0., 1.])
        assert self.va_rand.normal == approx(
            [-0.80739002,  0.22175107, -0.54676121])

    def test_pyr_distance(self):
        assert self.va_sp3.pyr_distance == approx(
            -self.CC_bond * np.cos(self.theta_sp3))
        assert self.va_sq.pyr_distance == approx(
            -self.CC_bond * np.cos(np.radians(100)))

    def test_from_pyramid(self):
        # sp3 pyramid
        va = VertexAtom.from_pyramid(self.CC_bond, self.theta_sp3,
                                     radians=True)
        assert va.a.shape == (3, )
        assert self.va_sp3.a == approx(va.a)
        assert va.star_a.shape == (3, 3)
        assert self.va_sp3.star_a.flatten() == approx(va.star_a.flatten())

        va = VertexAtom.from_pyramid(self.CC_bond, 100, 4, radians=False)
        assert va.a.shape == (3,)
        assert va.star_a.shape == (4, 3)
        assert self.va_sq.a == approx(va.a)
        assert self.va_sq.star_a.flatten() == approx(va.star_a.flatten())


class TestTrivalentVertex:
    """ Tests about the pychemcurv.TrivalentVertex class """

    def setup_method(self):

        self.theta_sp3 = np.arccos(-1 / 3)
        self.theta_sp2 = np.pi / 2
        self.CC_bond = 1.3

        # sp3 pyramid
        coords = [[0, 0, -self.CC_bond * np.cos(self.theta_sp3)]]
        IB = self.CC_bond * np.sin(self.theta_sp3)
        for angle in [0, 2 * np.pi / 3, 4 * np.pi / 3]:
            coords.append([IB * np.cos(angle), IB * np.sin(angle), 0])
        coords = np.array(coords, dtype=np.float64)
        self.va_sp3 = TrivalentVertex(coords[0], coords[1:])

        # sp2 case
        coords = [[0, 0, 0]]
        for angle in [0, 2 * np.pi / 3, 4 * np.pi / 3]:
            coords.append([self.CC_bond * np.cos(angle), 
                           self.CC_bond * np.sin(angle),
                           0])
        coords = np.array(coords, dtype=np.float64)
        self.va_sp2 = TrivalentVertex(coords[0], coords[1:])

        # random case:
        coords = [[-2.62985741,  6.99670582, -2.89817324],
                  [-2.32058737,  5.49122664, -3.13957301],
                  [-2.92519373,  6.96241176, -1.65009278],
                  [-1.62640146,  7.93539179, -3.17337668]]
        self.va_rand = TrivalentVertex(coords[0], coords[1:])

    def test_pyrA(self):
        assert self.va_sp3.pyrA == approx(np.degrees(self.theta_sp3) - 90.)
        assert self.va_sp3.pyrA_r == approx(self.theta_sp3 - np.pi / 2)
        assert self.va_sp2.pyrA == approx(0)
        assert self.va_rand.pyrA == approx(18.7104053164)
        assert self.va_rand.pyrA_r == approx(np.radians(18.7104053164))

    def test_angular_defect(self):
        assert self.va_sp3.angular_defect == approx(
            2 * np.pi - 3 * np.arccos(- 1 / 3))
        assert self.va_sp2.angular_defect == approx(0)
        assert np.degrees(self.va_rand.angular_defect) == approx(29.83127456)

    def test_spherical_curvature(self):
        assert self.va_sp3.spherical_curvature == approx(0.5128205128205)
        assert np.isnan(self.va_sp2.spherical_curvature)
        assert self.va_rand.spherical_curvature == approx(0.4523719038)

    def test_improper(self):
        assert self.va_sp3.improper == approx(np.radians(-35.2643896828))
        assert self.va_sp2.improper == approx(0.)
        assert self.va_rand.improper == approx(np.radians(-30.021240733))

    def test_pyr_distance(self):
        assert self.va_sp2.pyr_distance == approx(0.)
        dist = self.CC_bond * np.sin(self.theta_sp3 - np.pi / 2)
        assert self.va_sp3.pyr_distance == approx(dist)
        assert self.va_rand.pyr_distance == approx(0.4515551342307116)

    def test_from_pyramid(self):
        # sp3 pyramid
        va = VertexAtom.from_pyramid(self.CC_bond, self.theta_sp3,
                                     radians=True)
        assert self.va_sp3.a.shape == (3, )
        assert self.va_sp3.a == approx(va.a)
        assert self.va_sp3.star_a.shape == (3, 3)
        assert self.va_sp3.star_a.flatten() == approx(va.star_a.flatten())

        # sp2 case
        va = VertexAtom.from_pyramid(self.CC_bond, self.theta_sp2,
                                     radians=True)
        assert self.va_sp2.a.shape == (3, )
        assert self.va_sp2.a == approx(va.a)
        assert self.va_sp2.star_a.shape == (3, 3)
        assert self.va_sp2.star_a.flatten() == approx(va.star_a.flatten())


class TestPOAV1:
    """ Tests about the pychemcurv.core.POAV1 class """

    def setup_method(self):
        theta = np.degrees(np.arccos(-1 / 3))
        self.pyrA_sp3 = theta - 90
        v_sp3 = TrivalentVertex.from_pyramid(1.3, theta)
        v_sp2 = TrivalentVertex.from_pyramid(1.3, 90.)
        v_r = TrivalentVertex.from_pyramid(1.3, 108.2)

        self.poav_sp3 = POAV1(v_sp3)
        self.poav_sp2 = POAV1(v_sp2)
        self.poav_a = POAV1(v_r)

        # random coordinates of a pyramid with pyrA = 0.326558177 radians
        coords = [[-2.62985741,  6.99670582, -2.89817324],
                  [-2.32058737,  5.49122664, -3.13957301],
                  [-2.92519373,  6.96241176, -1.65009278],
                  [-1.62640146,  7.93539179, -3.17337668]]
        self.poav_b1 = POAV1(VertexAtom(coords[0], coords[1:]))

    def test_pyrA(self):
        assert self.poav_sp3.pyrA == approx(self.pyrA_sp3)
        assert self.poav_sp2.pyrA == approx(0.)
        assert self.poav_a.pyrA == approx(18.2)

        assert self.poav_sp3.pyrA_r == approx(np.arccos(-1 / 3) - np.pi / 2)
        assert self.poav_sp2.pyrA_r == approx(0.)
        assert self.poav_a.pyrA_r == approx(np.radians(18.2))

        assert self.poav_b1.pyrA_r == approx(0.326558177)
        assert self.poav_b1.pyrA == approx(18.71040530758611)

    def test_coeffs(self):
        poavs = [self.poav_sp3, self.poav_sp2, self.poav_a]
        c_pis = [1/2, 0., 0.46496976]
        for poav, c_pi in zip(poavs, c_pis):
            assert poav.c_pi == approx(c_pi)

        lambda_pi2s = [3/4, 1., 0.78380312]
        for poav, l_pi2 in zip(poavs, lambda_pi2s):
            assert poav.lambda_pi**2 == approx(l_pi2)

        ns = [3, 2, 2.82749177]
        ms = [1/3, 0, 0.27583059]
        for poav, n, m in zip(poavs, ns, ms):
            assert poav.n == approx(n)
            assert poav.m == approx(m)

    def test_POAV(self):
        assert self.poav_sp3.poav == approx([0., 0., 1.])
        assert self.poav_b1.poav == approx(
            [-0.81987531,  0.24597365, -0.51701203])


class TestPOAV2:
    """ Tests about the pychemcurv.core.POAV2 class """

    def setup_method(self):
        # random coordinates of a pyramid with pyrA = 0.326558177 radians
        coords = [[-2.62985741,  6.99670582, -2.89817324],
                  [-2.32058737,  5.49122664, -3.13957301],
                  [-2.92519373,  6.96241176, -1.65009278],
                  [-1.62640146,  7.93539179, -3.17337668]]
        self.poav = POAV2(TrivalentVertex(coords[0], coords[1:]))

    def test_matrix(self):
        m = np.array([[-0.23175585, -0.12713934,  0.49598426],
                      [0.04802619,  0.47612631,  0.02444729],
                      [0.18372966, -0.34898698, -0.52043155]])
        assert self.poav.matrix.flatten() == approx(m.flatten())

    def test_pyrA(self):
        assert self.poav.pyrA_r == approx([1.80097239, 1.75115937, 2.09343774])
        assert self.poav.pyrA == approx(
            [103.18811722, 100.33404089, 119.94514689])

    def test_u_pi(self):
        assert self.poav.u_pi == approx(
            [-0.91097524,  0.11226749, -0.39688806])

    def test_sigma_hyb_nbrs(self):
        assert self.poav.sigma_hyb_nbrs == approx((4.602500609983161,
                                                   7.444750673988855,
                                                   0.961463263653626))

    def test_pi_hyb_nbr(self):
        assert self.poav.pi_hyb_nbr == approx(0.23956910356339622)
