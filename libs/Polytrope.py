#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Polytropic profiles"""

__author__ = "Felipe G. Goicovic"
__copyright__ = "Copyright 2016, Felipe G. Goicovic"
__licence__ = "MIT"
__email__ = "fagarri1@uc.cl"

import numpy as np
import sys
import matplotlib.pyplot as plt

class PolytropicProfile:
    """ Class that computes the polytropic profile with a given index by
    solving the Lane-Emden equation

    Atributes
    ---------
    n: polytropic index
    gamma: polytropic exponent (1+1/n)
    x: dimensionless radial coordinate
    theta: dimensionless solution of Lane-Emden equation
    rho: density profile
    P: pressure profile
    u: internal energy profile
    """

    def __init__(self, n=3., gamma=None):

        #To initialize the profile it is necessary to specify either 'n' or
        #'gamma'. If both are given, the script will use the latter.

        if gamma != None:
            if gamma <= 1.2:
                print("ERROR: gamma = %.2f is not allowed, please choose a \
                      value larger than 1.2"%(gamma))
                sys.exit()
            n = 1 / (gamma - 1.)
        else:
            gamma = 1. + 1. / n

        if n <= 0 or n >= 5.:
            print("ERROR: n = %.2f is not allowed, please choose between 0 \
                  and 5"%(n))
            sys.exit()

        self.n = float(n)
        self.gamma = gamma

        self.theta = np.array([0.])
        self.rho = np.array([0.])

    def solve_lane_emden(self, h=0.01, method="runge"):
        """
        Arguments
        ---------
        h: step of the integration
        method: Only allowed 'runge' (4th order Runge-kutta)
        """

        def dtheta_dx(x, mu):
            if x > 0.:
                return -mu / x**2.
            else:
                return 0.
        def dmu_dx(x, theta):
            if theta >= 0:
                return x**2 * theta**self.n
            else:
                return 0.

        n_points = int(np.floor(10. / h))

        theta = np.full(n_points, -1.)
        mu    = np.zeros(n_points)
        x     = np.zeros(n_points)
        theta[0] = 1.
        print("Solving Lane-Emden equation...")
        if method == "runge":
            print("Method: {0} with integration step of {1}".format("Runge-Kutta",
                                                                    h))
            i = 0
            while(theta[i] > 0.):
                # only add more elements if needed
                if i >= len(theta) - 1:
                    theta = np.append(theta, np.zeros(n_points)-1.)
                    mu = np.append(mu, np.zeros(n_points))
                    x = np.append(x, np.zeros(n_points))
                k1 = dtheta_dx(x[i], mu[i])
                c1 = dmu_dx(x[i], theta[i])
                theta_aux = theta[i] + k1 * h / 2.
                mu_aux = mu[i] + c1 * h / 2.
                x_aux = x[i] + h / 2.
                k2 = dtheta_dx(x_aux, mu_aux)
                c2 = dmu_dx(x_aux, theta_aux)
                theta_aux = theta[i] + k2 * h / 2.
                mu_aux = mu[i] + c2 * h / 2.
                k3 = dtheta_dx(x_aux, mu_aux)
                c3 = dmu_dx(x_aux, theta_aux)
                theta_aux = theta[i] + k3 * h / 2.
                mu_aux = mu[i] + c3 * h / 2.
                x_aux = x[i] + h
                k4 = dtheta_dx(x_aux, mu_aux)
                c4 = dmu_dx(x_aux, theta_aux)
                k = (k1 + 2. * k2 + 2. * k3 + k4) / 6.
                c = (c1 + 2. * c2 + 2. * c3 + c4) / 6.
                theta[i+1] = theta[i] + k * h
                mu[i+1] = mu[i] + c * h
                x[i+1] = x_aux
                i += 1
        else:
            print("ERROR: Method %s not allowed"%(method))
            sys.exit()
        idx = np.where(theta >= 0.)[0]
        self.theta = theta[idx]
        self.mu    = mu[idx]
        self.x     = x[idx]
        print("Done with Lane-Emden.")

    def compute_profile(self, G=1., M=1., R=1.):
        """
        Arguments
        ---------
        G: gravitational constant
        M: total mass of the sphere
        R: radius of the sphere
        """
        if self.theta[0] == 0.:
            print("Warning: Running 'solve_lane_emden()' with default parameters.")
            self.solve_lane_emden()

        theta = self.theta
        mu = self.mu
        x = self.x
        m_n = M / mu[-1]
        r_n = R / x[-1]
        rho_c = m_n / (4*np.pi*r_n**3)
        P_c = G*m_n**2 / (4*np.pi*(self.n+1)*r_n**4)

        self.rho = rho_c * theta**self.n
        self.P   = P_c * theta**(self.n + 1.)
        self.u   = 1 / (self.gamma - 1.) * self.P / self.rho

    def plot_profile(self, field='rho', G=1., M=1., R=1., axes=None, **kargs):
        if self.rho[0] == 0.:
            compute_profile(G=G, M=M, R=R)

        x   = self.x
        x1  = x[-1]
        rho = self.rho
        P   = self.P
        u   = self.u

        if axes == None:
            fig = plt.figure()
            axes = fig.add_subplot(111)
        if field == "rho":
            axes.plot(x/x1, rho, '-', **kargs)
        elif field == "P":
            axes.plot(x/x1, P, '-', **kargs)
        elif field == "u":
            axes.plot(x/x1, u, '-', **kargs)

        return axes

