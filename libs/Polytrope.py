# POLYTROPIC PROFILES
#
# Developed by Felipe G. Goicoivc (fagarri1@uc.cl)
# Nov 2016
#
import numpy as np
import sys
import matplotlib.pyplot as plt

class profile:

  def __init__(self, n = 3.,gamma=None):
    if gamma != None:
      if gamma <= 1.2:
        print("ERROR: gamma = %.2f is not allowed, please choose a value larger than 1.2"%(gamma))
        sys.exit()
      n = 1 / (gamma - 1.)
    else:
      gamma = 1. + 1. / n
    if n <= 0 or n >= 5.:
      print("ERROR: n = %.2f is not allowed, please choose between 0 and 5"%(n))
      sys.exit()
    self.n = float(n)
    self.gamma = gamma 

    self.theta = np.array([0.])

  def solve_lane_emden(self, h = 0.01, method = "runge"):

    def dtheta_dx(x,mu):
      if x > 0.:
        return -mu/x**2.
      else:
        return 0.
    def dmu_dx(x,theta):
      if theta >= 0:
        return x**2*theta**self.n
      else:
        return 0.

    n_points = int(np.floor(10. / h))

    theta = np.zeros(n_points) - 1.
    mu = np.zeros(n_points)
    x = np.zeros(n_points)
    theta[0] = 1.
    mu[0] = 0.
    x[0] = 0.
    print("Solving Lane-Emden equations...")
    if method == "runge":
      i = 0
      while(theta[i] > 0.):
        if i >= len(theta)-1:
          theta = np.append(theta,np.zeros(n_points)-1.)
          mu = np.append(mu,np.zeros(n_points))
          x = np.append(x,np.zeros(n_points))
        k1 = dtheta_dx(x[i],mu[i])
        c1 = dmu_dx(x[i],theta[i])
        theta_aux = theta[i] + k1 * h / 2.
        mu_aux = mu[i] + c1 * h / 2.
        x_aux = x[i] + h / 2.
        k2 = dtheta_dx(x_aux,mu_aux)
        c2 = dmu_dx(x_aux,theta_aux)
        theta_aux = theta[i] + k2 * h / 2.
        mu_aux = mu[i] + c2 * h / 2.
        k3 = dtheta_dx(x_aux,mu_aux)
        c3 = dmu_dx(x_aux,theta_aux)
        theta_aux = theta[i] + k3 * h / 2.
        mu_aux = mu[i] + c3 * h / 2.
        x_aux = x[i] + h 
        k4 = dtheta_dx(x_aux,mu_aux)
        c4 = dmu_dx(x_aux,theta_aux)
        k = (k1 + 2.*k2 + 2.*k3 + k4) / 6.
        c = (c1 + 2.*c2 + 2.*c3 + c4) / 6.
        theta[i+1] = theta[i] + k*h
        mu[i+1] = mu[i] + c*h
        x[i+1] = x_aux
        i += 1
    else:
      print("ERROR: Method %s not allowed"%(method))
      sys.exit()
    idx = np.where(theta >= 0.)
    theta = theta[idx[0]]
    mu = mu[idx[0]]
    x = x[idx[0]]
    self.theta = theta
    self.mu = mu
    self.x = x 

  def get_profile(self, G=1., M=1., R=1.):
    if self.theta[0] == 0.:
      print("ERROR: Please run 'solve_lane_emden()' first.")
      sys.exit()

    theta = self.theta
    x = self.x
    dtheta_x1 = np.fabs((theta[-2]-theta[-1])/(x[-2]-x[-1]))
    x1 = x[-1]
    rho_c = M*x1/(4*np.pi*R**3*dtheta_x1)
    P_c = G*M**2/(4*np.pi*(self.n+1)*dtheta_x1**2*R**4)
    
    rho = rho_c*theta**self.n
    P = P_c*theta**(self.n+1.)
    return x/x1, rho, P

  def plot_profile(self, field='rho', G=1., M=1., R=1., axes=None, **kargs):
    if self.theta[0] == 0.:
      print("ERROR: Please run 'solve_lane_emden()' first.")
      sys.exit()

    theta = self.theta
    x = self.x
    dtheta_x1 = np.fabs((theta[-2]-theta[-1])/(x[-2]-x[-1]))
    x1 = x[-1]
    rho_c = M*x1/(4*np.pi*R**3*dtheta_x1)
    P_c = G*M**2/(4*np.pi*(self.n+1)*dtheta_x1**2*R**4)
    
    #rho = rho_c*theta**self.n
    rho = rho_c*theta**self.n
    P = P_c*theta**(self.n+1.)
    u  = 1 / (self.gamma - 1.) * P/rho
    
    if axes == None:
      fig = plt.figure()
      axes = fig.add_subplot(111)
    if field == "rho":
      axes.plot(x/x1,rho,'-',**kargs)
    elif field == "P":
      axes.plot(x/x1,P,'-',**kargs)
    elif field == "u":
      axes.plot(x/x1, u,'-',**kargs)

    return axes
