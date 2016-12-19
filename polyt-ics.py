#!/usr/bin/env python

import numpy as np
import struct
import matplotlib.pyplot as plt

from libs.Polytrope import *
from libs.OptionsParser import *
from libs.Utils import *

M = 1.
R = 1. 
G = 1.

op = OptionsParser()
op = op.get_args()
N = op.num

#solve the Lane-Emden equations
prof = profile(gamma=op.gamma)
prof.solve_lane_emden()
x, rho, P = prof.get_profile()

if op.plot:
    ax1 = prof.plot_profile(label='Analytic profile',lw=2)

Nside = int(np.floor( N**(1./3.) ))
Naux = Nside**3
pos = np.zeros( (2 * Naux, 3) ) #+100
h = 2 * R / Nside
for i in range(Nside):
    for j in range(Nside):
        for k in range(Nside):
            l = Nside*Nside*i + Nside*j + k
            pos[l,0]      = -R + h * (i + 0.25)
            pos[l,1]      = -R + h * (j + 0.25)
            pos[l,2]      = -R + h * (k + 0.25)
            pos[l+Naux,0] = -R + h * (i + 0.75)
            pos[l+Naux,1] = -R + h * (j + 0.75)
            pos[l+Naux,2] = -R + h * (k + 0.75)

r = np.linalg.norm(pos,axis=1)
ig = np.where(r<=R)
pos=pos[ig[0],:]
r = r[ig[0]]
Ngas = len(r)
print "Placed ",Ngas," particles in close-packed uniform sphere"
m = M/Ngas
vel = np.zeros((Ngas,3))
u = np.zeros(Ngas)
ids = np.arange(1,Ngas+1,step=1)
mass = m*np.ones_like(r)

r_n = np.array([0.])
i = 0
if Ngas > 100:
    n_bin = 0.01*Ngas
else:
    n_bin = 1.
r_n0 = n_bin / Ngas
while r_n[i] <= 1.:
    r_n = np.append(r_n,(r_n0+r_n[i]**3)**(1./3))
    i += 1

print "Stretching particles..."
r_n_new = np.zeros(len(r_n))
rho_n = np.zeros(len(r_n)-1)
for i in range(len(r_n)-1):
    idx = np.where((r<r_n[i+1]) & (r>=r_n[i]))
    ngas = len(idx[0])
    #print ngas
    ig = np.where(np.fabs(x-r_n_new[i]) == np.min(np.fabs(x-r_n_new[i])))
    vol = 4*np.pi/3.*(r_n[i+1]**3-r_n[i]**3)
    #print ngas, vol
    rho0 = m * ngas / vol
    r_n_new[i+1] = (rho0/rho[ig[0]]*(r_n[i+1]**3-r_n[i]**3)+r_n_new[i]**3)**(1./3)
    vol_new = 4*np.pi/3.*(r_n_new[i+1]**3-r_n_new[i]**3)
    if vol_new > 0.:
        rho_n[i] = m * ngas / vol_new

for i in range(len(r_n)-1):
    idx = np.where((r<r_n[i+1]) & (r>=r_n[i]))
    for j in idx[0]:
        r_new = r_n_new[i] + (r[j]-r_n[i])*(r_n_new[i+1]-r_n_new[i])/(r_n[i+1]-r_n[i])
        pos[j,0] *= r_new/r[j]
        pos[j,1] *= r_new/r[j]
        pos[j,2] *= r_new/r[j]

r_new = np.linalg.norm(pos,axis=1)

if op.plot:
    hist,edges = np.histogram(r_new, bins='auto')
    rho_c = np.zeros(len(hist))
    r_c = np.zeros(len(hist))
    for i in range(len(hist)):
        vol = 4*np.pi/3.*(edges[i+1]**3-edges[i]**3)
        rho_c[i] = m*hist[i]/vol
        r_c[i] = (edges[i+1]+edges[i])/2.
    
    ax1.plot(r_c,rho_c,'go',label='Stretched particles')
    ax1.legend(loc='best',numpoints=1)
    ax1.set_xlabel(r"$r/R$",fontsize=16)
    ax1.set_ylabel(r"$\rho$",fontsize=16)
    plt.show()


u_prof = 1 / (op.gamma-1) * P / rho 
for i in range(len(r_new)):
    ig = np.where(np.fabs(x-r_new[i]) == np.min(np.fabs(x-r_new[i])))
    u[i] = u_prof[ig[0]]

print "Writing output file..."
save_particles(pos, vel, mass, u, op.outfile, op.format)

print "done...bye!"

#if __name__=="__main__":
#    main()
