import numpy as np
import sys
import logging
import struct
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import rcParams



def display_profiles(gamma, m, pos, u, ax):
    fontsize=18
    
    font = {'family' : 'serif',
            'serif': ['Times'],
            'size'   : fontsize}
    
    rcParams['text.latex.preamble'] = \
        [r"\usepackage{amsmath,bm}\newcommand{\sbf}[1]{\textsf{\textbf{#1}}}"]
    rc('text', usetex=True)
    rc('font', **font)

    r = np.linalg.norm(pos, axis=1)
    hist, edges = np.histogram(r, bins='auto')
    nbins = len(hist)
    rho_c = np.zeros(nbins)
    u_c   = np.zeros(nbins)
    r_c   = np.zeros(nbins)
    for i in range(len(hist)):
        vol      = 4 * np.pi / 3. * (edges[i+1]**3 - edges[i]**3)
        rho_c[i] = m * hist[i] / vol
        r_c[i]   = (edges[i+1] + edges[i]) / 2.
        ig = np.where((r < edges[i+1]) & (r > edges[i]))
        u_c[i] = np.mean(u[ig[0]])

    ax1 = ax[0]
    ax2 = ax[1]
    ax1.minorticks_on()
    ax2.minorticks_on()
    ax1.plot(r_c, rho_c, 'go', label='Stretched particles')
    ax1.legend(loc='best', numpoints=1, title=r'$\gamma=$%.1f'%(gamma))
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylabel(r"$\rho$")
    ax2.plot(r_c, u_c, 'go', label='Stretched particles')
    ax2.set_xlabel(r"$r$")
    ax2.set_ylabel(r"$u$")


def save_particles(ids, pos, vel, mass, u, outfile, format):

    if format == 0:
        # Openning file
        try:
            ofile = open(outfile,'w')
        except IOError as e:
            msg = "IO Error({0}): {1}".format(e.errno, e.strerror)
            logging.warning(msg)
        except:
            print "Unexpected error:", sys.exc_info()[0]
            raise

        # Preparing every line to print to the file
        for i in range(len(pos)):
            # Formatting particle attributes
            m  = '% 3.8e' % mass[i]
            rx = '% 3.8e' % pos[i][0]
            ry = '% 3.8e' % pos[i][1]
            rz = '% 3.8e' % pos[i][2]
            ie = '% 3.8e' % u[i]

            # Right-align the strings
            outstring = "{0} {1} {2} {3} {4}\n".format( m.rjust(12),\
                                                    rx.rjust(12),\
                                                    ry.rjust(12),\
                                                    rz.rjust(12),\
                                                    ie.rjust(12))
            # Write to file
            ofile.write(outstring)

        # Closing the file
        ofile.close()

    elif format == 1:
        Ngas = len(mass)
        Npart = np.array([Ngas, 0, 0, 0, 0, 0])
        Nmass = np.array([0, 0, 0, 0, 0, 0])
        # Linearizing the 3D-array of the position and velocity
        pos = pos.ravel()
        vel = vel.ravel()
        
        dummy = np.zeros(Npart[0])
        time = 0.
        redshift = 0.0  # double
        flag_sfr = 0  # long
        flag_feedback = 0  # long
        bytesleft = 256 - 6*4 - 6*8 - 8 - 8 - 2*4 - 6*4
        fill = np.zeros(int(bytesleft/4.0), dtype=np.int)  # int
        
        with open(outfile, 'wb') as f:
            nbytes = 256
            # Header
            f.write(struct.pack('i', nbytes))
            f.write(struct.pack('i' * len(Npart), *Npart))
            f.write(struct.pack('d' * len(Nmass), *Nmass))
            f.write(struct.pack('d', time))
            f.write(struct.pack('d', redshift))
            f.write(struct.pack('i', flag_sfr))
            f.write(struct.pack('i', flag_feedback))
            f.write(struct.pack('i' * len(Npart), *Npart))
            f.write(struct.pack('i' * len(fill), *fill))
            f.write(struct.pack('i', nbytes))
        
            # Positions
            nbytes = int(len(pos) * 4)
            f.write(struct.pack('i', nbytes))
            f.write(struct.pack('f' * len(pos), *pos))
            f.write(struct.pack('i', nbytes))
        
            # Velocities
            f.write(struct.pack('i', nbytes))
            f.write(struct.pack('f' * len(vel), *vel))
            f.write(struct.pack('i', nbytes))
        
            # Ids
            nbytes = int(len(ids) * 4)
            f.write(struct.pack('i', nbytes))
            f.write(struct.pack('i' * len(ids), *ids))
            f.write(struct.pack('i', nbytes))
        
            # Masses
            nbytes = len(mass) * 4
            f.write(struct.pack('i', nbytes))
            f.write(struct.pack('f' * len(mass), *mass))
            f.write(struct.pack('i', nbytes))
        
            # Energy
            nbytes = len(u) * 4
            f.write(struct.pack('i', nbytes))
            f.write(struct.pack('f' * len(u), *u))
            f.write(struct.pack('i', nbytes))




