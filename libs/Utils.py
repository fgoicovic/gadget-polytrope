import numpy as np
import sys
import logging
import struct

def save_particles(pos, vel, mass, u, outfile, format):

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
            vx = '% 3.8e' % vel[i][0]
            vy = '% 3.8e' % vel[i][1]
            vz = '% 3.8e' % vel[i][2]

            # Right-align the strings
            outstring = "{0} {1} {2} {3} {4} {5} {6}\n".format( m.rjust(12),\
                                                               rx.rjust(12),\
                                                               ry.rjust(12),\
                                                               rz.rjust(12),\
                                                               vx.rjust(12),\
                                                               vy.rjust(12),\
                                                               vz.rjust(12))
            # Write to file
            ofile.write(outstring)

        # Closing the file
        ofile.close()

    elif format == 1:
        Ngas = len(mass)
        ids = np.arange(Ngas)
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




