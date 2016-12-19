import numpy as np


BOLTZMANN = 1.38066e-16
PROTON = 1.6726e-24
mean_weigth = 4./(8.-5.*(1.-0.76))
Unit_Mass = 1.989e33
Unit_Length = 6.957e10
Unit_Velocity = np.sqrt(6.67e-8 / Unit_Length * Unit_Mass)
Unit_Time = Unit_Length / Unit_Velocity
