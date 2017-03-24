# gadget-polytrope

Generates a spherical particle distribution with a polytropic profile. 
It is oriented for producing initial conditions for GADGET, although
it is possible to obtain the positions in ASCII format as well.

The code first solves the Lane-Emden equations with a given index and 
then it radially stretches a close-packed sphere of particles to match
the corresponding polytropic density profile. Based on the final radial
position, the code also assigns an internal energy to each
particle to reproduce the polytropic temperature profile.

This script is **under development**, thus comments and contributions
are very much welcome.

# Author

Felipe G. Goicovic

# Usage

The basic usage is
```bash
./polyt-ics -g GAMMA -n NUM
```
where GAMMA is the polytropic index and NUM is the number of particles desired.
This will produce a binary file called 'Polytrope.dat' that works as initial
conditions for GADGET (format 1).

For a complete description of parameters use
```bash
./polyt-ics -h
```


# Licence 

MIT Licence

