# gadget-polytrope

Generates a spherical particle distribution with a polytropic profile. 
It is oriented for producing initial conditions for GADGET, although
it is possible to obtain the positions in ASCII format as well.


This script is **under development**.

# Author

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

