import argparse
import textwrap

class OptionsParser:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description='Generate a sphere of particles with a polytropic profile.',
                                              formatter_class=argparse.RawTextHelpFormatter)

        self.parser.add_argument("-g", "-gamma",
                            dest     = "gamma",
                            type     = float,
                            help     = "Polytropic index",
                            required = True)

        self.parser.add_argument("-N", "-n",
                            dest     = "num",
                            type     = int,
                            help     = "Number of particles",
                            required = True)

        self.parser.add_argument("-o",
                            metavar = "outfile",
                            dest    = "outfile",
                            help    = "Name of output file",
                            default = "Polytrope.dat")

        self.parser.add_argument("--plot",
                            dest    = "plot",
                            help    = "Plot analytic and actual radial profiles (histogram)",
                            action  = "store_true")

        self.parser.add_argument("-format",
                            dest    = "format",
                            type    = int,
                            help    = "Format of output file. 0 = ASCII, 1 = Gadget binary (format 1)",
                            default = 1 )

        self.parser.add_argument("-m", "-mass",
                            dest     = "mass",
                            type     = float,
                            help     = "Total gas mass",
                            default  = 1.)

        self.parser.add_argument("-r", "-radius",
                            dest     = "radius",
                            type     = float,
                            help     = "Radius of the sphere",
                            default  = 1.)

    def get_args(self):
        return self.parser.parse_args()
