import argparse
import textwrap

class OptionsParser:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description='Generate Gadget ICs of a polytropic profile.',
                                              formatter_class=argparse.RawTextHelpFormatter)
        self.parser.add_argument("-g", "--gamma",
                            dest     = "gamma",
                            type     = float,
                            help     = "Polytropic index",
                            required = True)

        self.parser.add_argument("-N", "-n",
                            dest     = "num",
                            type     = int,
                            help     = "Number of gaseous particles",
                            required = True)

        self.parser.add_argument("-o",
                            metavar = "outfile",
                            dest    = "outfile",
                            help    = "Name of output file",
                            default = "Polytrope.dat")

    def get_args(self):
        return self.parser.parse_args()
