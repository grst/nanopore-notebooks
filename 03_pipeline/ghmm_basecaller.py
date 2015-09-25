#!/home/ibis/gregor.sturm/bin/anaconda/bin/python2.7
# -*- coding: utf-8 -*-

import sys
from nbwrapper import Nbwrapper
import argparse
import os

if __name__ == "__main__":
    argp = argparse.ArgumentParser()
    argp.add_argument("--events", nargs="?", help="event-pickle")
    argp.add_argument("--out_basename", nargs="?", help="path to output file")
    argp.add_argument("--ncores", nargs="?", default="4", help="number of CPU cores to use")
    argp.add_argument("--ref", nargs="?")
    argp.add_argument("--hmm_params", nargs="?")
    argp.add_argument("--nmers", nargs="?")
    argp.add_argument("--multivariate", nargs="?")

    args =  argp.parse_args()
    print(args)
    path = os.path.dirname(os.path.realpath(__file__))
    nb = Nbwrapper(args, path + "/ghmm_basecaller_multivariate.ipynb")
    nb.run()


