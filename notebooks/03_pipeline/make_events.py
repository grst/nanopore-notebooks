#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, "../../nbwrapper")
from nbwrapper import Nbwrapper
import argparse


if __name__ == "__main__":
    argp = argparse.ArgumentParser()
    argp.add_argument("--f5_path", nargs="?", help="path to folder containing processed fast5-files")
    argp.add_argument("--output", nargs="?", help="path to output file")
    argp.add_argument("--ncores", nargs="?", default="4", help="number of CPU cores to use")

    args =  argp.parse_args()
    print(args)
    nb = Nbwrapper(args, "make_events.ipynb")
    nb.run()


