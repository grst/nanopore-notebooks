#!/usr/bin/env python
# -*- coding: utf-8 -*-
from nbwrapper import nbwrapper
import argparse


argp = argparse.ArgumentParser()
argp.add_argument("--foo", '-a', help="foo arg")
argp.add_argument("--bar", '-b', help="bar arg")

args = argp.parse_args()
nb = nbwrapper(args, "test.ipynb")
nb.run()

