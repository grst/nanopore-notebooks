#!/usr/bin/env python
# -*- coding: utf-8 -*-
ENV_NAME = "ipythonargs" #name of environment variable

import warnings

with warnings.catch_warnings():
    # the current version of runipy uses deprecated methods.
    # supress warnings until new version of runipy is available.
    warnings.simplefilter("ignore")
    from runipy.notebook_runner import NotebookRunner
    from IPython.nbformat.current import read

from os import environ
import json

def getargs():
    """
    Retrieve the arguments within ipython notebook.

    Returns:
        dict containing the arguments.
    """

    try:
        args = json.loads(environ[ENV_NAME])
    except KeyError:
        warnings.warn("no arguments passed!", RuntimeWarning)
        args = {}
    return args

class nbwrapper:
    """ run a ipython notebook with arguments """
    def __init__(self, args, notebookpath):
        """
        Args:
            args: Namespace containing the arguments to pass to the notebook
            notebookpath: path to the ipynb file
        """
        self._args = vars(args)
        self._notebook = read(open(notebookpath), 'json')

    def run(self):
        environ[ENV_NAME] = json.dumps(self._args)
        r = NotebookRunner(self._notebook)
        r.run_notebook()
        for i, cell in enumerate(r.nb["worksheets"][0]["cells"]):
            print(">IN [{0}]:".format(i))
            print (cell["input"])
            print()
            print(">OUT [{0}]:".format(i))
            try:
                print(cell["outputs"][0]["text"])
            except IndexError:
                print()
            print("------------------")



