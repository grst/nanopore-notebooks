#!/usr/bin/env python
# -*- coding: utf-8 -*-
ENV_NAME = "ipythonargs" #name of environment variable

import warnings

with warnings.catch_warnings():
    # the current version of runipy uses deprecated methods.
    # supress warnings until new version of runipy is available.
    warnings.simplefilter("ignore")
    from runipy.notebook_runner import NotebookRunner, NotebookError
    from IPython.nbformat.current import read

from os import environ
import json
import logging
import atexit


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

class CustomNotebookRunner(NotebookRunner):
    """
    extends the NotebookRunner by a cell-specific-callback
    """
    def run_notebook(self, skip_exceptions=False,
            progress_callback=None, cell_callback=None):
        '''
        Run all the cells of a notebook in order and update
        the outputs in-place.
        If ``skip_exceptions`` is set, then if exceptions occur in a cell, the
        subsequent cells are run (by default, the notebook execution stops).
        '''
        for i, cell in enumerate(self.iter_code_cells()):
            try:
                self.run_cell(cell)
            except NotebookError:
                if not skip_exceptions:
                    raise
            if progress_callback:
                progress_callback(i)
            if cell_callback:
                cell_callback(cell)


class Nbwrapper:
    """ run a ipython notebook with arguments """
    def __init__(self, args, notebookpath):
        """
        Args:
            args: Namespace containing the arguments to pass to the notebook
            notebookpath: path to the ipynb file
        """
        log_format = '%(asctime)s %(levelname)s: %(message)s'
        log_datefmt = '%m/%d/%Y %I:%M:%S %p'
        logging.basicConfig(level=logging.INFO, format=log_format, datefmt=log_datefmt)
        self._args = vars(args)
        self._notebook = read(open(notebookpath), 'json')

    def run(self):
        def cell_callback(cell):
            try:
                out = "OUTPUT:\n" + cell["outputs"][0]["text"]
            except IndexError:
                out = "\n"
            logging.info(out)

        def cleanup(r):
            r.shutdown_kernel()

        environ[ENV_NAME] = json.dumps(self._args)
        r = CustomNotebookRunner(self._notebook)
        atexit.register(cleanup, r)
        r.run_notebook(cell_callback=cell_callback)
        cleanup(r)



