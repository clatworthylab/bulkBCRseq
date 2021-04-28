#!/usr/bin/env python
# basic requirements for test data
import os
from subprocess import run
from typing import Union


def test_print_bsub_cmd(option, metadata = None, bsub = True, verbose = True, execute = False):
    """
    Test to print the bsub commands

    Parameters
    ----------
    option : int
        options for running the main script
    metadata : str, optional
        path to metadata.txt file
    bsub : bool
        whether or not use bsub command
    verbose : bool
        whether or not to print command
    execute : bool
        where or not to run command

    """
    env = os.environ.copy()
    run(cmd, env=env)

    if metadata is None:
        meta = 'tests/data/Sample_metadata.txt'
    else:
        meta = metadata

    if bsub:
      bsub_ = 'Y'
    else:
      bsub_ = 'N'

    if verbose:
      verbose_ = 'Y'
    else:
      verbose_ = 'N'

    if execute:
      execute_ = 'Y'
    else:
      execute_ = 'N'

    cmd = ['python',
           'Processing_sequences_large_scale.py',
           meta,
           str(option),
           bsub_,
           verbose_,
           execute_]
    run(cmd, env=env)


if __name__ == "__main__":
    test_print_bsub_cmd(1)
    test_print_bsub_cmd(2)
    test_print_bsub_cmd(3)
    test_print_bsub_cmd(4)
