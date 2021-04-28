#!/usr/bin/env python
import os
import pytest
from subprocess import call  # it's equivalent to run in python>=3.5


@pytest.fixture()
def argument_printer():
    def _test_script(option=None, metadata=None, bsub=True, verbose=True, execute=False):
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

        if option is None:
            opt = 1
        else:
            opt = option

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
               str(opt),
               bsub_,
               verbose_,
               execute_]
        call(cmd, env=env)

    return _test_script

def test_script_1_bsub(argument_printer):
    argument_printer(1, None, True, True, False)

def test_script_2_bsub(argument_printer):
    argument_printer(2, None, True, True, False)

def test_script_3_bsub(argument_printer):
    argument_printer(3, None, True, True, False)

def test_script_4_bsub(argument_printer):
    argument_printer(4, None, True, True, False)

def test_script_1(argument_printer):
    argument_printer(1, None, False, True, False)

def test_script_2(argument_printer):
    argument_printer(2, None, False, True, False)

def test_script_3(argument_printer):
    argument_printer(3, None, False, True, False)

def test_script_4(argument_printer):
    argument_printer(4, None, False, True, False)

if __name__ == "__main__":
    # print bsub commands
    test_script_1_bsub()
    test_script_2_bsub()
    test_script_3_bsub()
    test_script_4_bsub()
    # print non-bsub commands
    test_script_1()
    test_script_2()
    test_script_3()
    test_script_4()
