#!/usr/bin/env python
import os
from subprocess import call

import pytest


class Tester:

    """Summary

    Attributes
    ----------
    bsub : TYPE
        Description
    execute : TYPE
        Description
    metadata : TYPE
        Description
    option : TYPE
        Description
    verbose : TYPE
        Description
    """

    __test__ = False

    def __init__(self, option, metadata, bsub, verbose, execute):
        """Summary

        Parameters
        ----------
        option : TYPE
            Description
        metadata : TYPE
            Description
        bsub : TYPE
            Description
        verbose : TYPE
            Description
        execute : TYPE
            Description
        """
        self.option = option
        self.metadata = metadata
        self.bsub = bsub
        self.verbose = verbose
        self.execute = execute

    def run_test(self):
        """Summary"""
        env = os.environ.copy()
        if self.option is None:
            opt = 1
        else:
            opt = self.option
        if self.metadata is None:
            meta = "tests/data/Sample_metadata.txt"
        else:
            meta = self.metadata
        if self.bsub:
            bsub_ = "Y"
        else:
            bsub_ = "N"
        if self.verbose:
            verbose_ = "Y"
        else:
            verbose_ = "N"
        if self.execute:
            execute_ = "Y"
        else:
            execute_ = "N"
        cmd = [
            "python",
            "Processing_sequences_large_scale.py",
            meta,
            str(opt),
            bsub_,
            verbose_,
            execute_,
        ]
        call(cmd, env=env)


@pytest.fixture
def tester_standard(request):
    """Summary

    Parameters
    ----------
    request : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    return Tester(request.param, None, False, True, False)


@pytest.fixture
def tester_bsub(request):
    """Summary

    Parameters
    ----------
    request : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    return Tester(request.param, None, True, True, False)


class TestRun:

    """Summary"""

    @pytest.mark.parametrize(
        "tester_standard", [1, 2, 3, 4], indirect=["tester_standard"]
    )
    def test_run1(self, tester_standard):
        """Summary

        Parameters
        ----------
        tester_standard : TYPE
            Description
        """
        tester_standard.run_test()

    @pytest.mark.parametrize(
        "tester_bsub", [1, 2, 3, 4], indirect=["tester_bsub"]
    )
    def test_run2(self, tester_bsub):
        """Summary

        Parameters
        ----------
        tester_bsub : TYPE
            Description
        """
        tester_bsub.run_test()
