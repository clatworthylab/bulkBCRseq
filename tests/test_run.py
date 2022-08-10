#!/usr/bin/env python
import os
from subprocess import call

import pytest


class Tester:

    """Summary

    Attributes
    ----------
    base : TYPE
        Description
    bsub : TYPE
        Description
    execute : TYPE
        Description
    metadata : TYPE
        Description
    option : TYPE
        Description
    sub : TYPE
        Description
    verbose : TYPE
        Description
    """

    __test__ = False

    def __init__(self, option, metadata, bsub, verbose, execute, base, sub):
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
        base : TYPE
            Description
        sub : TYPE
            Description
        """
        self.option = option
        self.metadata = metadata
        self.bsub = bsub
        self.verbose = verbose
        self.execute = execute
        self.base = base
        self.sub = sub

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
        cmd1 = [
            "python",
            "Processing_sequences_large_scale.py",
            meta,
            str(opt),
            bsub_,
            verbose_,
            execute_,
        ]
        call(cmd1, env=env)

    def check_output(self):
        """Summary"""
        for f in self.sub:
            cmd2 = ["ls", self.base + f]
            call(cmd2)


@pytest.fixture
def tester(request):
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
    return Tester(
        request.param,
        None,
        False,
        True,
        True,
        "tests/output/",
        [
            "FASTQ_FILES",
            "ORIENTATED_SEQUENCES",
            "ORIENTATED_SEQUENCES/NETWORKS",
            "ORIENTATED_SEQUENCES/TMP",
        ],
    )


class TestRun:

    """Summary"""

    @pytest.mark.parametrize("tester", [1, 2, 3, 4], indirect=["tester"])
    def test_run(self, tester):
        """Summary

        Parameters
        ----------
        tester : TYPE
            Description
        """
        tester.run_test()
        tester.check_output()
