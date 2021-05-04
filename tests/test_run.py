#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2021-04-29 20:57:15
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-05-04 11:07:18

import os
from subprocess import call

import pytest


class Tester:
    __test__ = False

    def __init__(self, option, metadata, bsub, verbose, execute, base, sub):
        self.option = option
        self.metadata = metadata
        self.bsub = bsub
        self.verbose = verbose
        self.execute = execute
        self.base = base
        self.sub = sub

    def run_test(self):
        env = os.environ.copy()
        if self.option is None:
            opt = 1
        else:
            opt = self.option
        if self.metadata is None:
            meta = 'tests/data/Sample_metadata.txt'
        else:
            meta = self.metadata
        if self.bsub:
            bsub_ = 'Y'
        else:
            bsub_ = 'N'
        if self.verbose:
            verbose_ = 'Y'
        else:
            verbose_ = 'N'
        if self.execute:
            execute_ = 'Y'
        else:
            execute_ = 'N'
        cmd1 = [
            'python', 'Processing_sequences_large_scale.py', meta,
            str(opt), bsub_, verbose_, execute_
        ]
        call(cmd1, env=env)

    def check_output(self):
        for f in self.sub:
            cmd2 = ['ls', self.base + f]
            call(cmd2)


@pytest.fixture
def tester(request):
    return Tester(request.param, None, False, True, True, 'tests/output/', [
        'FASTQ_FILES', 'ORIENTATED_SEQUENCES', 'ORIENTATED_SEQUENCES/NETWORKS',
        'ORIENTATED_SEQUENCES/TMP'
    ])


class TestRun:
    @pytest.mark.parametrize('tester', [1, 2, 3, 4], indirect=['tester'])
    def test_run(self, tester):
        tester.run_test()
        tester.check_output()
