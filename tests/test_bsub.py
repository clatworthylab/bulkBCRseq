#!/usr/bin/env python
# basic requirements for test data
import os
from subprocess import run
from typing import Union


def test_print_bsub_cmd(step: int, metadata: Union[None, str] = None, bsub: bool = True, print: bool = True, run: bool = False):
    env = os.environ.copy()
    run(cmd, env=env)

    if metadata is None:
        meta = 'tests/data/Sample_metadata.txt'
    else:
        meta = metadata

    cmd = ['python',
           'Processing_sequences_large_scale.py',
           meta,
           str(step),
           'Y',
           'Y',
           'N']
    run(cmd, env=env)


if __name__ == "__main__":
    test_print_bsub_cmd(1)
    test_print_bsub_cmd(2)
    test_print_bsub_cmd(3)
    test_print_bsub_cmd(4)
