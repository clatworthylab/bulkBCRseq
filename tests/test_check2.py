#!/usr/bin/env python
from subprocess import run

import pytest


@pytest.mark.parametrize("option", [1, 2, 3, 4])
def test_call_script(option):
    """Basic call script."""
    cmd = [
        "python",
        "Processing_sequences_large_scale.py",
        "tests/data/Sample_metadata.txt",
        str(option),
        "N",
        "Y",
        "Y",
    ]
    run(cmd)
