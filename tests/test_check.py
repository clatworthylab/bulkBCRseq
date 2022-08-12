#!/usr/bin/env python
import pytest
import shutil

from pathlib import Path
from subprocess import run
from glob import glob

TESTFOLDER = Path("tests")
TESTDATFOLDER = TESTFOLDER / "data"
TESTOUTFOLDER = TESTFOLDER / "output"
TESTORTSEQFOLDER = TESTOUTFOLDER / "ORIENTATED_SEQUENCES"
TESTFASTQFOLDER = TESTOUTFOLDER / "FASTQ_FILES"
TESTNET = TESTORTSEQFOLDER / "NETWORKS"


@pytest.mark.parametrize(
    "option,expected2,expected3,expected4",
    [
        pytest.param(1, 0, 0, 0),
        pytest.param(2, 1, 0, 0),
        pytest.param(3, 1, 1, 0),
        pytest.param(4, 1, 1, 2),
    ],
)
def test_call_script(option, expected2, expected3, expected4):
    """Basic call script."""
    cmd = [
        "python",
        "isotyper.py",
        "-i",
        str(TESTDATFOLDER / "Sample_metadata1.txt"),
        "-s",
        str(option),
    ]
    run(cmd)
    assert len(glob(str(TESTFASTQFOLDER))) == 1
    assert len(glob(str(TESTORTSEQFOLDER))) == 1
    assert len(glob(str(TESTFASTQFOLDER / "*.qc.fq"))) == 2
    assert len(glob(str(TESTFASTQFOLDER / "Fail*.fasta"))) == expected2
    assert len(glob(str(TESTNET / "Fully_reduced*.fasta"))) == expected3
    assert (
        len(glob(str(TESTOUTFOLDER / "Network_statistics*.txt"))) == expected4
    )


def test_clean_up():
    """clean up after test run."""
    for out in [TESTORTSEQFOLDER, TESTFASTQFOLDER]:
        shutil.rmtree(out)
