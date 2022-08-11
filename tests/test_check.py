#!/usr/bin/env python
import pytest
import shutil

from pathlib import Path
from subprocess import run
from glob import glob

TESTFOLDER = Path("tests")
TESTDATFOLDER = TESTFOLDER / "data"
TESTOUTFOLDER = TESTFOLDER / "output"


@pytest.mark.parametrize(
    "option,expected2,expected3,expected4",
    [
        pytest.param(1, 0, 0, 0),
        pytest.param(2, 1, 0, 0),
        pytest.param(3, 1, 1, 0),
        # 4
    ],
)
def test_call_script(option, expected2, expected3, expected4):
    """Basic call script."""
    cmd = [
        "python",
        "Processing_sequences_large_scale.py",
        str(TESTDATFOLDER / "Sample_metadata.txt"),
        str(option),
        "N",
        "Y",
        "Y",
    ]
    run(cmd)
    assert len(glob(str(TESTOUTFOLDER / "FASTQ_FILES"))) == 1
    assert len(glob(str(TESTOUTFOLDER / "ORIENTATED_SEQUENCES"))) == 1
    assert len(glob(str(TESTOUTFOLDER / "FASTQ_FILES" / "*.qc.fq"))) == 2
    assert (
        len(glob(str(TESTOUTFOLDER / "FASTQ_FILES" / "Fail*.fasta")))
        == expected2
    )
    assert (
        len(glob(str(TESTOUTFOLDER / "NETWORKS" / "Fully_reduced*")))
        == expected3
    )


def test_clean_up():
    for out in ["FASTQ_FILES", "ORIENTATED_SEQUENCES"]:
        shutil.rmtree(TESTOUTFOLDER / out)


# @pytest.mark.parametrize("option", [1, 2, 3, 4])
# def test_call_script_fastq_gz(option):
#     """Basic call script on fastq.gz."""
#     cmd = [
#         "python",
#         "Processing_sequences_large_scale.py",
#         str(TESTDATFOLDER / "Sample_metadata2a.txt"),
#         str(option),
#         "N",
#         "Y",
#         "Y",
#     ]
#     run(cmd)
#     assert len(glob(str(TESTOUTFOLDER / "FASTQ_FILES" / "*test1*.qc.fq"))) > 0


# def test_prep_fastq():
#     """unzip gz file."""
#     for file in ["test1_R1_001.fastq.gz", "test1_R2_001.fastq.gz"]:
#         zipped = TESTDATFOLDER / file
#         run(["gunzip", zipped])


# @pytest.mark.parametrize("option", [1, 2, 3, 4])
# def test_call_script_fastq(option):
#     """Basic call script on fastq."""
#     cmd = [
#         "python",
#         "Processing_sequences_large_scale.py",
#         str(TESTDATFOLDER / "Sample_metadata2b.txt"),
#         str(option),
#         "N",
#         "Y",
#         "Y",
#     ]
#     run(cmd)
