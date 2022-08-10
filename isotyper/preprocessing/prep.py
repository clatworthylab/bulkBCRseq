#!/usr/bin/env python
import re
import shutil
import subprocess

from glob import glob
from pathlib import Path
from typing import Union

from isotyper.utilities._settings import R1PATTERN, R2PATTERN, MIN_QUAL, EXTPATH
from isotyper.utilities._args import (
    LENGTH,
    OUT_FASTQ,
    OUT_NET,
    OUT_ORTSEQ,
    OUT_ORTSEQ_TMP,
    OUT_PATH,
    SAMPLE_ID,
    SOURCE,
)


def bam_to_fastq(source_path: Path):
    """Convert bam to fastq,

    Parameters
    ----------
    source_path : Path
        path of input file (cram or bam).
    """
    if (
        len(
            glob(str(source_path.parent / "*.bam"))
            + glob(str(source_path.parent / "*.cram"))
        )
        != 0
    ):
        pre_qc_fastq1 = OUT_FASTQ / f"Sequences_{SAMPLE_ID}_1.fastq"
        pre_qc_fastq2 = OUT_FASTQ / f"Sequences_{SAMPLE_ID}_2.fastq"
        pre_qc_bam = OUT_FASTQ / f"Sequences_{SAMPLE_ID}.bam"
        if len(glob(source_path.parent / "*.cram")) != 0:
            cram_to_bam(cram_path=source_path, out_pre_qc_bam_path=pre_qc_bam)
        command1 = [
            "picard",
            "SamToFastq",
            f"I={str(pre_qc_bam)}",
            f"FASTQ={pre_qc_fastq1}",
            f"SECOND_END_FASTQ={pre_qc_fastq2}",
        ]
        subprocess.run(command1)


def copy_prepped_fastq(fastq_path: Path, read_num: str):
    """Copy and unzip fastq files.

    Parameters
    ----------
    fastq_path : Path
        location of original fastq file.
    read_num : str
        `1` or `2` to specify read 1 or read 2.
    """
    if fastq_path.suffix == ".gz":
        extension = ".gz"
    else:
        extension = ""
    new_fastq = OUT_FASTQ / f"Sequences_{SAMPLE_ID}_{read_num}.fastq{extension}"
    shutil.copy(fastq_path, new_fastq)
    if fastq_path.suffix == ".gz":
        cmdg = ["gunzip", "-f", f"{new_fastq}"]
        subprocess.run(cmdg)


def cram_to_bam(cram_path: Path, out_pre_qc_bam_path: Path):
    """Convert cram to bam.

    Parameters
    ----------
    cram_path : Path
        path to cram file.
    out_pre_qc_bam_path : Path
        path to output bam file.
    """
    cmd1 = [
        "samtools",
        "view",
        "-b",
        "-o",
        str(out_pre_qc_bam_path),
        str(cram_path),
    ]
    print(cmd)
    subprocess.run(cmd)


def intialise_files():
    """Initialise to output folder."""
    for d in [
        OUT_FASTQ,
        OUT_ORTSEQ,
        OUT_ORTSEQ_TMP,
        OUT_NET,
    ]:
        d.mkdir(exist_ok=True, parents=True)


def prep_fastqs(source_path: Path, r1pattern: str, r2pattern: str):
    """Prepare fastqs for input into the script.

    Parameters
    ----------
    source_path : Path
        location of R1 of fastq file.
    r1pattern : str
        suffix pattern before .fastq to try and match for R1
    r2pattern : str
        suffix pattern before .fastq to try and match for R2
    """
    if re.search(r1pattern, str(source_path)):
        r1_original = Path(source_path)
        r2_original = Path(re.sub(r1pattern, r2pattern, str(source_path)))
    else:
        raise ValueError(
            "Your input file {} does not contain the {} pattern.".format(
                str(source_path), r1pattern
            )
        )
    copy_prepped_fastq(
        fastq_path=r1_original,
        read_num="1",
    )
    copy_prepped_fastq(
        fastq_path=r2_original,
        read_num="2",
    )


def qc_samples(
    out_path: Path,
    min_length: Union[str, int] = 100,
    min_threshold: Union[str, int] = 32,
):
    """Perform QC on samples using QUASR and convert to fasta file with perl script.

    Parameters
    ----------
    out_path : Path
        path to output folder.
    min_length : Union[str, int], optional
        minimum read length cutoff.
    min_threshold : Union[str, int], optional
        minimum median-read-quality cutoff.
    """
    reads1 = out_path / f"Sequences_{SAMPLE_ID}_1.fastq"
    reads2 = out_path / f"Sequences_{SAMPLE_ID}_2.fastq"
    quasr_qc_jar_path = EXTPATH / "QUASR_v7.01" / "qualityControl.jar"
    # see https://github.com/andrewjpage/QUASR for updated version
    cmd1 = [
        "java",
        "-jar",
        str(quasr_qc_jar_path),
        "-f",
        str(reads1),
        "-o",
        str(reads1.parent / reads1.stem),
        "-m",
        str(min_threshold),
        "-l",
        str(min_length),
    ]
    cmd2 = [
        "java",
        "-jar",
        str(quasr_qc_jar_path),
        "-f",
        str(reads2),
        "-o",
        str(reads2.parent / reads2.stem),
        "-m",
        str(min_threshold),
        "-l",
        str(min_length),
    ]
    cmd3 = [
        "perl",
        "-e",
        PERLCMD,
    ]
    subprocess.run(cmd1)
    subprocess.run(cmd2)
    subprocess.run(
        cmd3,
        stdin=open(reads1.with_suffix(".qc.fq"), "r"),
        stdout=open(reads1.with_suffix(".fasta"), "w"),
    )
    subprocess.run(
        cmd3,
        stdin=open(reads2.with_suffix(".qc.fq"), "r"),
        stdout=open(reads2.with_suffix(".fasta"), "w"),
    )


def main():
    """main function to prepare step 1."""
    intialise_files()
    if (
        len(
            glob(str(SOURCE.parent / "*.bam*"))
            + glob(str(SOURCE.parent / "*.cram"))
        )
        != 0
    ):
        bam_to_fastq(source_path=SOURCE)
    elif (
        len(
            glob(str(SOURCE.parent / "*.fastq"))
            + glob(str(SOURCE.parent / "*.fastq.gz"))
            + glob(str(SOURCE.parent / "*.fq*"))
            + glob(str(SOURCE.parent / "*.fq.gz"))
        )
        != 0
    ):
        # rename them to Seqeuence_{SAMPLE_ID}_1.fastq Seqeuence_{SAMPLE_ID}_2.fastq
        prep_fastqs(
            source_path=SOURCE,
            r1pattern=R1PATTERN,
            r2pattern=R2PATTERN,
        )
    qc_samples(
        out_path=OUT_FASTQ,
        min_length=LENGTH,
        min_threshold=MIN_QUAL,
    )


if __name__ == "__main__":
    main()
