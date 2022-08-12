#!/usr/bin/env python
import subprocess
import textwrap
import pandas as pd
from pathlib import Path
from typing import Literal, Tuple, List
import argparse

from isotyper.utilities._settings import ISOTYPERPREFIX

Path("logs").mkdir(exist_ok=True)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter
    )
    # main arguments
    main_group = parser.add_argument_group(title="main arguments")
    main_group.add_argument(
        "-i",
        "--input",
        help=textwrap.dedent(
            """\
            input meta.txt file to run isotyper.
            file must contain the following four columns:
                1st column - name of sample.
                2nd column - path to input file. Either .cram file or read 1 fastq(.gz) file.
                3rd column - path to output folder.
                4th column - organism. Either HOMO_SAPIENS or MUS_MUSCULUS.
                no column names allowed.
            """
        ),
    )
    main_group.add_argument(
        "-s",
        "--step",
        help=textwrap.dedent(
            """\
            step to perform:
                1 - Convert raw sequencing files to fastq and perform QC.
                2 - Trim and filter reads.
                3 - Generate networks.
                4 - Generate network statistics.
            """
        ),
    )
    main_group.add_argument(
        "-l",
        "--length",
        default=100,
        help=("minimum length of reads to keep. [Default 100]"),
    )
    main_group.add_argument(
        "-dr",
        "--dryrun",
        action="store_true",
        help=("if passed, prints commands but don't actually run."),
    )
    # bsub arguments
    bsub_group = parser.add_argument_group(title="bsub arguments")
    bsub_group.add_argument(
        "-b",
        "--bsub",
        action="store_true",
        help=("if passed, submits each row in meta.txt file as a job to bsub."),
    )
    bsub_group.add_argument(
        "-c",
        "--cores",
        default=10,
        help=("number of cores to run this on. [Default 10]"),
    )
    bsub_group.add_argument(
        "-m", "--mem", default=8000, help=("job memory request. [Default 8000]")
    )
    bsub_group.add_argument(
        "-q",
        "--queue",
        default="normal",
        help=("job queue to submit to. [Default normal]"),
    )
    bsub_group.add_argument(
        "-p",
        "--project",
        default="team205",
        help=("sanger project to send as job. [Default team205]"),
    )
    bsub_group.add_argument(
        "-g",
        "--group",
        default="teichlab",
        help=("sanger group to send as job. [Default teichlab]"),
    )
    args = parser.parse_args()
    return args


def get_info(file: Path) -> Tuple[List[str], List[Path], List[Path], List[str]]:
    """Summary

    Parameters
    ----------
    file : Path
        path to input meta.txt file.

    Returns
    -------
    Tuple[List[str], List[Path], List[Path], List[str]]
        lists of sample ids, input paths, output paths and organism/species.
    """
    sep = "\t" if file.suffix != ".csv" else ","
    meta = pd.read_csv(file, sep=sep, header=None)
    sample_ids = list(meta[0])
    input_paths = [Path(i) for i in meta[1]]
    output_paths = [Path(o) for o in meta[2]]
    orgs = list(meta[3])
    return (
        sample_ids,
        input_paths,
        output_paths,
        orgs,
    )


def bsub_job_options(
    sample_id: str,
    mem: int = 8000,
    queue: Literal["normal", "long"] = "normal",
    ncores: int = 1,
    project: str = "team205",
    group: str = "teichlab",
) -> List:
    """Generate bsub commands.

    Parameters
    ----------
    sample_id : str
        name of sample to run.
    mem : int, optional
        memory required.
    queue : Literal["normal", "long"], optional
        job queue.
    ncores : int, optional
        number of cores
    project : str, optional
        project name on farm
    group : str, optional
        group name on farm
    """
    prog_arg = ["-P", project]
    grp_arg = ["-G", group]
    log_arg = ["-o", f"logs/log_{sample_id}"]
    core_arg = ["-n", str(ncores)]
    mem_arg = [
        f"-R'select[mem>{str(mem)}] rusage[mem={str(mem)}] span[hosts=1]'",
        "-M",
        f"{str(mem)}",
    ]
    queue_arg = ["-q", queue]
    cmd = (
        ["bsub"] + prog_arg + grp_arg + log_arg + queue_arg + core_arg + mem_arg
    )
    return cmd


def main():
    """Main script."""
    args = parse_args()
    sample_ids, input_paths, output_paths, orgs = get_info(
        file=Path(args.input)
    )
    for i in range(0, len(sample_ids)):
        if args.bsub:
            job_arg = bsub_job_options(
                sample_id=sample_ids[i],
                mem=args.mem,
                queue=args.queue,
                ncores=args.cores,
                project=args.project,
                group=args.group,
            )
        else:
            job_arg = []
        main_cmd = job_arg + [
            "python",
            str(ISOTYPERPREFIX / "read_processing_and_quality.py"),
            str(sample_ids[i]),
            str(input_paths[i]),
            str(output_paths[i]),
            str(orgs[i]),
            str(args.length),
            str(args.step),
        ]
        print(" ".join(main_cmd))
        if not args.dryrun:
            subprocess.run(main_cmd)


if __name__ == "__main__":
    main()
