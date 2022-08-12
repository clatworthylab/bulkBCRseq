[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5717959.svg)](https://doi.org/10.5281/zenodo.5717959)
[![codecov](https://codecov.io/gh/clatworthylab/bulkBCRseq/branch/master/graph/badge.svg?token=I6APMCARTA)](https://codecov.io/gh/clatworthylab/bulkBCRseq)

# bulk_BCR_analysis
Bulk BCR-seq processing scripts use in Fitzpatrick et al., Nature (2020). Package belongs to Rachael Bashford-Rogers.

This repo is an older version of what seems to be now at https://github.com/rbr1/BCR_TCR_PROCESSING_PIPELINE.

Requires python>=3.8 (or 2.7 if using the legacy branch). Currently only works when cloned onto farm with all paths set up pointing to this folder properly.

## Citation
Please cite the following papers:

*Fitzpatrick, Z., Frazer, G., Ferro, A., Clare, S., Bouladoux, N., Ferdinand, J., Tuong, Z.K., Negro-Demontel, M.L., Kumar, N., Suchanek, O. and Tajsic, T., 2020. Gut-educated IgA plasma cells defend the meningeal venous sinuses. Nature, 587(7834), pp.472-476.*

*Bashford-Rogers, R.J., Palser, A.L., Huntly, B.J., Rance, R., Vassiliou, G.S., Follows, G.A. and Kellam, P., 2013. Network properties derived from deep sequencing of human B-cell receptor repertoires delineate B-cell populations. Genome research, 23(11), pp.1874-1884.*

*Bashford-Rogers, R.J.M., Bergamaschi, L., McKinney, E.F., Pombal, D.C., Mescia, F., Lee, J.C., Thomas, D.C., Flint, S.M., Kellam, P., Jayne, D.R.W. and Lyons, P.A., 2019. Analysis of the B cell receptor repertoire in six immune-mediated diseases. Nature, 574(7776), pp.122-126.*


## Pre-requisites:
```bash
# create a conda virtual environment
# sample for python 3 set up, switch to python 2 where appropriate
# install miniconda
# see https://docs.conda.io/en/latest/miniconda.html#linux-installers
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh
eval "$(/path/to/miniconda2/bin/conda shell.bash hook)"
conda init
conda create --name isotyper python=3.9

# clone this repo
git clone https://github.com/clatworthylab/bulkBCRseq

# change into the directory and install dependencies
cd bulkBCRseq
conda env update --name isotyper --file environment.yml
```

Usage instructions on Farm:
```bash
conda activate isotyper
```

## Note!
If you are starting from fastq files directly, please change the 2nd column in the `.txt` file (path to `.cram`) to path to `_R1_001.fastq.gz` (read1) instead. If your read1 suffix isn't this pattern, please modify the `R1PATTERN` variable after cloning this repository, in here directly:
https://github.com/clatworthylab/bulkBCRseq/blob/3d17a2752a6b482f50c0b8d211db94ddf5e655d1/BIN/Read_processing_and_quality.py#L3641-L3643


## Basic usage:
```
usage: isotyper.py [-h] [-i INPUT] [-s STEP] [-c CORES] [-l LENGTH] [-dr] [-b] [-m MEM] [-q QUEUE] [-p PROJECT] [-g GROUP]

optional arguments:
  -h, --help            show this help message and exit

main arguments:
  -i INPUT, --input INPUT
                        input meta.txt file to run isotyper.
                        file must contain the following four columns:
                            1st column - name of sample.
                            2nd column - path to input file. Either .cram file or read 1 fastq(.gz) file.
                            3rd column - path to output folder.
                            4th column - organism. Either HOMO_SAPIENS or MUS_MUSCULUS.
                            no column names allowed.
  -s STEP, --step STEP  step to perform:
                            1 - Convert raw sequencing files to fastq and perform QC.
                            2 - Trim and filter reads.
                            3 - Generate networks.
                            4 - Generate network statistics.
  -c CORES, --cores CORES
                        number of cores to run this on. [Default 1]
  -l LENGTH, --length LENGTH
                        minimum length of reads to keep. [Default 100]
  -dr, --dryrun         if passed, prints commands but don't actually run.

bsub arguments:
  -b, --bsub            if passed, submits each row in meta.txt file as a job to bsub.
  -m MEM, --mem MEM     job memory request. [Default 8000]
  -q QUEUE, --queue QUEUE
                        job queue to submit to. [Default normal]
  -p PROJECT, --project PROJECT
                        sanger project to send as job. [Default team205]
  -g GROUP, --group GROUP
                        sanger group to send as job. [Default teichlab]
```

### Basic usage
```bash
# initial QC
python isotyper.py -i meta.txt -s 1
# trimming
python isotyper.py -i meta.txt -s 2  --cores 10
# generate network
python isotyper.py -i meta.txt -s 3
# generate network statistic
python isotyper.py -i meta.txt -s 4
```

If using Sanger's farm:
```bash
# initial QC
python isotyper.py -i meta.txt -s 1 --bsub
# trimming
python isotyper.py -i meta.txt -s 2 --bsub --cores 10
# generate network
python isotyper.py -i meta.txt -s 3 --bsub
# generate network statistic
python isotyper.py -i meta.txt -s 4 --bsub
```

Take a look [here](https://github.com/clatworthylab/bulkBCRseq/tree/master/tests/data) for example files to provide to the tool.