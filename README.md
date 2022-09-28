[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5717959.svg)](https://doi.org/10.5281/zenodo.5717959)
[![codecov](https://codecov.io/gh/clatworthylab/bulkBCRseq/branch/master/graph/badge.svg?token=I6APMCARTA)](https://codecov.io/gh/clatworthylab/bulkBCRseq)

# bulkBCRseq : isotyper
Bulk BCR-seq processing package used in `Fitzpatrick et al., Nature (2020)`. The original (legacy) package/scripts was provided by Dr. Rachael Bashford-Rogers (Oxford).

This repository is a `python3` reimplementation of the original `python2` scripts (found in [legacy branch](https://github.com/clatworthylab/bulkBCRseq/tree/legacy)); the original script is an older version of what seems to be now at https://github.com/rbr1/BCR_TCR_PROCESSING_PIPELINE.

Requires `python>=3.8` (or `python==2.7.9` if using the [legacy branch](https://github.com/clatworthylab/bulkBCRseq/tree/legacy)).

## Citation
Please cite the following papers:
```
Fitzpatrick, Z., Frazer, G., Ferro, A., Clare, S., Bouladoux, N., Ferdinand, J., Tuong, Z.K., Negro-Demontel, M.L., Kumar, N., Suchanek, O. and Tajsic, T., 2020. Gut-educated IgA plasma cells defend the meningeal venous sinuses. Nature, 587(7834), pp.472-476.

Bashford-Rogers, R.J., Palser, A.L., Huntly, B.J., Rance, R., Vassiliou, G.S., Follows, G.A. and Kellam, P., 2013. Network properties derived from deep sequencing of human B-cell receptor repertoires delineate B-cell populations. Genome research, 23(11), pp.1874-1884.

Bashford-Rogers, R.J.M., Bergamaschi, L., McKinney, E.F., Pombal, D.C., Mescia, F., Lee, J.C., Thomas, D.C., Flint, S.M., Kellam, P., Jayne, D.R.W. and Lyons, P.A., 2019. Analysis of the B cell receptor repertoire in six immune-mediated diseases. Nature, 574(7776), pp.122-126.
```

## Setup:
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

# clone this repository
git clone https://github.com/clatworthylab/bulkBCRseq

# change into the directory and install dependencies
cd bulkBCRseq
conda env update --name isotyper --file environment.yml
```

```bash
# either run this everytime or just
# export to your ~/.bashrc or ~/.bash_profile
export PYTHONPATH=/path/to/bulkBCRseq:$PYTHONPATH
export REF_PATH=/lustre/scratch117/core/sciops_repository/cram_cache/%2s/%2s/%s:/lustre/scratch118/core/sciops_repository/cram_cache/%2s/%2s/%s:URL=http:://refcache.dnapipelines.sanger.ac.uk::8000/%s

# always activate the environment before proceeding
conda activate isotyper
# main usage
python /path/to/bulkBCRseq/isotyper.py [options]
```

```
usage: isotyper.py [-h] [-i INPUT] [-s STEP] [-l LENGTH] [-dr] [-b] [-c CORES] [-m MEM] [-q QUEUE] [-p PROJECT] [-g GROUP]

options:
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
  -l LENGTH, --length LENGTH
                        minimum length of reads to keep. [Default 100]
  -dr, --dryrun         if passed, prints commands but don't actually run.

bsub arguments:
  -b, --bsub            if passed, submits each row in meta.txt file as a job to bsub.
  -c CORES, --cores CORES
                        number of cores to run this on. [Default 10]
  -m MEM, --mem MEM     job memory request. [Default 8000]
  -q QUEUE, --queue QUEUE
                        job queue to submit to. [Default normal]
  -p PROJECT, --project PROJECT
                        sanger project to send as job. [Default team205]
  -g GROUP, --group GROUP
                        sanger group to send as job. [Default teichlab]
```

## Note!
If you are starting from fastq files directly, please change the 2nd column in the `.txt` file (path to `.cram`) to path to `_R1_001.fastq.gz` (read1) instead. If your read1/read2 suffix isn't this pattern, please modify the `R1PATTERN` and `R2PATTERN` variables file after cloning this repository, in the `_settings.py` directly:
https://github.com/clatworthylab/bulkBCRseq/blob/5d310de8863b64352d68230977c6e7e62d5c0b8f/isotyper/utilities/_settings.py#L25-L27


### Basic usage
```bash
# initial QC
python isotyper.py -i meta.txt -s 1
# trimming
python isotyper.py -i meta.txt -s 2
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
python isotyper.py -i meta.txt -s 2 --bsub
# generate network
python isotyper.py -i meta.txt -s 3 --bsub
# generate network statistic
python isotyper.py -i meta.txt -s 4 --bsub
```

Take a look [here](https://github.com/clatworthylab/bulkBCRseq/tree/master/tests/data) for example files to provide to the tool.


### Post-processing

After running steps `1` to `4`, please annotate the `Fully_reduced_{sample_id}.fasta` file for downstream analysis. You can annotate with [IMGT/HighV-QUEST](https://imgt.org/HighV-QUEST/home.action) or via other software e.g. [MiXCR](https://mixcr.readthedocs.io/en/latest/) in shotgun mode.

```bash
mixcr analyze shotgun -s hsa --starting-material rna --receptor-type igh Fully_reduced_{sample_id}.fasta {sample_id} 
# export to AIRR format
mixcr exportAirr --imgt-gaps in.[vdjca|clns|clna] out.tsv
```

To generate the network plots, you would use the node table (`Att_{sample_id}.txt`) and edge table (`Edges_{sample_id}.txt`) and feed it into a graphing software e.g. `networkx`/`igraph` and continue as per normal. The `orphan` folder has example scripts (probably buggy) on how to use `python-igraph` to generate the plots.

