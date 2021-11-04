# bulk_BCR_analysis
Package belongs to Rachael Bashford-Rogers. See also https://github.com/rbr1/BCR_TCR_PROCESSING_PIPELINE
Requires python 2.7 (or 3.7 if using the python3 branch). Currently only works when cloned onto farm with all paths set up pointing to this folder properly.

## Pre-requisites:
```bash
# create a conda virtual environment
# sample for python 2 set up, swtich to python 3 where appropriate
# install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
eval "$(/path/to/miniconda2/bin/conda shell.bash hook)"
conda init
conda create --name py2 python=2.7 pip numpy pandas networkx
```

Usage instructions on Farm:
```bash
conda activate py2
```
## Basic usage:
```bash
python Processing_sequences_large_scale.py [sample file list] [commands (comma separated list)] [bsub command: Y/N] [print commands: Y/N] [run commands: Y/N]
```
Available commands: 1, 2, 3, 4

### Basic analysis: 1 - Converting raw sequencing files to fastq, QC
```bash
python Processing_sequences_large_scale.py Samples_Mouse_Zach.txt 1 Y Y Y
```

### Basic analysis: 2 - Trimming and filtering reads
```bash
python Processing_sequences_large_scale.py Samples_Mouse_Zach.txt 2 Y Y Y
```
### Basic analysis: 3 - Network generation
```bash
python Processing_sequences_large_scale.py Samples_Mouse_Zach.txt 3 Y Y Y
```
### Basic analysis: 4 - Generating network and population statistics
```bash
python Processing_sequences_large_scale.py Samples_Mouse_Zach.txt 4 Y Y Y
```

## Advanced usage - some private adjustments - not complete!:
```bash
python Processing_sequences_large_scale.py [sample file list] [concat file list] [commands (comma separated list)] [bsub command: Y/N] [print commands: Y/N] [run commands: Y/N]
```
Available commands: 3.5, 3.51
### Create the network from fully reduced fasta sequences: 3.5
```bash
python Processing_sequences_large_scale.py Samples_Mouse_DSS_2020.txt Samples_Mouse_DSS_2020_combined.txt 3.5 Y Y Y
```
### rerun the network generation pipeline using AIRR files: 3.51
```bash
python Processing_sequences_large_scale.py Samples_Mouse_DSS_2020.txt Samples_Mouse_DSS_2020_combined.txt 3.51 Y Y Y
```