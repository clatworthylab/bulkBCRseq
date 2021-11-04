# bulk_BCR_analysis
Bulk BCR-seq processing scripts use in Fitzpatrick et al., Nature (2020). Package belongs to Rachael Bashford-Rogers.

This repo is an older version of what seems to be now at https://github.com/rbr1/BCR_TCR_PROCESSING_PIPELINE.

Requires python 2.7 (or 3.7 if using the python3 branch). Currently only works when cloned onto farm with all paths set up pointing to this folder properly.

## Citation
Please cite the following papers:

```
Fitzpatrick, Z., Frazer, G., Ferro, A., Clare, S., Bouladoux, N., Ferdinand, J., Tuong, Z.K., Negro-Demontel, M.L., Kumar, N., Suchanek, O. and Tajsic, T., 2020. Gut-educated IgA plasma cells defend the meningeal venous sinuses. Nature, 587(7834), pp.472-476.
```

```
Bashford-Rogers, R.J., Palser, A.L., Huntly, B.J., Rance, R., Vassiliou, G.S., Follows, G.A. and Kellam, P., 2013. Network properties derived from deep sequencing of human B-cell receptor repertoires delineate B-cell populations. Genome research, 23(11), pp.1874-1884.
```

```
Bashford-Rogers, R.J.M., Bergamaschi, L., McKinney, E.F., Pombal, D.C., Mescia, F., Lee, J.C., Thomas, D.C., Flint, S.M., Kellam, P., Jayne, D.R.W. and Lyons, P.A., 2019. Analysis of the B cell receptor repertoire in six immune-mediated diseases. Nature, 574(7776), pp.122-126.
```


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