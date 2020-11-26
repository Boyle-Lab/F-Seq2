## Installation

Binary installers for the latest released version are available at the [python package index](https://pypi.org/project/pandas) and on conda.  
The source code is currently hosted on GitHub at: https://github.com/Boyle-Lab/F-Seq2

### Prerequisites
- python >= 3.6
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) >= 2.29.0 (required by pybedtools)


### Prepare a conda environment (optional)
```
$ conda create --name fseq2_env python=3.6
$ conda activate fseq2_env
```

### Install from binary installer

```
# PyPI
$ pip install fseq2

# or conda
$ conda install -c bioconda -c dmentipl -c samzhao fseq2
```

### Install from source

```
$ wget https://github.com/Boyle-Lab/F-Seq2/archive/master.zip

$ unzip master.zip

$ cd F-Seq2-master

# install all dependencies listed below

$ python setup.py install
```

### Integration test
```
$ cd Fseq2

$ python setup.py test
```
This may take 1 min


### Packages installed after installation

- numpy >= 1.15.4
- scipy >= 1.1.0
- pandas >= 0.24
- statsmodels
- KDEpy
- pybedtools
- pyBigWig
- h5py
- matplotlib
