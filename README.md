# AbNumber

[![Build & Test](https://github.com/prihoda/AbNumber/actions/workflows/build-test.yml/badge.svg)](https://github.com/prihoda/AbNumber/actions/workflows/build-test.yml)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/abnumber.svg?style=flag&label=BioConda%20install&color=green)](https://anaconda.org/bioconda/abnumber) 
[![Docs](https://readthedocs.org/projects/pip/badge/?version=latest&style=flat)](http://abnumber.readthedocs.org/)

Convenience Python APIs for antibody numbering and alignment using [ANARCI](https://github.com/oxpig/ANARCI)

Try it out in your browser using Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/prihoda/AbNumber/HEAD?filepath=examples%2FAbNumber_getting_started.ipynb)

## Features:

- **Streamlined Python API** using `Chain` object
- **Identifying CDR regions**, e.g. using `chain.regions` or `chain.cdr3_seq`
- **Indexing and slicing**, e.g. using `chain['5']` or `chain['H2':'H5']`
- **Pairwise and multiple sequence alignment** in the given numbering using `chain.align(another)`
- **Alignment to nearest human germline** using `chain.align(chain.find_merged_human_germline())`
- **Humanization using CDR grafting** by `chain.graft_cdrs_onto_human_germline()`

See [AbNumber Documentation](https://abnumber.readthedocs.io/en/latest/) for the full reference.

## Installation

Install using Bioconda:
```bash
conda install -c bioconda abnumber
```

Note: **Windows is not supported** due to [HMMER](https://github.com/EddyRivasLab/hmmer) dependency. AbNumber is currently only available on UNIX & MacOS. 

## Examples

```python
from abnumber import Chain

seq = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSSAKTTAPSVYPLA'
chain = Chain(seq, scheme='imgt')

chain
# QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
#                          ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^

chain.cdr3_seq
# ARYYDDHYCLDY

chain.print(numbering=True)
# 0        1        2         3     4         5         6       7        8         9         10        11       12       
# 12345678912345678901234567890567890123456789012345678923456789012456789012345678901234567890123456789023456789012345678
# QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
#                          ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^           
```

Chain can be iterated:

```python
for pos, aa in chain:
    print(pos, aa)
# H1  Q
# H2  V
# H3  Q
# H4  L
# H5  Q
```

Chain can also be indexed and sliced using scheme numbering:

```python
chain['5']
# 'Q'
for pos, aa in chain['H2':'H5']:
    print(pos, aa)
# H2  V
# H3  Q
# H4  L
# H5  Q
```

For all methods see [AbNumber Documentation](https://abnumber.readthedocs.io/en/latest/)

## Credits

See [ANARCI on GitHub](https://github.com/oxpig/ANARCI) and the ANARCI paper: [ANARCI: antigen receptor numbering and receptor classification](https://doi.org/10.1093/bioinformatics/btv552)
