# AbNumber

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/abnumber.svg?style=flag&label=BioConda%20install&color=green)](https://anaconda.org/bioconda/abnumber) 
[![Docs](https://readthedocs.org/projects/pip/badge/?version=latest&style=flat)](http://abnumber.readthedocs.org/)

Convenience Python APIs for antibody numbering using [ANARCI](https://github.com/oxpig/ANARCI)

## Installation

Install using Bioconda:
```bash
conda install -c bioconda abnumber
```

## Documentation

See [AbNumber Documentation](https://abnumber.readthedocs.io/en/latest/)

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
