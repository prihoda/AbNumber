# AbNumber

Convenience Python APIs for antibody numbering using [ANARCI](https://github.com/oxpig/ANARCI)

Try it out in your browser using Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/prihoda/AbNumber/HEAD?filepath=examples%2FAbNumber_getting_started.ipynb)

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
