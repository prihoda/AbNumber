# AbNumber

Convenience Python APIs for antibody numbering using [ANARCI](https://github.com/oxpig/ANARCI)

## Installation

Install using Bioconda:
```bash
conda install -c bioconda abnumber
```

## Usage

```python
from abnumber import Chain

chain = Chain('ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIKRTV')

print(chain.format())
# ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIK
#                           ^^^^^^                 ^^^                                    ^^^^^^^^^          

print(chain.cdr3_seq)
# QQFNSYPLT 
```

## Credits

See [ANARCI on GitHub](https://github.com/oxpig/ANARCI) and the ANARCI paper: [ANARCI: antigen receptor numbering and receptor classification](https://doi.org/10.1093/bioinformatics/btv552)