import sys
from typing import List, Tuple
import re
import numpy as np
from abnumber.exceptions import ChainParseError
try:
    from anarci.anarci import anarci
except ImportError:
    # Only print the error without failing - required to import
    print('ANARCI module not available. Please install it separately or install AbNumber through Bioconda')
    print('See: https://abnumber.readthedocs.io/')
    sys.exit(1)

POS_REGEX = re.compile(r'([HL]?)(\d+)([A-Z]?)')
WHITESPACE = re.compile(r'\s+')


def _validate_chain_type(chain_type):
    assert chain_type in ['H', 'L', 'K'], \
        f'Invalid chain type "{chain_type}", it should be "H" (heavy),  "L" (lambda light chian) or "K" (kappa light chain)'


def _anarci_align(sequence, scheme, allowed_species, assign_germline=False) -> List[Tuple]:
    from abnumber.position import Position
    sequence = re.sub(WHITESPACE, '', sequence)
    all_numbered, all_ali, all_hits = anarci(
        [('id', sequence)],
        scheme=scheme,
        allowed_species=allowed_species,
        assign_germline=assign_germline
    )
    seq_numbered = all_numbered[0]
    seq_ali = all_ali[0]
    if seq_numbered is None:
        raise ChainParseError(f'Variable chain sequence not recognized: "{sequence}"')
    assert len(seq_numbered) == len(seq_ali), 'Unexpected ANARCI output'
    results = []
    for (positions, start, end), ali in zip(seq_numbered, seq_ali):
        chain_type = ali['chain_type']
        species = ali['species']
        v_gene = ali['germlines']['v_gene'][0][1] if assign_germline else None
        j_gene = ali['germlines']['j_gene'][0][1] if assign_germline else None
        aa_dict = {Position(chain_type=chain_type, number=num, letter=letter, scheme=scheme): aa
                   for (num, letter), aa in positions if aa != '-'}
        tail = sequence[end+1:]
        results.append((aa_dict, chain_type, tail, species, v_gene, j_gene))
    return results


def _get_unique_chains(chains):
    seqs = set()
    chains_filtered = []
    for chain in chains:
        if chain.seq in seqs:
            continue
        seqs.add(chain.seq)
        chains_filtered.append(chain)
    return chains_filtered


# Based on positive score in Blosum62
SIMILAR_PAIRS = {'AA', 'AS', 'CC', 'DD', 'DE', 'DN', 'ED', 'EE', 'EK', 'EQ', 'FF', 'FW', 'FY', 'GG', 'HH', 'HN', 'HY',
                 'II', 'IL', 'IM', 'IV', 'KE', 'KK', 'KQ', 'KR', 'LI', 'LL', 'LM', 'LV', 'MI', 'ML', 'MM', 'MV', 'ND',
                 'NH', 'NN', 'NS', 'PP', 'QE', 'QK', 'QQ', 'QR', 'RK', 'RQ', 'RR', 'SA', 'SN', 'SS', 'ST', 'TS', 'TT',
                 'VI', 'VL', 'VM', 'VV', 'WF', 'WW', 'WY', 'YF', 'YH', 'YW', 'YY'}


def is_similar_residue(a, b):
    if a == '-' or b == '-':
        return a == b
    return a+b in SIMILAR_PAIRS


def is_integer(object):
    return isinstance(object, int) or isinstance(object, np.integer)


SUPPORTED_SCHEMES = ['imgt', 'aho', 'chothia', 'kabat']
SUPPORTED_CDR_DEFINITIONS = ['imgt', 'chothia', 'kabat', 'north']

SCHEME_BORDERS = {
               # Start coordinates
               # CDR1, FR2, CDR2, FR3, CDR3, FR4
         'imgt': [27,  39,  56,   66,  105,  118, 129],
      'kabat_H': [31,  36,  50,   66,  95,   103, 114],
      'kabat_K': [24,  35,  50,   57,  89,    98, 108],
      'kabat_L': [24,  35,  50,   57,  89,    98, 108],
    'chothia_H': [26,  33,  52,   57,  95,   103, 114],
    'chothia_K': [24,  35,  50,   57,  89,    98, 108],
    'chothia_L': [24,  35,  50,   57,  89,    98, 108],
      'north_H': [23,  36,  50,   59,  93,   103, 114],
      'north_K': [24,  35,  49,   57,  89,    98, 108],
      'north_L': [24,  35,  49,   57,  89,    98, 108],
}

# { scheme -> { region -> list of position numbers } }
SCHEME_REGIONS = {
    scheme: {
        'FR1': list(range(1, borders[0])),
        'CDR1': list(range(borders[0], borders[1])),
        'FR2': list(range(borders[1], borders[2])),
        'CDR2': list(range(borders[2], borders[3])),
        'FR3': list(range(borders[3], borders[4])),
        'CDR3': list(range(borders[4], borders[5])),
        'FR4': list(range(borders[5], borders[6])),
    } for scheme, borders in SCHEME_BORDERS.items()
}

# { scheme -> { position number -> region } }
SCHEME_POSITION_TO_REGION = {
    scheme: {pos_num: region for region, positions in regions.items() for pos_num in positions} \
    for scheme, regions in SCHEME_REGIONS.items()
}

# { scheme -> set of vernier position numbers }
SCHEME_VERNIER = {
    #    'imgt_H': frozenset([2,                 52, 53, 54, 76, 78, 80, 82, 87,         118]),
    # 'chothia_H': frozenset([2,                 47, 48, 49, 67, 69, 71, 73, 78, 93, 94, 103]),
    #   'north_H': frozenset([2,                 47, 48, 49, 67, 69, 71, 73, 78, 93, 94, 103]),
      'kabat_H': frozenset([2, 27, 28, 29, 30, 47, 48, 49, 67, 69, 71, 73, 78, 93, 94, 103]),

    #    'imgt_K': frozenset([2, 4, 41, 42, 52, 53, 54, 55, 78, 80, 84, 85, 87, 118]),
    #    'imgt_L': frozenset([2, 4, 41, 42, 52, 53, 54, 55, 78, 80, 84, 85, 87, 118]),
    # 'chothia_K': frozenset([2, 4, 35, 36, 46, 47, 48, 49, 64, 66, 68, 69, 71, 98]),
    # 'chothia_L': frozenset([2, 4, 35, 36, 46, 47, 48, 49, 64, 66, 68, 69, 71, 98]),
    #   'north_K': frozenset([2, 4, 35, 36, 46, 47, 48, 49, 64, 66, 68, 69, 71, 98]),
    #   'north_L': frozenset([2, 4, 35, 36, 46, 47, 48, 49, 64, 66, 68, 69, 71, 98]),
      'kabat_K': frozenset([2, 4, 35, 36, 46, 47, 48, 49, 64, 66, 68, 69, 71, 98]),
      'kabat_L': frozenset([2, 4, 35, 36, 46, 47, 48, 49, 64, 66, 68, 69, 71, 98]),
}

#'kabat_H': 31-35, 50-65, 95-102
#'kabat_K': 24-34, 50-56, 89-97
