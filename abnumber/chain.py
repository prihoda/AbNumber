from collections import OrderedDict
from Bio.SubsMat import MatrixInfo
import sys
try:
    from anarci.anarci import anarci
except ImportError:
    # Only print the error without failing - required to import
    print('ANARCI module not available. Please install it separately or install AbNumber through Bioconda')
    print('See: https://abnumber.readthedocs.io/')
    sys.exit(1)
from termcolor import colored as colored_fn
from Bio.Seq import Seq


class Chain:
    """
    Structure storing an aligned antibody chain using a chosen antibody numbering scheme

    Use `Chain('EVQLQV...')` to create a Chain from an unaligned amino acid sequence
    """

    def __init__(self, sequence, name=None, scheme='imgt', allowed_species=None, aa_dict=None, chain_type=None):
        """
        Create a Chain object from an unaligned string sequence

        :param sequence: Unaligned string sequence
        :param name: Optional sequence identifier
        :param scheme: Numbering scheme: FIXME available numbering schemes
        :param allowed_species: None or one or more of: 'human', 'mouse','rat','rabbit','rhesus','pig','alpaca'
        :param aa_dict: Create Chain object directly from dictionary of region objects (internal use)
        :param chain_type: Explicitly assign chain type, used together with aa_dict= (internal use)
       """

        assert scheme in ['imgt']

        if aa_dict is not None:
            if sequence is not None:
                raise ValueError('Only one of aa_dict= and sequence= can be provided')
            assert isinstance(aa_dict, dict), f'Expected dict, got: {type(aa_dict)}'
        else:
            if chain_type is not None:
                raise ValueError('Do not use chain_type= when providing sequence=, it will be inferred automatically')
            if isinstance(sequence, Seq):
                sequence = str(sequence)
            aa_dict, chain_type = _anarci_align(sequence, scheme=scheme, allowed_species=allowed_species)

        _validate_chain_type(chain_type)

        self.name = name
        self.chain_type = chain_type
        self.scheme = scheme

        self.fw1_dict = OrderedDict()
        self.cdr1_dict = OrderedDict()
        self.fw2_dict = OrderedDict()
        self.cdr2_dict = OrderedDict()
        self.fw3_dict = OrderedDict()
        self.cdr3_dict = OrderedDict()
        self.fw4_dict = OrderedDict()

        self._init_from_dict(aa_dict)

    def _init_from_dict(self, aa_dict):
        regions_list = [self.fw1_dict, self.cdr1_dict, self.fw2_dict, self.cdr2_dict, self.fw3_dict, self.cdr3_dict, self.fw4_dict]
        region = 0
        for pos in sorted(aa_dict.keys()):
            while pos.number >= IMGT_BORDERS[region]:
                region += 1
            aa = aa_dict[pos].upper()
            # TODO validate amino acid aa
            # TODO use Region class instead of str?
            regions_list[region][pos] = aa

    def __repr__(self):
        return self.seq

    def format(self, method='wide'):
        if method == 'wide':
            return self.to_wide_string()
        elif method == 'tall':
            return self.to_tall_string()
        raise ValueError(f'Use method="wide" or method="tall", unknown method: "{method}"')

    def to_tall_string(self):
        seq = []
        for region, aa_dict in self.regions.items():
            for pos, aa in aa_dict.items():
                seq.append(f'{region: >4} {pos: <4} {aa}')
        return '\n'.join(seq)

    def to_wide_string(self):
        annot = ' ' * len(self.fw1_dict)
        annot += '^' * len(self.cdr1_dict)
        annot += ' ' * len(self.fw2_dict)
        annot += '^' * len(self.cdr2_dict)
        annot += ' ' * len(self.fw3_dict)
        annot += '^' * len(self.cdr3_dict)
        annot += ' ' * len(self.fw4_dict)
        return self.seq + '\n' + annot

    def is_heavy_chain(self):
        return self.chain_type == 'H'

    def is_light_chain(self):
        return self.is_lambda_light_chain() or self.is_kappa_light_chain()

    def is_lambda_light_chain(self):
        return self.chain_type == 'L'

    def is_kappa_light_chain(self):
        return self.chain_type == 'K'

    def has_same_cdr_positions(self, other):
        if len(self.cdr1_dict) != len(other.cdr1) or self.cdr1_dict.keys() != other.cdr1.keys():
            return False
        if len(self.cdr2_dict) != len(other.cdr2) or self.cdr2_dict.keys() != other.cdr2.keys():
            return False
        if len(self.cdr3_dict) != len(other.cdr3) or self.cdr3_dict.keys() != other.cdr3.keys():
            return False
        return True

    def get_fw1_matches(self, other):
        return sum(aa == other.fw1.get(pos) for pos, aa in self.fw1_dict.items())

    def get_cdr1_matches(self, other):
        return sum(aa == other.cdr1.get(pos) for pos, aa in self.cdr1_dict.items())

    def get_fw2_matches(self, other):
        return sum(aa == other.fw2.get(pos) for pos, aa in self.fw2_dict.items())

    def get_cdr2_matches(self, other):
        return sum(aa == other.cdr2.get(pos) for pos, aa in self.cdr2_dict.items())

    def get_fw3_matches(self, other):
        return sum(aa == other.fw3.get(pos) for pos, aa in self.fw3_dict.items())

    def get_cdr3_matches(self, other):
        return sum(aa == other.cdr3.get(pos) for pos, aa in self.cdr3_dict.items())

    def get_fw4_matches(self, other):
        return sum(aa == other.fw4.get(pos) for pos, aa in self.fw4_dict.items())

    def get_fw_matches(self, other):
        return self.get_fw1_matches(other) + self.get_fw2_matches(other) + self.get_fw3_matches(
            other) + self.get_fw4_matches(other)

    def get_cdr_matches(self, other):
        return self.get_cdr1_matches(other) + self.get_cdr2_matches(other) + self.get_cdr3_matches(other)

    def get_matches(self, other):
        return self.get_fw_matches(other) + self.get_cdr_matches(other)

    def align(self, *others) -> 'Alignment':
        """
        Align this chain to other chains by using their existing numbering
        :param others: Chains to align
        :return: Alignment object
        """
        pos_dicts = [self.positions]
        for other in others:
            assert isinstance(other, Chain), f'Expected Chain object, got {type(other)}: {other}'
            pos_dicts.append(other.positions)
        shared_pos = sorted(set(pos for pos_dict in pos_dicts for pos in pos_dict.keys()))
        residues = [tuple(pos_dict.get(pos, '-') for pos_dict in pos_dicts) for pos in shared_pos]
        return Alignment(shared_pos, residues)

    def clone(self, replace_seq: str = None):
        """
        Create a copy of this chain, optionally with a replacement sequence that is placed into the same numbering
        :param replace_seq: Optional replacement sequence, needs to be the same length
        :return: new Chain object
        """
        aa_dict = {}
        positions = self.positions
        if replace_seq is not None:
            assert len(replace_seq) == len(positions), 'Sequence needs to be the same length'
        for i, (pos, aa) in enumerate(positions.items()):
            aa_dict[pos] = replace_seq[i] if replace_seq is not None else aa
        return Chain(sequence=None, aa_dict=aa_dict, name=self.name, scheme=self.scheme, chain_type=self.chain_type)

    def graft_cdrs_onto(self, other: 'Chain', name: str = None) -> 'Chain':
        """
        Graft CDRs from this Chain onto another chain
        :param other: Chain to graft CDRs into
        :param name: Name of new Chain
        :return: Chain with CDRs grafted from this chain and frameworks from the given chain
        """
        assert self.scheme == other.scheme, \
            f'Sequences need to have the same numbering scheme, got {self.scheme} and {other.scheme}'
        assert self.chain_type == other.chain_type, \
            f'Sequences need to have the same chain type, got {self.chain_type} and {other.chain_type}'

        grafted_dict = {}
        for (self_region, self_dict), (other_region, other_dict) in zip(self.regions.items(), other.regions.items()):
            assert self_region == other_region
            # TODO use Region class instead of str?
            if self_region.upper().startswith('CDR'):
                grafted_dict.update(self_dict)
            else:
                grafted_dict.update(other_dict)

        return Chain(sequence=None, aa_dict=grafted_dict, name=name, chain_type=self.chain_type, scheme=self.scheme)

    @property
    def regions(self):
        """
        Get dictionary of region dictionaries (Position -> Amino acid)
        :return: Dictionary of Region name -> dictionary of (Position -> Amino acid)
        """
        return OrderedDict(
            fw1=self.fw1_dict,
            cdr1=self.cdr1_dict,
            fw2=self.fw2_dict,
            cdr2=self.cdr2_dict,
            fw3=self.fw3_dict,
            cdr3=self.cdr3_dict,
            fw4=self.fw4_dict
        )

    @property
    def positions(self):
        """
        Get dictionary of Position -> Amino acid
        :return: dictionary of Position -> Amino acid
        """
        positions = OrderedDict()
        for region, aa_dict in self.regions.items():
            for pos, aa in aa_dict.items():
                positions[pos] = aa
        return positions

    @property
    def seq(self):
        """
        Get unaligned string representation of the variable chain sequence
        :return: Unaligned string representation of the variable chain sequence
        """
        return ''.join(self.positions.values())

    @property
    def fw1_seq(self):
        return ''.join(self.fw1_dict.values())

    @property
    def cdr1_seq(self):
        return ''.join(self.cdr1_dict.values())

    @property
    def fw2_seq(self):
        return ''.join(self.fw2_dict.values())

    @property
    def cdr2_seq(self):
        return ''.join(self.cdr2_dict.values())

    @property
    def fw3_seq(self):
        return ''.join(self.fw3_dict.values())

    @property
    def cdr3_seq(self):
        return ''.join(self.cdr3_dict.values())

    @property
    def fw4_seq(self):
        return ''.join(self.fw4_dict.values())


class Alignment:
    def __init__(self, positions, residues):
        assert len(positions) == len(residues)
        self.positions = positions
        self.residues = residues
        self._zipped = list(zip(self.positions, self.residues))

    def __repr__(self):
        return self.format(colored=False)

    def __iter__(self):
        yield from self._zipped.__iter__()

    def __getitem__(self, item):
        args = (self.positions[item], self.residues[item])
        if isinstance(item, slice):
            return Alignment(*args)
        return args

    def __len__(self):
        return len(self.positions)

    def format(self, mark_identity=True, mark_cdrs=True, colored=False):
        seq1 = ''
        identity = ''
        seq2 = ''
        cdrs = ''
        # TODO support multiple sequence alignment
        for pos, (a, b) in self:
            if not colored or a == b:
                seq1 += a
                seq2 += b
            elif is_similar_residue(a, b):
                seq1 += colored_fn(a, 'white', 'on_yellow', attrs=['bold'])
                seq2 += colored_fn(b, 'white', 'on_yellow', attrs=['bold'])
            else:
                seq1 += colored_fn(a, 'white', 'on_red', attrs=['bold'])
                seq2 += colored_fn(b, 'white', 'on_red', attrs=['bold'])

            if mark_identity:
                identity += '|' if a == b else ('+' if is_similar_residue(a, b) else '.')
            if mark_cdrs:
                cdrs += '^' if pos.is_in_cdr() else ' '
        return seq1 + (('\n' + identity) if mark_identity else '') + '\n' + seq2 + (
            ('\n' + cdrs) if mark_cdrs else '')

    def has_mutation(self):
        return any(len(set(aas)) != 1 for aas in self.residues)


class Position:
    def __init__(self, chain_type: str, number: int, letter: str, scheme: str):
        _validate_chain_type(chain_type)
        self.chain_type: str = chain_type
        self.number: int = int(number)
        self.letter: str = letter
        self.scheme: str = scheme

    def __repr__(self):
        return f'{self.chain_type_format()}{self.number}{self.letter}({self.scheme})'

    def __str__(self):
        return f'{self.chain_type_format()}{self.number}{self.letter}'

    def __hash__(self):
        return self.__repr__().__hash__()

    def __eq__(self, other):
        return self.__repr__() == other.__repr__()

    def __lt__(self, other):
        assert self.chain_type == other.chain_type, f'Positions do not come from the same chain: {self}, {other}'
        assert self.scheme == other.scheme, 'Comparing positions in different schemes is not implemented'
        return self._sort_key() < other._sort_key()

    def chain_type_format(self):
        if self.chain_type == 'H':
            return 'H'
        if self.chain_type in ['K', 'L']:
            return 'L'
        raise NotImplementedError(f'Unknown chain type "{self.chain_type}"')

    def _sort_key(self):
        if self.scheme == 'imgt':
            letter_ord = ord(self.letter) if self.letter else 0
            if self.number == 112:
                # position 112 is sorted in reverse
                letter_ord = -letter_ord
        else:
            raise NotImplementedError(f'Cannot compare positions of scheme: {self.scheme}')
        return self.number, letter_ord

    def get_region(self):
        if self.scheme == 'imgt':
            return IMGT_POS_DICT[self.number]
        else:
            raise NotImplementedError(f'Not supported scheme: {self.scheme}')

    def is_in_cdr(self):
        # FIXME
        return self.get_region().lower().startswith('cdr')


def _validate_chain_type(chain_type):
    assert chain_type in ['H', 'L', 'K'], \
        f'Invalid chain type "{chain_type}", it should be "H" (heavy),  "L" (lambda light chian) or "K" (kappa light chain)'


def _anarci_align(sequence, scheme, allowed_species):
    all_numbered, all_ali, all_hits = anarci([('id', sequence)], scheme=scheme, allowed_species=allowed_species)
    # We only have one sequence
    numbered = all_numbered[0]
    ali = all_ali[0]
    hits = all_hits[0]
    if numbered is None:
        raise ValueError(f'No alignment found for sequence: "{sequence}"')
    if len(numbered) != 1:
        raise NotImplementedError(f'Unsupported: Multiple ANARCI domains found in sequence: "{sequence}"')
    positions, start, end = numbered[0]
    chain_type = ali[0]['chain_type']
    aa_dict = {Position(chain_type=chain_type, number=num, letter=letter, scheme=scheme): aa for (num, letter), aa in
               positions if aa != '-'}
    return aa_dict, chain_type


def is_similar_residue(a, b, matrix=MatrixInfo.blosum62):
    if a == '-' or b == '-':
        return a == b
    pair = (a, b) if (a, b) in matrix else (b, a)
    return matrix[pair] > 0


IMGT_BORDERS = [27, 39, 56, 66, 105, 118, 129]

IMGT_FW1 = list(range(1, IMGT_BORDERS[0]))
IMGT_CDR1 = list(range(IMGT_BORDERS[0], IMGT_BORDERS[1]))
IMGT_FW2 = list(range(IMGT_BORDERS[1], IMGT_BORDERS[2]))
IMGT_CDR2 = list(range(IMGT_BORDERS[2], IMGT_BORDERS[3]))
IMGT_FW3 = list(range(IMGT_BORDERS[3], IMGT_BORDERS[4]))
IMGT_CDR3 = list(range(IMGT_BORDERS[4], IMGT_BORDERS[5]))
IMGT_FW4 = list(range(IMGT_BORDERS[5], IMGT_BORDERS[6]))

IMGT_CDR = IMGT_CDR1 + IMGT_CDR2 + IMGT_CDR3
IMGT_FW = IMGT_FW1 + IMGT_FW2 + IMGT_FW3 + IMGT_FW4

IMGT_CDR_DICT = {'cdr1': IMGT_CDR1, 'cdr2': IMGT_CDR2, 'cdr3': IMGT_CDR3}
IMGT_FW_DICT = {'fw1': IMGT_FW1, 'fw2': IMGT_FW2, 'fw3': IMGT_FW3, 'fw4': IMGT_FW4}
IMGT_REGION_DICT = {**IMGT_CDR_DICT, **IMGT_FW_DICT}
IMGT_POS_DICT = {pos_num: region for region, positions in IMGT_REGION_DICT.items() for pos_num in positions}
