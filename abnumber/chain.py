from collections import OrderedDict
from Bio.SubsMat import MatrixInfo
from anarci.anarci import anarci
from termcolor import colored as colored_fn


class Chain:

    def __init__(self, aa_dict, name=None):
        self.name = name
        self.fw1 = OrderedDict()
        self.cdr1 = OrderedDict()
        self.fw2 = OrderedDict()
        self.cdr2 = OrderedDict()
        self.fw3 = OrderedDict()
        self.cdr3 = OrderedDict()
        self.fw4 = OrderedDict()

        regions_list = [self.fw1, self.cdr1, self.fw2, self.cdr2, self.fw3, self.cdr3, self.fw4]
        region = 0
        for pos in sorted(aa_dict.keys()):
            while pos.number >= IMGT_BORDERS[region]:
                region += 1
            aa = aa_dict[pos].upper()
            # TODO validate amino acid aa
            regions_list[region][pos] = aa

    def clone(self, replace_seq=None):
        aa_dict = {}
        positions = self.positions
        if replace_seq is not None:
            assert len(replace_seq) == len(positions), 'Sequence needs to be the same length'
        for i, (pos, aa) in enumerate(positions.items()):
            aa_dict[pos] = replace_seq[i] if replace_seq is not None else aa
        return Chain(aa_dict=aa_dict, name=self.name)

    def __repr__(self):
        return self.seq

    def get(self, position, default=None):
        for region, aa_dict in self.regions.items():
            if position in aa_dict:
                return aa_dict[position]
        return default

    def __getitem__(self, position):
        aa = self.get(position)
        if aa is None:
            raise IndexError(f'Position "{position}" not found in chain')
        return aa

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
        annot = ' ' * len(self.fw1)
        annot += '^' * len(self.cdr1)
        annot += ' ' * len(self.fw2)
        annot += '^' * len(self.cdr2)
        annot += ' ' * len(self.fw3)
        annot += '^' * len(self.cdr3)
        annot += ' ' * len(self.fw4)
        return self.seq + '\n' + annot

    def has_same_cdr_positions(self, other):
        if len(self.cdr1) != len(other.cdr1) or self.cdr1.keys() != other.cdr1.keys():
            return False
        if len(self.cdr2) != len(other.cdr2) or self.cdr2.keys() != other.cdr2.keys():
            return False
        if len(self.cdr3) != len(other.cdr3) or self.cdr3.keys() != other.cdr3.keys():
            return False
        return True

    def get_fw1_matches(self, other):
        return sum(aa == other.fw1.get(pos) for pos, aa in self.fw1.items())

    def get_cdr1_matches(self, other):
        return sum(aa == other.cdr1.get(pos) for pos, aa in self.cdr1.items())

    def get_fw2_matches(self, other):
        return sum(aa == other.fw2.get(pos) for pos, aa in self.fw2.items())

    def get_cdr2_matches(self, other):
        return sum(aa == other.cdr2.get(pos) for pos, aa in self.cdr2.items())

    def get_fw3_matches(self, other):
        return sum(aa == other.fw3.get(pos) for pos, aa in self.fw3.items())

    def get_cdr3_matches(self, other):
        return sum(aa == other.cdr3.get(pos) for pos, aa in self.cdr3.items())

    def get_fw4_matches(self, other):
        return sum(aa == other.fw4.get(pos) for pos, aa in self.fw4.items())

    def get_fw_matches(self, other):
        return self.get_fw1_matches(other) + self.get_fw2_matches(other) + self.get_fw3_matches(
            other) + self.get_fw4_matches(other)

    def get_cdr_matches(self, other):
        return self.get_cdr1_matches(other) + self.get_cdr2_matches(other) + self.get_cdr3_matches(other)

    def get_matches(self, other):
        return self.get_fw_matches(other) + self.get_cdr_matches(other)

    def align(self, *others):
        pos_dicts = [self.positions]
        for other in others:
            pos_dicts.append(other.positions)
        shared_pos = sorted(set(pos for pos_dict in pos_dicts for pos in pos_dict.keys()))
        residues = [tuple(pos_dict.get(pos, '-') for pos_dict in pos_dicts) for pos in shared_pos]
        return Alignment(shared_pos, residues)

    @classmethod
    def from_str(cls, seq_str, scheme='imgt', allowed_species=None, name=None):
        # Allowed species: ['human', 'mouse','rat','rabbit','rhesus','pig','alpaca']
        assert scheme in ['imgt']
        all_numbered, all_ali, all_hits = anarci([('id', seq_str)], scheme=scheme, allowed_species=allowed_species)
        # We only have one sequence
        numbered = all_numbered[0]
        ali = all_ali[0]
        hits = all_hits[0]
        if numbered is None:
            raise ValueError(f'No alignment found for sequence: "{seq_str}"')
        if len(numbered) != 1:
            raise NotImplementedError(f'Unsupported: Multiple ANARCI domains found in sequence: "{seq_str}"')
        positions, start, end = numbered[0]
        # FIXME assign chain to cls.chain
        return cls({Position(chain='?', number=num, letter=letter, scheme=scheme): aa for (num, letter), aa in positions if aa != '-'}, name=name)

    @property
    def regions(self):
        return OrderedDict(
            fw1=self.fw1,
            cdr1=self.cdr1,
            fw2=self.fw2,
            cdr2=self.cdr2,
            fw3=self.fw3,
            cdr3=self.cdr3,
            fw4=self.fw4
        )

    @property
    def positions(self):
        positions = OrderedDict()
        for region, aa_dict in self.regions.items():
            for pos, aa in aa_dict.items():
                positions[pos] = aa
        return positions

    @property
    def seq(self):
        return ''.join(self.positions.values())

    @property
    def fw1_seq(self):
        return ''.join(self.fw1.values())

    @property
    def cdr1_seq(self):
        return ''.join(self.cdr1.values())

    @property
    def fw2_seq(self):
        return ''.join(self.fw2.values())

    @property
    def cdr2_seq(self):
        return ''.join(self.cdr2.values())

    @property
    def fw3_seq(self):
        return ''.join(self.fw3.values())

    @property
    def cdr3_seq(self):
        return ''.join(self.cdr3.values())

    @property
    def fw4_seq(self):
        return ''.join(self.fw4.values())


class Alignment:
    def __init__(self, positions, residues):
        self.positions = positions
        self.residues = residues
        self._zipped = list(zip(self.positions, self.residues))

    def __repr__(self):
        return self.format(colored=False)

    def __iter__(self):
        yield from self._zipped.__iter__()

    def __getitem__(self, item):
        return self._zipped.__getitem__(item)

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


class Position:
    def __init__(self, chain, number, letter, scheme):
        self.chain = chain # TODO type
        self.number: str = number
        self.letter: str = letter
        self.scheme = scheme # TODO type

    def __repr__(self):
        return f'{self.chain}{self.number}{self.letter}({self.scheme})'

    def __str__(self):
        return f'{self.chain}{self.number}{self.letter}'

    def __hash__(self):
        return self.__repr__().__hash__()

    def __eq__(self, other):
        return self.__repr__() == other.__repr__()

    def __lt__(self, other):
        assert self.chain == other.chain, f'Positions do not come from the same chain: {self}, {other}'
        assert self.scheme == other.scheme, 'Comparing positions in different schemes is not implemented'
        return self.sort_key < other.sort_key

    @property
    def sort_key(self):
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
