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
    Antibody chain aligned to a chosen antibody numbering scheme

    :example:

    >>> from abnumber import Chain
    >>>
    >>> seq = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSSAKTTAPSVYPLA'
    >>> chain = Chain(seq, scheme='imgt')
    >>>
    >>> print(chain.format())
    QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
                             ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^

    :param sequence: Unaligned string sequence
    :param name: Optional sequence identifier
    :param scheme: Numbering scheme: FIXME available numbering schemes
    :param allowed_species: None or one or more of: 'human', 'mouse','rat','rabbit','rhesus','pig','alpaca'
    :param aa_dict: Create Chain object directly from dictionary of region objects (internal use)
    :param chain_type: Explicitly assign chain type, used together with aa_dict= (internal use)
    """

    def __init__(self, sequence, scheme, name=None, allowed_species=None, aa_dict=None, chain_type=None, tail=None):

        assert scheme in ['imgt']

        if aa_dict is not None:
            if sequence is not None:
                raise ValueError('Only one of aa_dict= and sequence= can be provided')
            assert isinstance(aa_dict, dict), f'Expected dict, got: {type(aa_dict)}'
        else:
            if chain_type is not None:
                raise ValueError('Do not use chain_type= when providing sequence=, it will be inferred automatically')
            if tail is not None:
                raise ValueError('Do not use tail= when providing sequence=, it will be inferred automatically')
            if isinstance(sequence, Seq):
                sequence = str(sequence)
            aa_dict, chain_type, tail = _anarci_align(sequence, scheme=scheme, allowed_species=allowed_species)

        _validate_chain_type(chain_type)

        self.name = name
        """User-provided sequence identifier
        """
        self.chain_type = chain_type
        """Chain type as identified by ANARCI: ``H`` (heavy), ``K`` (kappa light) or ``L`` (lambda light)
        
        See also :meth:`Chain.is_heavy_chain` and :meth:`Chain.is_light_chain`.
        """
        self.scheme = scheme
        """Numbering scheme used to align the sequence
        """
        self.tail = tail
        """Constant region sequence
        """

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
        """Format sequence to string

        By default, produces "wide" format with sequence on first line and CDR regions higlighted with ``^`` on second line

        >>> print(chain.format())
        QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
                                 ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^

        :param method: use ``"wide"`` for :meth:`Chain.wide_format` or ``"tall"`` for :meth:`Chain.tall_format()`
        :return: formatted string
        """
        if method == 'wide':
            return self.wide_format()
        elif method == 'tall':
            return self.tall_format()
        raise ValueError(f'Use method="wide" or method="tall", unknown method: "{method}"')

    def tall_format(self):
        """Create string with one position per line, showing position numbers and amino acids

        >>> print(chain.tall_format())
        fw1 H1    Q
        fw1 H2    V
        fw1 H3    Q
        fw1 H4    L
        fw1 H5    Q
        fw1 H6    Q
        fw1 H7    S
        ...

        :return: formatted string
        """
        seq = []
        for region, aa_dict in self.regions.items():
            for pos, aa in aa_dict.items():
                seq.append(f'{region: >4} {str(pos): <5} {aa}')
        return '\n'.join(seq)

    def wide_format(self):
        """Create string with sequence on first line and CDR regions higlighted with `^` on second line

        >>> print(chain.wide_format())
        QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
                                 ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^

        :return: formatted string
        """
        annot = ' ' * len(self.fw1_dict)
        annot += '^' * len(self.cdr1_dict)
        annot += ' ' * len(self.fw2_dict)
        annot += '^' * len(self.cdr2_dict)
        annot += ' ' * len(self.fw3_dict)
        annot += '^' * len(self.cdr3_dict)
        annot += ' ' * len(self.fw4_dict)
        return self.seq + '\n' + annot

    def is_heavy_chain(self):
        """Check if this chain is heavy chain (``chain_type=="H"``)"""
        return self.chain_type == 'H'

    def is_light_chain(self):
        """Check if this chain is light chain (``chain_type=="K" or chain_type=="L"``)"""
        return self.is_lambda_light_chain() or self.is_kappa_light_chain()

    def is_lambda_light_chain(self):
        """Check if this chain is lambda light chain (``chain_type=="L"``)"""
        return self.chain_type == 'L'

    def is_kappa_light_chain(self):
        """Check if this chain is kappa light chain (``chain_type=="K"``)"""
        return self.chain_type == 'K'

    def has_same_cdr_positions(self, other):
        """Check if this chain has the same set of CDR position numbers

        Used to filter pairs of chains that could be functionally related
        """
        if len(self.cdr1_dict) != len(other.cdr1) or self.cdr1_dict.keys() != other.cdr1.keys():
            return False
        if len(self.cdr2_dict) != len(other.cdr2) or self.cdr2_dict.keys() != other.cdr2.keys():
            return False
        if len(self.cdr3_dict) != len(other.cdr3) or self.cdr3_dict.keys() != other.cdr3.keys():
            return False
        return True

    def get_fw1_matches(self, other):
        """Get number of identical residues at corresponding positions in Framework 1 region"""
        return sum(aa == other.fw1.get(pos) for pos, aa in self.fw1_dict.items())

    def get_cdr1_matches(self, other):
        """Get number of identical residues at corresponding positions in the CDR 1 region"""
        return sum(aa == other.cdr1.get(pos) for pos, aa in self.cdr1_dict.items())

    def get_fw2_matches(self, other):
        """Get number of identical residues at corresponding positions in the Framework 2 region"""
        return sum(aa == other.fw2.get(pos) for pos, aa in self.fw2_dict.items())

    def get_cdr2_matches(self, other):
        """Get number of identical residues at corresponding positions in the Framework 4 region"""
        return sum(aa == other.cdr2.get(pos) for pos, aa in self.cdr2_dict.items())

    def get_fw3_matches(self, other):
        """Get number of identical residues at corresponding positions in the Framework 3 region"""
        return sum(aa == other.fw3.get(pos) for pos, aa in self.fw3_dict.items())

    def get_cdr3_matches(self, other):
        """Get number of identical residues at corresponding positions in the Framework 4 region"""
        return sum(aa == other.cdr3.get(pos) for pos, aa in self.cdr3_dict.items())

    def get_fw4_matches(self, other):
        """Get number of identical residues at corresponding positions in the Framework 4 region"""
        return sum(aa == other.fw4.get(pos) for pos, aa in self.fw4_dict.items())

    def get_fw_matches(self, other):
        """Get number of identical residues at corresponding positions in the Framework regions"""
        return self.get_fw1_matches(other) + self.get_fw2_matches(other) + self.get_fw3_matches(
            other) + self.get_fw4_matches(other)

    def get_cdr_matches(self, other):
        """Get number of identical residues at corresponding positions in all the CDR regions"""
        return self.get_cdr1_matches(other) + self.get_cdr2_matches(other) + self.get_cdr3_matches(other)

    def get_matches(self, other):
        """Get number of identical residues at corresponding positions in the full variable region"""
        return self.get_fw_matches(other) + self.get_cdr_matches(other)

    def align(self, *other) -> 'Alignment':
        """Align this chain to other chains by using their existing numbering

        >>> from abnumber import Chain
        >>>
        >>> seq1 = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSSAKTTAP'
        >>> chain1 = Chain(seq1, scheme='imgt')
        >>>
        >>> seq2 = 'QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYDDYLDRWGQGTTLTVSSAKTTAP'
        >>> chain2 = Chain(seq2, scheme='imgt')
        >>>
        >>> alignment = chain1.align(chain2)
        >>> print(alignment.format())
        QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPS-RGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
        ||||.||||||.||||+|||||||||||.||||||||||||||||+||||||||.|.||||||||||||||||||||||||||.+|||||||||||||||||....||.|||||||||||
        QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYD--DYLDRWGQGTTLTVSS
                                 ^^^^^^^^                 ^^^^^^^^^                                      ^^^^^^^^^^^^

        :param other: The :class:`Chain` object to align, can be repeated to create a multiple sequence alignment
        :return: :class:`Alignment` object
        """
        pos_dicts = [self.positions]
        for chain in other:
            assert isinstance(chain, Chain), f'Expected Chain object, got {type(chain)}: {chain}'
            pos_dicts.append(chain.positions)
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
        """Dictionary of region dictionaries (Position -> Amino acid)
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
        """Dictionary of Position -> Amino acid
        :return: dictionary of Position -> Amino acid
        """
        positions = OrderedDict()
        for region, aa_dict in self.regions.items():
            for pos, aa in aa_dict.items():
                positions[pos] = aa
        return positions

    @property
    def seq(self):
        """Unaligned string representation of the variable chain sequence

        :return: Unaligned string representation of the variable chain sequence
        """
        return ''.join(self.positions.values())

    @property
    def fw1_seq(self):
        """Unaligned string representation of the Framework 1 region sequence"""
        return ''.join(self.fw1_dict.values())

    @property
    def cdr1_seq(self):
        """Unaligned string representation of the CDR 1 region sequence"""
        return ''.join(self.cdr1_dict.values())

    @property
    def fw2_seq(self):
        """Unaligned string representation of the Framework 2 region sequence"""
        return ''.join(self.fw2_dict.values())

    @property
    def cdr2_seq(self):
        """Unaligned string representation of the CDR 2 region sequence"""
        return ''.join(self.cdr2_dict.values())

    @property
    def fw3_seq(self):
        """Unaligned string representation of the Framework 3 region sequence"""
        return ''.join(self.fw3_dict.values())

    @property
    def cdr3_seq(self):
        """Unaligned string representation of the CDR 3 region sequence"""
        return ''.join(self.cdr3_dict.values())

    @property
    def fw4_seq(self):
        """Unaligned string representation of the Framework 4 region sequence"""
        return ''.join(self.fw4_dict.values())


class Alignment:
    """Antibody chain alignment of two or more chains

    >>> from abnumber import Chain
    >>>
    >>> seq1 = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSSAKTTAP'
    >>> chain1 = Chain(seq1, scheme='imgt')
    >>>
    >>> seq2 = 'QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYDDYLDRWGQGTTLTVSSAKTTAP'
    >>> chain2 = Chain(seq2, scheme='imgt')
    >>> alignment = chain1.align(chain2)
    >>>
    >>> # Alignment can be sliced and iterated
    >>> for pos, (aa, bb) in alignment[:5]:
    >>>     print(pos, aa, bb)
    H1  Q Q
    H2  V V
    H3  Q Q
    H4  L L
    H5  Q V

    """
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
        """Format alignment to string

        >>> print(alignment.format())
        QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPS-RGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
        ||||.||||||.||||+|||||||||||.||||||||||||||||+||||||||.|.||||||||||||||||||||||||||.+|||||||||||||||||....||.|||||||||||
        QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYD--DYLDRWGQGTTLTVSS
                                 ^^^^^^^^                 ^^^^^^^^^                                      ^^^^^^^^^^^^
        >>> print(alignment.format(mark_identity=False, mark_cdrs=False))
        QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPS-RGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
        QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYD--DYLDRWGQGTTLTVSS

        :param mark_identity: Add BLAST style middle line showing identity (``|``), similar residue (``+``) or different residue (``.``)
        :param mark_cdrs: Add line highlighting CDR regions using ``^``
        :param colored: Highlight mutations using shell colors
        :return: formatted string
        """
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
        """Check if there is a mutation in the alignment or not"""
        return any(len(set(aas)) != 1 for aas in self.residues)


class Position:
    """Numbered position using a given numbering scheme

    Used as a key to store Position -> Amino acid information.

    Position objects are sortable according to the schema simply using ``sorted()``:
    """
    def __init__(self, chain_type: str, number: int, letter: str, scheme: str):
        _validate_chain_type(chain_type)
        self.chain_type: str = chain_type
        self.number: int = int(number)
        self.letter: str = letter
        self.scheme: str = scheme

    def __repr__(self):
        return f'{self.chain_type_prefix()}{self.number}{self.letter}({self.scheme})'

    def __str__(self):
        return f'{self.chain_type_prefix()}{self.number}{self.letter}'

    def __hash__(self):
        return self.__repr__().__hash__()

    def __eq__(self, other):
        return self.__repr__() == other.__repr__()

    def __lt__(self, other):
        assert self.chain_type == other.chain_type, f'Positions do not come from the same chain: {self}, {other}'
        assert self.scheme == other.scheme, 'Comparing positions in different schemes is not implemented'
        return self._sort_key() < other._sort_key()

    def chain_type_prefix(self):
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
        """Get string name of this position's region"""
        if self.scheme == 'imgt':
            return IMGT_POS_DICT[self.number]
        else:
            raise NotImplementedError(f'Not supported scheme: {self.scheme}')

    def is_in_cdr(self):
        """Check if given position is found in the CDR regions"""
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
    tail = sequence[end+1:]
    return aa_dict, chain_type, tail


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
