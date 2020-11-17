from collections import OrderedDict
from typing import Union, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import pandas as pd
from abnumber.exceptions import ChainParseError
import numpy as np

try:
    from anarci.anarci import anarci
except ImportError:
    # Only print the error without failing - required to import
    print('ANARCI module not available. Please install it separately or install AbNumber through Bioconda')
    print('See: https://abnumber.readthedocs.io/')
    sys.exit(1)
from Bio.Seq import Seq
import re

POS_REGEX = re.compile(r'([HL]?)(\d+)([A-Z]?)')


class Chain:
    """
    Antibody chain aligned to a chosen antibody numbering scheme

    :example:

    >>> from abnumber import Chain
    >>>
    >>> seq = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSSAKTTAPSVYPLA'
    >>> chain = Chain(seq, scheme='imgt')
    >>> chain
    QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
                             ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^

    Chain can be iterated:

    >>> for pos, aa in chain:
    >>>     print(pos, aa)
    H1  Q
    H2  V
    H3  Q
    H4  L
    H5  Q
    ...

    Chain can also be indexed and sliced using scheme numbering:

    >>> chain['5']
    'Q'
    >>> for pos, aa in chain['H2':'H5']:
    >>>     print(pos, aa)
    H2  V
    H3  Q
    H4  L
    H5  Q

    :param sequence: Unaligned string sequence
    :param name: Optional sequence identifier
    :param scheme: Numbering scheme: FIXME available numbering schemes
    :param allowed_species: ``None`` to allow all species, or one or more of: ``'human', 'mouse','rat','rabbit','rhesus','pig','alpaca'``
    :param aa_dict: Create Chain object directly from dictionary of region objects (internal use)
    :param tail: Constant region sequence
    :param species: Species as identified by ANARCI
    """

    def __init__(self, sequence, scheme, cdr_scheme=None, name=None, allowed_species=None, aa_dict=None, chain_type=None, tail=None, species=None):

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
            aa_dict, chain_type, tail, species = _anarci_align(sequence, scheme=scheme, allowed_species=allowed_species)

        _validate_chain_type(chain_type)

        self.name: str = name
        """User-provided sequence identifier"""
        self.chain_type: str = chain_type
        """Chain type as identified by ANARCI: ``H`` (heavy), ``K`` (kappa light) or ``L`` (lambda light)
        
        See also :meth:`Chain.is_heavy_chain` and :meth:`Chain.is_light_chain`.
        """
        self.scheme: str = scheme
        """Numbering scheme used to align the sequence"""
        self.cdr_scheme: str = cdr_scheme or scheme
        """Numbering scheme to be used for definition of CDR regions (same as ``scheme`` by default)"""
        self.tail: str = tail
        """Constant region sequence"""
        self.species: str = species
        """Species as identified by ANARCI"""

        self.fw1_dict = OrderedDict()
        self.cdr1_dict = OrderedDict()
        self.fw2_dict = OrderedDict()
        self.cdr2_dict = OrderedDict()
        self.fw3_dict = OrderedDict()
        self.cdr3_dict = OrderedDict()
        self.fw4_dict = OrderedDict()

        self._init_from_dict(aa_dict)

    def _init_from_dict(self, aa_dict):
        if self.scheme not in SUPPORTED_SCHEMES:
            raise NotImplementedError(f'Scheme "{self.scheme}" is not supported. Available schemes: {", ".join(SUPPORTED_SCHEMES)}')
        if self.cdr_scheme in ['aho']:
            raise ValueError('CDR regions are not defined for AHo, '
                             'you need to specify cdr_scheme="chothia" or another scheme for CDR extraction.')
        # list of region start positions
        borders = SCHEME_BORDERS[self.scheme] if self.scheme in SCHEME_BORDERS else SCHEME_BORDERS[f'{self.scheme}_{self.chain_type}']

        regions_list = [self.fw1_dict, self.cdr1_dict, self.fw2_dict, self.cdr2_dict, self.fw3_dict, self.cdr3_dict, self.fw4_dict]
        region_idx = 0

        for pos in sorted(aa_dict.keys()):
            while pos.number >= borders[region_idx]:
                region_idx += 1
            aa = aa_dict[pos].upper().strip()
            if aa in ['*', '-', None, '']:
                continue
            regions_list[region_idx][pos] = aa

    def __repr__(self):
        return self.format()

    def __str__(self):
        return self.seq

    def __iter__(self):
        yield from self.positions.items().__iter__()

    def __getitem__(self, item):
        if isinstance(item, slice):
            if item.step is not None and item.step != 1:
                raise IndexError(f'Slicing with step != 1 is not implemented, got: {item}')
            return self.slice(start=item.start, stop=item.stop)
        pos = self.parse_position(item)
        return self.positions[pos]

    def __len__(self):
        return len(self.positions)

    @classmethod
    def to_fasta(cls, chains, path_or_fd, keep_tail=False):
        if isinstance(chains, Chain):
            records = chains.to_seq_record(keep_tail=keep_tail)
        else:
            records = (chain.to_seq_record(keep_tail=keep_tail) for chain in chains)
        return SeqIO.write(records, path_or_fd, 'fasta-2line')

    def to_seq_record(self, keep_tail=False):
        if not self.name:
            raise ValueError('Name needs to be present to convert to a SeqRecord')
        seq = Seq(self.seq + self.tail if keep_tail else self.seq)
        return SeqRecord(seq, id=self.name, description='')

    @classmethod
    def to_anarci_csv(cls, chains: List['Chain'], path):
        df = cls.to_dataframe(chains)
        df.to_csv(path)

    @classmethod
    def to_dataframe(cls, chains: List['Chain']):
        df = pd.DataFrame([chain.to_series() for chain in chains]).fillna('-')
        df.index.name = 'Id'

        # Each chain can have a different set of positions
        # so we need to sort the columns to make sure they are in the right order
        # this is using the correct Position sorting
        prop_columns = [c for c in df.columns if not isinstance(c, Position)]
        position_columns = sorted([c for c in df.columns if isinstance(c, Position)])
        df = df[prop_columns + position_columns]

        # Finally convert position columns to string
        df.columns = df.columns.map(lambda pos: pos.format(chain_type=False))

        return df

    def to_series(self):
        props = {
            'chain_type': self.chain_type,
            'species': self.species
        }
        return pd.Series({**props, **self.positions}, name=self.name)

    @classmethod
    def from_series(cls, series, scheme) -> 'Chain':
        chain_type = series['chain_type']
        species = series['species']
        position_index = [c for c in series.index if c[:1].isnumeric()]
        aa_dict = {Position.from_string(pos, chain_type=chain_type, scheme=scheme): aa
                   for pos, aa in series[position_index].items() if aa != '-' and not pd.isna(aa)}
        return cls(sequence=None, aa_dict=aa_dict, name=series.name, scheme=scheme, chain_type=chain_type, species=species)

    @classmethod
    def from_anarci_csv(cls, path, scheme, as_series=False) -> Union[List['Chain'], pd.Series]:
        df = pd.read_csv(path, index_col=0)
        return cls.from_dataframe(df, scheme=scheme, as_series=as_series)

    @classmethod
    def from_dataframe(cls, df, scheme, as_series=False) -> Union[List['Chain'], pd.Series]:
        chains = [cls.from_series(series, scheme=scheme) for i, series in df.iterrows()]
        if as_series:
            return pd.Series(chains, index=[c.name for c in chains])
        return chains

    def format(self, method='wide', **kwargs):
        """Format sequence to string

        :param method: use ``"wide"`` for :meth:`Chain.format_wide` or ``"tall"`` for :meth:`Chain.format_tall()`
        :return: formatted string
        """
        if method == 'wide':
            return self.format_wide(**kwargs)
        elif method == 'tall':
            return self.format_tall(**kwargs)
        raise ValueError(f'Use method="wide" or method="tall", unknown method: "{method}"')
    
    def print(self, method='wide'):
        """Print string representation using :meth:`Chain.format`

        By default, produces "wide" format with sequence on first line and CDR regions higlighted with ``^`` on second line:

        >>> chain.print()
        QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
                                 ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^

        :param method: use ``"wide"`` for :meth:`Chain.format_wide` or ``"tall"`` for :meth:`Chain.format_tall()`
        """
        print(self.format(method=method))

    def format_tall(self, columns=5):
        """Create string with one position per line, showing position numbers and amino acids

        :return: formatted string
        """
        height = int(np.ceil(len(self) / columns))
        rows = [''] * height
        for column, start in enumerate(range(0, len(self), height)):
            chain_slice = self.raw[start:start+height]
            for row, (pos, aa) in enumerate(chain_slice):
                rows[row] = rows[row].ljust(column * 15)
                pos_format = (pos.get_region() + ' ' if pos.is_in_cdr() else '') + pos.format()
                rows[row] += f'{pos_format.rjust(9)} {aa}'

        return '\n'.join(rows)

    def print_tall(self, columns=5):
        """Print string representation using :meth:`Chain.format_tall`

        >>> chain.print_tall()
        fw1 H1    Q
        fw1 H2    V
        fw1 H3    Q
        fw1 H4    L
        fw1 H5    Q
        fw1 H6    Q
        fw1 H7    S
        ...
        """
        print(self.format_tall(columns=columns))

    def format_wide(self):
        """Create string with sequence on first line and CDR regions higlighted with `^` on second line

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

    def print_wide(self):
        """Print string representation using :meth:`Chain.format_wide`

        >>> chain.print_wide()
        QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
                                 ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^
        """
        print(self.print_wide())

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
        return Alignment(shared_pos, residues, chain_type=self.chain_type, scheme=self.scheme)

    def clone(self, replace_seq: str = None):
        """Create a copy of this chain, optionally with a replacement sequence
        
        :param replace_seq: Optional replacement sequence, needs to be the same length
        :return: new Chain object
        """
        return self.slice(replace_seq=replace_seq)

    def slice(self, replace_seq: str = None, start: Union[str, int, 'Position'] = None,
              stop: Union[str, int, 'Position'] = None, stop_inclusive: bool = True, allow_raw: bool = False):
        """Create a slice of this chain, optionally with a replacement sequence that is placed into the same numbering

        You can also slice directly using ``chain['111':'112A']`` or ``chain.raw[10:20]``.

        :param replace_seq: Optional replacement sequence, needs to be the same length
        :param start: Optional slice start position (inclusive), :class:`Position` or string (e.g. '111A')
        :param stop: Optional slice stop position (inclusive), :class:`Position` or string (e.g. '112A')
        :param stop_inclusive: Include stop position in slice
        :param allow_raw: Allow unaligned numeric indexing from 0 to length of sequence - 1
        :return: new Chain object
        """
        aa_dict = {}
        positions = self.positions
        if replace_seq is not None:
            assert len(replace_seq) == len(positions), 'Sequence needs to be the same length'

        start = self.parse_position(start, allow_raw=allow_raw) if start is not None else None
        stop = self.parse_position(stop, allow_raw=allow_raw) if stop is not None else None

        for i, (pos, aa) in enumerate(positions.items()):
            if start is not None and pos < start:
                continue
            if stop is not None and (pos > stop or (not stop_inclusive and pos >= stop)):
                break
            aa_dict[pos] = replace_seq[i] if replace_seq is not None else aa

        return Chain(sequence=None, aa_dict=aa_dict, name=self.name, scheme=self.scheme, chain_type=self.chain_type, species=self.species)

    def graft_cdrs_onto(self, other: 'Chain', name: str = None) -> 'Chain':
        """Graft CDRs from this Chain onto another chain

        :param other: Chain to graft CDRs into
        :param name: Name of new Chain. If not provided, use name of this chain.
        :return: Chain with CDRs grafted from this chain and frameworks from the given chain
        """
        assert self.scheme == other.scheme, \
            f'Sequences need to have the same numbering scheme, got {self.scheme} and {other.scheme}'
        assert self.chain_type == other.chain_type, \
            f'Sequences need to have the same chain type, got {self.chain_type} and {other.chain_type}'

        grafted_dict = {}
        for (self_region, self_dict), (other_region, other_dict) in zip(self.regions.items(), other.regions.items()):
            assert self_region == other_region
            if self_region.lower().startswith('cdr'):
                grafted_dict.update(self_dict)
            else:
                grafted_dict.update(other_dict)

        return Chain(sequence=None, aa_dict=grafted_dict, name=name or self.name, chain_type=self.chain_type, scheme=self.scheme)

    def parse_position(self, position: Union[int, str, 'Position'], allow_raw=False):
        """Create :class:`Position` key object from string

        :param position: Numeric or string position representation
        :param allow_raw: Also allow unaligned numeric (int) indexing from 0 to length of sequence - 1
        :return: new Position object
        """
        if isinstance(position, str):
            return Position.from_string(position, chain_type=self.chain_type, scheme=self.scheme)
        if isinstance(position, Position):
            return position
        try:
            position = int(position)
        except TypeError:
            raise IndexError(f'Invalid position key, expected Position, string or integer, got {type(position)}: "{position}"')
        if not allow_raw:
            raise IndexError("Use chain.raw[i] for raw numeric indexing or pass allow_raw=True. "
                             "For named position indexing, use string (e.g. chain['111A'] or chain['H111A'])")
        if position >= len(self.positions):
            return None
        return list(self.positions.keys())[position]

    @property
    def raw(self):
        """Access raw representation of this chain to allow unaligned numeric indexing and slicing

        >>> # String numbering is based on schema numbering
        >>> chain['1']
        'QVQLQQSGAE'
        >>> # Numbering of ``chain.raw`` starts at 0
        >>> chain.raw[0]
        'QVQLQQSGAE'
        >>> # Slicing with string is based on schema numbering, the end is inclusive
        >>> chain['1':'10']
        'QVQLQQSGAE'
        >>> # Slicing with ``chain.raw`` starts at 0, the end is exclusive (Python style)
        >>> chain.raw[0:10]
        'QVQLQQSGAE'

        :return: Raw chain accessor that can be sliced or indexed to produce a new :class:`Chain` object
        """
        return RawChainAccessor(self)

    @property
    def regions(self):
        """Dictionary of region dictionaries

        Region is an uppercase string, one of: ``"FW1", "CDR1", "FW2", "CDR2", "FW3", "CDR3", "FW4"``

        :return: Dictionary of Region name -> Dictionary of (:class:`Position` -> Amino acid)
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
        """Dictionary of :class:`Position` -> Amino acid"""
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


class RawChainAccessor:
    def __init__(self, chain: Chain):
        self.chain = chain

    def __getitem__(self, item):
        if isinstance(item, slice):
            if item.step is not None and item.step != 1:
                raise IndexError(f'Slicing with step != 1 is not implemented, got: {item}')
            if item.start is not None and not isinstance(item.start, int):
                raise IndexError(f'Expected int start index for chain.raw, got {type(item.start)}: {item.start}')
            if item.stop is not None and not isinstance(item.stop, int):
                raise IndexError(f'Expected int end index for chain.raw, got {type(item.stop)}: {item.stop}')
            return self.chain.slice(start=item.start, stop=item.stop, stop_inclusive=False, allow_raw=True)
        if not isinstance(item, int):
            raise IndexError(f'Expected int indexing for chain.raw, got {type(item)}: {item}')
        pos = self.chain.parse_position(item, allow_raw=True)
        return self.chain[pos]


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

    Alignment can be sliced and iterated:

    >>> for pos, (aa, bb) in alignment[:'5']:
    >>>     print(pos, aa, bb)
    H1  Q Q
    H2  V V
    H3  Q Q
    H4  L L
    H5  Q V
    ...

    """
    def __init__(self, positions, residues, scheme, chain_type):
        assert isinstance(positions, list), 'Expected list of positions and residues. ' \
                                            'Use chain.align(other) to create an alignment.'
        assert len(positions) == len(residues)
        self.positions = positions
        self.residues = residues
        self.scheme = scheme
        self.chain_type = chain_type
        self._zipped = list(zip(self.positions, self.residues))

    def __repr__(self):
        return self.format()

    def __iter__(self):
        yield from self._zipped.__iter__()

    def __len__(self):
        return len(self.positions)

    def __getitem__(self, item):
        if isinstance(item, slice):
            if item.step is not None and item.step != 1:
                raise IndexError(f'Slicing with step != 1 is not implemented, got: {item}')
            return self.slice(start=item.start, stop=item.stop)
        pos = self.parse_position(item)
        raw_pos = self.positions.index(pos)
        return self.residues[raw_pos]

    def slice(self, start: Union[str, int, 'Position'] = None, stop: Union[str, int, 'Position'] = None,
              stop_inclusive: bool = True, allow_raw: bool = False):
        """Create a slice of this alignment

        You can also slice directly using ``alignment['111':'112A']`` or ``alignment.raw[10:20]``.

        :param start: Slice start position (inclusive), :class:`Position` or string (e.g. '111A')
        :param stop: Slice stop position (inclusive), :class:`Position` or string (e.g. '112A')
        :param stop_inclusive: Include stop position in slice
        :param allow_raw: Allow unaligned numeric indexing from 0 to length of sequence - 1
        :return: new sliced Alignment object
        """

        start = self.parse_position(start, allow_raw=allow_raw) if start is not None else None
        stop = self.parse_position(stop, allow_raw=allow_raw) if stop is not None else None

        new_positions = []
        new_residues = []
        for pos, residues in zip(self.positions, self.residues):
            if start is not None and pos < start:
                continue
            if stop is not None and (pos > stop or (not stop_inclusive and pos >= stop)):
                break
            new_positions.append(pos)
            new_residues.append(residues)

        return Alignment(positions=new_positions, residues=new_residues, scheme=self.scheme, chain_type=self.chain_type)

    def parse_position(self, position: Union[int, str, 'Position'], allow_raw=False):
        """Create :class:`Position` key object from string

        :param position: Numeric or string position representation
        :param allow_raw: Also allow unaligned numeric (int) indexing from 0 to length of sequence - 1
        :return: new Position object
        """
        if isinstance(position, str):
            return Position.from_string(position, chain_type=self.chain_type, scheme=self.scheme)
        if isinstance(position, Position):
            return position
        try:
            position = int(position)
        except TypeError:
            raise IndexError(f'Invalid position key, expected Position, string or integer, got {type(position)}: "{position}"')
        if not allow_raw:
            raise IndexError("Use chain.raw[i] for raw numeric indexing or pass allow_raw=True. "
                             "For named position indexing, use string (e.g. chain['111A'] or chain['H111A'])")
        if position >= len(self.positions):
            return None
        return self.positions[position]

    def format(self, mark_identity=True, mark_cdrs=True):
        """Format alignment to string

        :param mark_identity: Add BLAST style middle line showing identity (``|``), similar residue (``+``) or different residue (``.``)
        :param mark_cdrs: Add line highlighting CDR regions using ``^``
        :return: formatted string
        """

        def _identity_symbol(a, b):
            return '|' if a == b else ('+' if is_similar_residue(a, b) else '.')

        lines = []
        for i in range(len(self.residues[0])):
            if mark_identity and i != 0:
                lines.append(''.join(_identity_symbol(aas[i], aas[i-1]) for pos, aas in self))
            lines.append(''.join(aas[i] for pos, aas in self))
        if mark_cdrs:
            lines.append(''.join('^' if pos.is_in_cdr() else ' ' for pos, aas in self))
        return '\n'.join(lines)

    def print(self, mark_identity=True, mark_cdrs=True):
        """Print string representation of alignment created using :meth:`Alignment.format`

        >>> alignment.print()
        QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPS-RGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
        ||||.||||||.||||+|||||||||||.||||||||||||||||+||||||||.|.||||||||||||||||||||||||||.+|||||||||||||||||....||.|||||||||||
        QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYD--DYLDRWGQGTTLTVSS
                                 ^^^^^^^^                 ^^^^^^^^^                                      ^^^^^^^^^^^^
        >>> alignment.print(mark_identity=False, mark_cdrs=False)
        QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPS-RGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
        QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYD--DYLDRWGQGTTLTVSS

        :param mark_identity: Add BLAST style middle line showing identity (``|``), similar residue (``+``) or different residue (``.``)
        :param mark_cdrs: Add line highlighting CDR regions using ``^``
        """
        print(self.format(mark_identity=mark_identity, mark_cdrs=mark_cdrs))

    def has_mutation(self):
        """Check if there is a mutation in the alignment or not"""
        return any(len(set(aas)) != 1 for aas in self.residues)

    def num_mutations(self):
        """Get number of mutations (positions with more than one type of residue)"""
        return sum(len(set(aas)) != 1 for aas in self.residues)

    @property
    def raw(self):
        """Access raw representation of this alignment to allow unaligned numeric indexing and slicing

        >>> # Numbering of ``chain.raw`` starts at 0
        >>> alignment.raw[0]
        'H1'
        >>> # Slicing with string is based on schema numbering, the end is inclusive
        >>> chain['1':'10']
        'QVQLQQSGAE'
        >>> # Slicing with ``chain.raw`` starts at 0, the end is exclusive (Python style)
        >>> chain.raw[0:10]
        'QVQLQQSGAE'
        :return: Raw alignment accessor that can be sliced or indexed to produce a new :class:`Alignment` object
        """
        return RawAlignmentAccessor(self)


class RawAlignmentAccessor:
    def __init__(self, alignment: Alignment):
        self.alignment = alignment

    def __getitem__(self, item):
        if isinstance(item, slice):
            if item.step is not None and item.step != 1:
                raise IndexError(f'Slicing with step != 1 is not implemented, got: {item}')
            if item.start is not None and not isinstance(item.start, int):
                raise IndexError(f'Expected int start index for alignment.raw, got {type(item.start)}: {item.start}')
            if item.stop is not None and not isinstance(item.stop, int):
                raise IndexError(f'Expected int end index for alignment.raw, got {type(item.stop)}: {item.stop}')
            return self.alignment.slice(start=item.start, stop=item.stop, stop_inclusive=False, allow_raw=True)
        if not isinstance(item, int):
            raise IndexError(f'Expected int indexing for alignment.raw, got {type(item)}: {item}')
        pos = self.alignment.parse_position(item, allow_raw=True)
        return self.alignment[pos]


class Position:
    """Numbered position using a given numbering scheme

    Used as a key to store Position -> Amino acid information.

    Position objects are sortable according to the schema simply using ``sorted()``.
    """
    def __init__(self, chain_type: str, number: int, letter: str, scheme: str):
        _validate_chain_type(chain_type)
        self.chain_type: str = chain_type
        self.number: int = int(number)
        self.letter: str = letter.strip() if letter.strip() else ''
        self.scheme: str = scheme

    def __repr__(self):
        return f'{self.chain_type_prefix()}{self.number}{self.letter} ({self.scheme})'

    def __str__(self):
        return self.format()

    def format(self, chain_type=True, region=False, rjust=False, ljust=False, fillchar=' '):
        formatted = f'{self.number}{self.letter}'
        if chain_type:
            formatted = f'{self.chain_type_prefix()}{formatted}'
        if region:
            formatted = f'{self.get_region()} {formatted}'
        just = 4 + 1*int(chain_type) + 5*int(region)
        if rjust:
            formatted = formatted.rjust(just, fillchar)
        if ljust:
            formatted = formatted.ljust(just, fillchar)
        return formatted

    def __hash__(self):
        return self.__repr__().__hash__()

    def __eq__(self, other):
        return self.__repr__() == other.__repr__()

    def __ge__(self, other):
        return self == other or self > other

    def __le__(self, other):
        return self == other or self < other

    def __lt__(self, other):
        if not isinstance(other, Position):
            raise TypeError(f'Cannot compare Position object with {type(other)}: {other}')
        assert self.is_heavy_chain() == other.is_heavy_chain(), f'Positions do not come from the same chain: {self}, {other}'
        assert self.scheme == other.scheme, 'Comparing positions in different schemes is not implemented'
        return self._sort_key() < other._sort_key()

    def chain_type_prefix(self):
        if self.chain_type == 'H':
            return 'H'
        if self.chain_type in ['K', 'L']:
            return 'L'
        raise NotImplementedError(f'Unknown chain type "{self.chain_type}"')

    def _sort_key(self):
        letter_ord = ord(self.letter) if self.letter else 0
        if self.scheme == 'imgt':
            if self.number == 112:
                # position 112 is sorted in reverse
                letter_ord = -letter_ord
        elif self.scheme in ['chothia', 'kabat', 'aho']:
            # all letters are sorted alphabetically for these schemes
            pass
        else:
            raise NotImplementedError(f'Cannot compare positions of scheme: {self.scheme}')
        return self.number, letter_ord

    def get_region(self):
        """Get string name of this position's region

        :return: uppercase string, one of: ``"FW1", "CDR1", "FW2", "CDR2", "FW3", "CDR3", "FW4"``
        """
        if self.scheme in SCHEME_POSITION_TO_REGION:
            regions = SCHEME_POSITION_TO_REGION[self.scheme]
        else:
            regions = SCHEME_POSITION_TO_REGION[f'{self.scheme}_{self.chain_type}']
        return regions[self.number]

    def is_in_cdr(self):
        """Check if given position is found in the CDR regions"""
        # FIXME
        return self.get_region().lower().startswith('cdr')

    @classmethod
    def from_string(cls, position, chain_type, scheme):
        match = POS_REGEX.match(position.upper())
        expected_chain_prefix = 'H' if chain_type == 'H' else 'L'
        if match is None:
            raise IndexError(f'Expected position format chainNumberLetter '
                             f'(e.g. "{expected_chain_prefix}112A" or "112A"), got: "{position}"')
        chain_prefix, number, letter = match.groups()
        number = int(number)
        if chain_prefix and expected_chain_prefix != chain_prefix:
            raise IndexError(f'Use no prefix or "{expected_chain_prefix}" prefix for "{chain_type}" chain. '
                             f'Got: "{chain_prefix}".')
        return cls(chain_type=chain_type, number=number, letter=letter, scheme=scheme)

    def is_heavy_chain(self):
        return self.chain_type == 'H'

    def is_light_chain(self):
        return self.chain_type in 'KL'


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
        raise ChainParseError(f'Variable chain sequence not recognized: "{sequence}"')
    if len(numbered) != 1:
        raise NotImplementedError(f'Unsupported: Multiple ANARCI domains found in sequence: "{sequence}"')
    positions, start, end = numbered[0]
    chain_type = ali[0]['chain_type']
    species = ali[0]['species']
    aa_dict = {Position(chain_type=chain_type, number=num, letter=letter, scheme=scheme): aa for (num, letter), aa in
               positions if aa != '-'}
    tail = sequence[end+1:]
    return aa_dict, chain_type, tail, species


# Based on positive score in Blosum62
SIMILAR_PAIRS = {'AA', 'AS', 'CC', 'DD', 'DE', 'DN', 'ED', 'EE', 'EK', 'EQ', 'FF', 'FW', 'FY', 'GG', 'HH', 'HN', 'HY',
                 'II', 'IL', 'IM', 'IV', 'KE', 'KK', 'KQ', 'KR', 'LI', 'LL', 'LM', 'LV', 'MI', 'ML', 'MM', 'MV', 'ND',
                 'NH', 'NN', 'NS', 'PP', 'QE', 'QK', 'QQ', 'QR', 'RK', 'RQ', 'RR', 'SA', 'SN', 'SS', 'ST', 'TS', 'TT',
                 'VI', 'VL', 'VM', 'VV', 'WF', 'WW', 'WY', 'YF', 'YH', 'YW', 'YY'}


def is_similar_residue(a, b):
    if a == '-' or b == '-':
        return a == b
    return a+b in SIMILAR_PAIRS


SUPPORTED_SCHEMES = ['imgt', 'aho', 'chothia', 'kabat']

SCHEME_BORDERS = {
               # Start coordinates
               # CDR1, FW2, CDR2, FW3, CDR3, FW4
         'imgt': [27,  39,  56,   66,  105,  118, 129],
    'chothia_H': [26,  33,  52,   57,  95,   103, 114],
    'chothia_K': [24,  35,  50,   57,  89,    98, 108],
    'chothia_L': [24,  35,  50,   57,  89,    98, 108],
      'kabat_H': [31,  36,  50,   66,  95,   103, 114],
      'kabat_K': [24,  35,  50,   57,  89,    98, 108],
      'kabat_L': [24,  35,  50,   57,  89,    98, 108],
}

# { scheme -> { region -> list of position numbers } }
SCHEME_REGIONS = {
    scheme: {
        'FW1': list(range(1, borders[0])),
        'CDR1': list(range(borders[0], borders[1])),
        'FW2': list(range(borders[1], borders[2])),
        'CDR2': list(range(borders[2], borders[3])),
        'FW3': list(range(borders[3], borders[4])),
        'CDR3': list(range(borders[4], borders[5])),
        'FW4': list(range(borders[5], borders[6])),
    } for scheme, borders in SCHEME_BORDERS.items()
}

# { scheme -> { position number -> region } }
SCHEME_POSITION_TO_REGION = {
    scheme: {pos_num: region for region, positions in regions.items() for pos_num in positions} \
    for scheme, regions in SCHEME_REGIONS.items()
}