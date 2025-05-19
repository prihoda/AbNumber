import warnings
from collections import OrderedDict
from typing import Union, List, Generator, Tuple
import pandas as pd

from abnumber.alignment import Alignment
from abnumber.common import _anarci_align, _validate_chain_type, SUPPORTED_SCHEMES, SUPPORTED_CDR_DEFINITIONS, \
    is_integer, SCHEME_BORDERS, _get_unique_chains
from abnumber.exceptions import ChainParseError, MultipleDomainsChainParseError
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from abnumber.position import Position


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
    :param scheme: Numbering scheme: One of ``imgt``, ``chothia``, ``kabat``, ``aho``
    :param cdr_definition: Numbering scheme to be used for definition of CDR regions. Same as ``scheme`` by default.
                           One of ``imgt``, ``chothia``, ``kabat``, ``north``. Required for ``aho``.
    :param assign_germline: Assign germline name using ANARCI based on best sequence identity
    :param allowed_species: Allowed species for germline assignment. Use ``None`` to allow all species, or one or more of: ``'human', 'mouse','rat','rabbit','rhesus','pig','alpaca'``
    :param aa_dict: (Internal use only) Create Chain object directly from dictionary of region objects (internal use)
    :param tail: (Internal use only) Constant region sequence
    :param species: (Internal use only) Species as identified by ANARCI
    :param germline: (Internal use only) Germline as identified by ANARCI
    :param use_anarcii: Use ANARCII (2.0) for numbering, otherwise use ANARCI (1.0)
    :param anarcii_args: Keyword arguments for Anarcii()
    """

    def __init__(self, sequence, scheme, cdr_definition=None, name=None, assign_germline=False, allowed_species=None, use_anarcii=False, anarcii_args=None, **kwargs):
        aa_dict = kwargs.pop('aa_dict', None)
        chain_type = kwargs.pop('chain_type', None)
        tail = kwargs.pop('tail', None)
        species = kwargs.pop('species', None)
        v_gene = kwargs.pop('v_gene', None)
        j_gene = kwargs.pop('j_gene', None)
        if isinstance(allowed_species, str):
            allowed_species = [allowed_species]
        if len(kwargs):
            raise TypeError(f'Argument not recognized: {", ".join(kwargs)}')
        if aa_dict is not None:
            if sequence is not None:
                raise ChainParseError('Only one of aa_dict= and sequence= can be provided')
            assert isinstance(aa_dict, dict), f'Expected dict, got: {type(aa_dict)}'
            assert tail is not None
            assert chain_type is not None
        else:
            if sequence is None:
                raise ChainParseError('Expected sequence, got None')
            if isinstance(sequence, list):
                raise ChainParseError('Expected string or Seq, got list. Please use Chain.batch() to parse multiple sequences')
            if not isinstance(sequence, str) and not isinstance(sequence, Seq):
                raise ChainParseError(f'Expected string or Seq, got {type(sequence)}: {sequence}')
            if '-' in sequence:
                raise ChainParseError(f'Please provide an unaligned sequence, got: {sequence}')
            if chain_type is not None:
                raise ChainParseError('Do not use chain_type= when providing sequence=, it will be inferred automatically')
            if tail is not None:
                raise ChainParseError('Do not use tail= when providing sequence=, it will be inferred automatically')
            if isinstance(sequence, Seq):
                sequence = str(sequence)
            results = _anarci_align(
                [sequence],
                scheme=scheme,
                allowed_species=allowed_species,
                assign_germline=assign_germline,
                use_anarcii=use_anarcii,
                anarcii_args=anarcii_args,
            )[0]
            if not results:
                raise ChainParseError(f'Variable chain sequence not recognized: "{sequence}"')
            if len(results) > 1:
                warnings.warn('Use Chain.multiple_domains(seq) to parse ScFvs and other sequences with multiple antibody domains')
                raise MultipleDomainsChainParseError(f'Found {len(results)} antibody domains in sequence: "{sequence}"')
            aa_dict, chain_type, tail, species, v_gene, j_gene = results[0]

        _validate_chain_type(chain_type)

        self.name: str = name
        """User-provided sequence identifier"""
        self.chain_type: str = chain_type
        """Chain type as identified by ANARCI: ``H`` (heavy), ``K`` (kappa light) or ``L`` (lambda light)
        
        See also :meth:`Chain.is_heavy_chain` and :meth:`Chain.is_light_chain`.
        """
        self.scheme: str = scheme
        """Numbering scheme used to align the sequence"""
        self.cdr_definition: str = cdr_definition or scheme
        """Numbering scheme to be used for definition of CDR regions (same as ``scheme`` by default)"""
        self.tail: str = tail
        """Constant region sequence"""
        self.species: str = species
        """Species as identified by ANARCI"""
        self.v_gene: str = v_gene
        """V gene germline as identified by ANARCI (if assign_germline is True)"""
        self.j_gene: str = j_gene
        """J gene germline as identified by ANARCI (if assign_germline is True)"""

        self.fr1_dict = OrderedDict()
        self.cdr1_dict = OrderedDict()
        self.fr2_dict = OrderedDict()
        self.cdr2_dict = OrderedDict()
        self.fr3_dict = OrderedDict()
        self.cdr3_dict = OrderedDict()
        self.fr4_dict = OrderedDict()

        self._init_from_dict(aa_dict, allowed_species=allowed_species, use_anarcii=use_anarcii)

    def _init_from_dict(self, aa_dict, allowed_species, use_anarcii=False):
        if self.scheme not in SUPPORTED_SCHEMES:
            raise NotImplementedError(f'Scheme "{self.scheme}" is not supported. Available schemes: {", ".join(SUPPORTED_SCHEMES)}')
        if self.cdr_definition in ['aho']:
            raise ValueError('CDR regions are not defined for AHo, '
                             'you need to specify cdr_definition="chothia" or another scheme for CDR extraction.')
        if self.cdr_definition not in SUPPORTED_CDR_DEFINITIONS:
            raise NotImplementedError(f'CDR definition "{self.cdr_definition}" is not supported. Available definitions: {", ".join(SUPPORTED_CDR_DEFINITIONS)}')
        # list of region start positions
        borders = SCHEME_BORDERS[self.cdr_definition] if self.cdr_definition in SCHEME_BORDERS else SCHEME_BORDERS[f'{self.cdr_definition}_{self.chain_type}']

        regions_list = [self.fr1_dict, self.cdr1_dict, self.fr2_dict, self.cdr2_dict, self.fr3_dict, self.cdr3_dict, self.fr4_dict]
        region_idx = 0

        sorted_positions = sorted(aa_dict.keys())

        cdr_definition_ready = True
        for pos in sorted_positions:
            assert pos.scheme == self.scheme, f'Schemes of provided position ({pos.scheme}) does not match Chain scheme ({self.scheme})'
            if pos.cdr_definition != self.cdr_definition:
                cdr_definition_ready = False

        if cdr_definition_ready:
            combined_aa_dict = aa_dict
        else:
            seq = ''.join(aa_dict[pos] for pos in sorted_positions)
            renumbered_aa_dict = _anarci_align(
                [seq],
                scheme=self.cdr_definition if self.cdr_definition != 'north' else 'chothia',
                allowed_species=allowed_species,
                use_anarcii=use_anarcii
            )[0][0][0]
            cdr_definition_positions = [pos.number for pos in sorted(renumbered_aa_dict.keys())]
            combined_aa_dict = {}
            for orig_pos, cdr_definition_position in zip(sorted_positions, cdr_definition_positions):
                aa = aa_dict[orig_pos]
                pos = orig_pos.copy()
                pos.set_cdr_definition(self.cdr_definition, cdr_definition_position)
                combined_aa_dict[pos] = aa

        for pos in sorted(combined_aa_dict.keys()):
            assert isinstance(pos, Position), f'Expected Position object, got {type(pos)}: {pos}'
            aa = combined_aa_dict[pos].upper().strip()
            if aa in [None, '*', '-', '', '.']:
                continue
            while pos.cdr_definition_position >= borders[region_idx]:
                region_idx += 1
            regions_list[region_idx][pos] = aa

    @classmethod
    def batch(cls, seq_dict: dict, scheme: str, cdr_definition=None, assign_germline=False, allowed_species=None, multiple_domains=False, use_anarcii=False, anarcii_args=None):
        """Create multiple Chain objects from dict of sequences

        :param seq_dict: Dictionary of sequence strings, keys are sequence identifiers
        :param scheme: Numbering scheme to align the sequences
        :param cdr_definition: Numbering scheme to be used for definition of CDR regions. Same as ``scheme`` by default.
        :param assign_germline: Assign germline name using ANARCI based on best sequence identity
        :param allowed_species: Allowed species for germline assignment. Use ``None`` to allow all species, or one or more of: ``'human', 'mouse','rat','rabbit','rhesus','pig','alpaca'``
        :param multiple_domains: Allow parsing multiple domains in a sequence - return dict name -> list of one or more Chain items - will use ANARCI 1.0!
        :param use_anarcii: Use ANARCII (2.0) for numbering, otherwise use ANARCI (1.0)
        :param anarcii_args: Keyword arguments for Anarcii()
        :return: tuple with (dict of Chain objects, dict of error strings)
        """
        assert isinstance(seq_dict, dict), f'Expected dictionary of sequences, got: {type(seq_dict).__name__}'
        if not seq_dict:
            return {}, {}
        names = list(seq_dict.keys())
        seq_list = list(seq_dict.values())
        all_results = _anarci_align(
            seq_list,
            scheme=scheme,
            allowed_species=allowed_species,
            assign_germline=assign_germline,
            use_anarcii=use_anarcii,
            anarcii_args=anarcii_args,
        )
        names = names or ([None] * len(seq_list))
        chains = {}
        errors = {}
        for sequence, results, name in zip(seq_list, all_results, names):
            if not results:
                errors[name] = f'Variable chain sequence not recognized: "{sequence}"'
            elif len(results) > 1 and not multiple_domains:
                warnings.warn('Use multiple_domains=True to allow parsing ScFvs and other sequences with multiple antibody domains')
                errors[name] = f'Found {len(results)} antibody domains: "{sequence}"'
            else:
                found_chains = []
                for aa_dict, chain_type, tail, species, v_gene, j_gene in results:
                    try:
                        found_chains.append(Chain(
                            sequence=None,
                            aa_dict=aa_dict,
                            name=name,
                            scheme=scheme,
                            chain_type=chain_type,
                            cdr_definition=cdr_definition,
                            tail=tail,
                            species=species,
                            v_gene=v_gene,
                            j_gene=j_gene
                        ))
                    except Exception as e:
                        errors[name] = f'Unexpected error creating Chain: {e}'
                chains[name] = found_chains if multiple_domains else found_chains[0]
        return chains, errors

    @classmethod
    def multiple_domains(cls, sequence: str, scheme: str, cdr_definition=None, name=None, assign_germline=False, allowed_species=None, use_anarcii=False) -> 'Chain':
        """Parse multi-domain sequence into a list of Chain objects

        :param sequence: Unaligned string sequence
        :param scheme: Numbering scheme to align the sequences
        :param cdr_definition: Numbering scheme to be used for definition of CDR regions. Same as ``scheme`` by default.
        :param name: Optional sequence identifier
        :param assign_germline: Assign germline name using ANARCI based on best sequence identity
        :param allowed_species: Allowed species for germline assignment. Use ``None`` to allow all species, or one or more of: ``'human', 'mouse','rat','rabbit','rhesus','pig','alpaca'``
        :param use_anarcii: Use ANARCII (2.0) for numbering, otherwise use ANARCI (1.0)
        :return: tuple with (dict of Chain objects, dict of error strings)
        """
        chains, errors = cls.batch({name: sequence}, scheme=scheme, cdr_definition=cdr_definition, assign_germline=assign_germline, allowed_species=allowed_species, multiple_domains=True, use_anarcii=use_anarcii)
        if error := errors.get(name):
            raise ChainParseError(error)
        return chains[name]


    def __repr__(self):
        return self.format()

    def __str__(self):
        return self.seq

    def __iter__(self):
        yield from self.positions.items().__iter__()

    def __contains__(self, item):
        position = self._parse_position(item)
        assert position.scheme == self.scheme
        return position in self.positions

    def __setitem__(self, item, value):
        """Mutate single position by directly modifying this Chain object"""
        assert isinstance(value, str), f'Expected amino acid string, got: {type(value)} {value}'
        assert len(value) == 1, f'Expected amino acid string, got: "{value}"'
        position = self._parse_position(item)
        assert position.scheme == self.scheme, f'Expected {self.scheme} scheme, got {position.scheme}'
        assert position.cdr_definition == self.cdr_definition, f'Expected {self.cdr_definition} CDR definition, got {position.cdr_definition}, please pass a Position object with correct CDR definition'
        region_dict = self.regions[position.get_region()]
        if position in region_dict:
            # position is already present, this will replace it
            region_dict[position] = value
        else:
            # make sure to keep the order of positions
            region_dict[position] = value
            for pos in sorted(region_dict):
                region_dict.move_to_end(pos)

    def __getitem__(self, item):
        if isinstance(item, slice):
            if item.step is not None and item.step != 1:
                raise IndexError(f'Slicing with step != 1 is not implemented, got: {item}')
            return self.slice(start=item.start, stop=item.stop)
        pos = self._parse_position(item)
        return self.positions[pos]

    def __len__(self):
        return len(self.positions)

    def __hash__(self):
        return hash(self.positions)

    def __eq__(self, other):
        """Check chain equality. Only checks scheme, aligned sequence and tail sequence, ignores name, metadata and CDR definitions."""
        assert isinstance(other, Chain), f'Can only compare Chain to another Chain, got {type(other)}: {other}'
        return self.positions == other.positions and self.tail == other.tail

    @classmethod
    def to_fasta(cls, chains, path_or_fd, keep_tail=False, description=''):
        """Save multiple chains to FASTA"""
        if isinstance(chains, Chain):
            records = chains.to_seq_record(keep_tail=keep_tail, description=description)
        else:
            records = (chain.to_seq_record(keep_tail=keep_tail, description=description) for chain in chains)
        return SeqIO.write(records, path_or_fd, 'fasta-2line')

    @classmethod
    def from_fasta(cls, path_or_handle, scheme, cdr_definition=None, as_series=False, as_generator=False, **kwargs) -> Union[List['Chain'], pd.Series, Generator['Chain', None, None]]:
        """Read multiple chains from FASTA"""
        generator = (cls(record.seq, name=record.name, scheme=scheme, cdr_definition=cdr_definition, **kwargs)
                     for record in SeqIO.parse(path_or_handle, 'fasta'))
        if as_generator:
            return generator
        chains = list(generator)
        if as_series:
            return pd.Series(chains, index=[c.name for c in chains])
        return chains

    def to_seq_record(self, keep_tail=False, description=''):
        """Create BioPython SeqRecord object from this Chain"""
        if not self.name:
            raise ValueError('Name needs to be present to convert to a SeqRecord')
        seq = Seq(self.seq + self.tail if keep_tail else self.seq)
        return SeqRecord(seq, id=self.name, description=description)

    @classmethod
    def to_anarci_csv(cls, chains: List['Chain'], path):
        """Save multiple chains to ANARCI-like CSV"""
        df = cls.to_dataframe(chains)
        df.to_csv(path)

    @classmethod
    def to_dataframe(cls, chains: List['Chain']):
        """Produce a Pandas dataframe with aligned chain sequences in the columns

        Note: Contains only positions (columns) that are present in the provided chains,
        so number of columns can differ based on the input.
        """
        series_list = [chain.to_series() for chain in chains]

        # Each chain can have a different set of positions
        # so we need to sort the columns to make sure they are in the right order
        # this is using the correct Position sorting
        columns = set(c for series in series_list for c in series.index)
        prop_columns = [c for c in columns if not isinstance(c, Position)]
        position_columns = sorted([c for c in columns if isinstance(c, Position)])
        # Columns can come from K and L chain, so we need to convert them to string and remove duplicates here
        position_columns_str = pd.Series(
            [pos.format(chain_type=False) for pos in position_columns]
        ).drop_duplicates().to_list()

        # Get full list of string columns
        columns_str = prop_columns + position_columns_str

        # Reindex each series using ordered list of string columns
        series_list_ordered = []
        for series in series_list:
            series.index = series.index.map(lambda pos: pos.format(chain_type=False))
            series_list_ordered.append(series.reindex(columns_str))

        df = pd.DataFrame(series_list_ordered)[columns_str].fillna('-')
        df.index.name = 'Id'

        return df

    def to_series(self):
        props = {
            'chain_type': self.chain_type,
            'species': self.species
        }
        return pd.Series({**props, **self.positions}, name=self.name)

    @classmethod
    def from_series(cls, series, scheme, cdr_definition=None) -> 'Chain':
        chain_type = series['chain_type']
        species = series.get('species')
        position_index = [c for c in series.index if c[:1].isnumeric()]
        aa_dict = {Position.from_string(pos, chain_type=chain_type, scheme=scheme): aa
                   for pos, aa in series[position_index].items() if aa != '-' and not pd.isna(aa)}
        return cls(sequence=None, aa_dict=aa_dict, name=series.name, scheme=scheme, cdr_definition=cdr_definition,
                   chain_type=chain_type, species=species, tail='')

    @classmethod
    def from_anarci_csv(cls, path, scheme, cdr_definition=None, as_series=False) -> Union[List['Chain'], pd.Series]:
        df = pd.read_csv(path, index_col=0)
        return cls.from_dataframe(df, scheme=scheme, cdr_definition=cdr_definition, as_series=as_series)

    @classmethod
    def from_dataframe(cls, df, scheme, cdr_definition=None, as_series=False) -> Union[List['Chain'], pd.Series]:
        chains = [cls.from_series(series, scheme=scheme, cdr_definition=cdr_definition) for i, series in df.iterrows()]
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
    
    def print(self, method='wide', **kwargs):
        """Print string representation using :meth:`Chain.format`

        By default, produces "wide" format with sequence on first line and CDR regions higlighted with ``^`` on second line:

        >>> chain.print()
        QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
                                 ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^

        :param method: use ``"wide"`` for :meth:`Chain.format_wide` or ``"tall"`` for :meth:`Chain.format_tall()`
        """
        print(self.format(method=method, **kwargs))

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
        FR1 H1    Q
        FR1 H2    V
        FR1 H3    Q
        FR1 H4    L
        FR1 H5    Q
        FR1 H6    Q
        FR1 H7    S
        ...
        """
        print(self.format_tall(columns=columns))

    def format_wide(self, numbering=False):
        """Create string with sequence on first line and CDR regions higlighted with `^` on second line

        :param numbering: Add position numbers on top
        :return: formatted string
        """
        lines = []
        if numbering:

            first_order = ''
            prev_number = None
            after_double_digit = False
            for pos in self.positions:
                number = str(pos.number // 10)
                if number != prev_number:
                    if after_double_digit:
                        # Special case: when double digits follow another double digits, do not print the first digit
                        number = number[1:]
                    first_order += number
                    if len(number) > 1:
                        after_double_digit = True
                else:
                    if after_double_digit:
                        # Special case: After 10, 11, etc, skip adding the space
                        after_double_digit = False
                    else:
                        first_order += ' '
                prev_number = number

            lines.append(first_order)
            lines.append(''.join(str(pos.number % 10) for pos in self.positions))
            letters = ''.join(pos.letter or ' ' for pos in self.positions)
            if letters.strip():
                lines.append(letters)
        lines.append(self.seq)
        if self.cdr_definition == 'kabat':
            lines.append(''.join('^' if pos.is_in_cdr() else ("°" if pos.is_in_vernier() else ' ') for pos in self.positions))
        else:
            lines.append(''.join('^' if pos.is_in_cdr() else ' ' for pos in self.positions))
        return '\n'.join(lines)

    def print_wide(self, numbering=False):
        """Print string representation using :meth:`Chain.format_wide`

        >>> chain.print_wide()
        QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
                                 ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^
        """
        print(self.format_wide(numbering=numbering))

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

        unique_cdr_definitions = set(pos.cdr_definition for pos_dict in pos_dicts for pos in pos_dict.keys())
        assert len(unique_cdr_definitions) <= 1, f'Aligned chains should use the same CDR definitions, got: {unique_cdr_definitions}'

        shared_pos = sorted(set(pos for pos_dict in pos_dicts for pos in pos_dict.keys()))
        residues = [tuple(pos_dict.get(pos, '-') for pos_dict in pos_dicts) for pos in shared_pos]
        return Alignment(shared_pos, residues, chain_type=self.chain_type, scheme=self.scheme, names=[self.name] + [chain.name for chain in other])

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

        start = self._parse_position(start, allow_raw=allow_raw) if start is not None else None
        stop = self._parse_position(stop, allow_raw=allow_raw) if stop is not None else None

        for i, (pos, aa) in enumerate(positions.items()):
            if start is not None and pos < start:
                continue
            if stop is not None and (pos > stop or (not stop_inclusive and pos >= stop)):
                break
            aa_dict[pos] = replace_seq[i] if replace_seq is not None else aa

        return Chain(
            sequence=None,
            aa_dict=aa_dict,
            name=self.name,
            scheme=self.scheme,
            chain_type=self.chain_type,
            cdr_definition=self.cdr_definition,
            tail=self.tail,
            species=self.species,
            v_gene=self.v_gene,
            j_gene=self.j_gene
        )

    def renumber(self, scheme=None, cdr_definition=None, allowed_species=None, use_anarcii=False, anarcii_args=None):
        """Return copy of this chain aligned using a different numbering scheme or CDR definition

        :param scheme: Change numbering scheme: One of ``imgt``, ``chothia``, ``kabat``, ``aho``.
        :param cdr_definition: Change CDR definition scheme: One of ``imgt``, ``chothia``, ``kabat``, ``north``.
        :param allowed_species: ``None`` to allow all species, or one or more of: ``'human', 'mouse','rat','rabbit','rhesus','pig','alpaca'``
        :param use_anarcii: Use ANARCII (2.0) for numbering, otherwise use ANARCI (1.0)
        :param anarcii_args: Keyword arguments for Anarcii()
        """

        return Chain(
            self.seq + self.tail,
            name=self.name,
            allowed_species=allowed_species,
            scheme=scheme or self.scheme,
            cdr_definition=cdr_definition or scheme or self.cdr_definition,
            assign_germline=self.v_gene is not None,
            use_anarcii=use_anarcii,
            anarcii_args=anarcii_args,
        )

    def graft_cdrs_onto(self, other: 'Chain', backmutate_vernier=False, backmutate_vhh_hallmark=False, backmutations: List[Union['Position',str]] = [], name: str = None) -> 'Chain':
        """Graft CDRs from this Chain onto another chain

        :param other: Chain to graft CDRs into (source of frameworks and tail sequence)
        :param backmutate_vernier: Also graft all Kabat Vernier positions from this chain (perform backmutations)
        :param backmutate_vhh_hallmark: Also graft Kabat 37, 44, 45, and 47 from this chain (VHH Hallmark residues)
        :param backmutations: List of positions that should additionally be grafted from this chain (str or or :class:`Position`)
        :param name: Name of new Chain. If not provided, use name of this chain.
        :return: Chain with CDRs grafted from this chain and frameworks from the given chain
        """
        assert self.scheme == other.scheme, \
            f'Sequences need to have the same numbering scheme, got {self.scheme} and {other.scheme}'
        assert self.cdr_definition == other.cdr_definition, \
            f'Sequences need to have the same CDR definition, got {self.cdr_definition} and {other.cdr_definition}'
        assert self.chain_type == other.chain_type, \
            f'Sequences need to have the same chain type, got {self.chain_type} and {other.chain_type}'

        backmutations = [self._parse_position(pos) for pos in backmutations]

        grafted_dict = {pos: aa for pos, aa in other if not pos.is_in_cdr()}
        for pos, aa in self:
            if aa == 'X':
                # Do not preserve X placeholder positions, use the template instead
                continue
            if pos.is_in_cdr() \
                    or (backmutate_vernier and pos.is_in_vernier()) \
                    or (backmutate_vhh_hallmark and pos.is_vhh_hallmark()) \
                    or pos in backmutations \
                    or grafted_dict.get(pos) == 'X':
                grafted_dict[pos] = aa

        return Chain(sequence=None, aa_dict=grafted_dict, name=name or self.name, chain_type=self.chain_type,
                     scheme=self.scheme, cdr_definition=self.cdr_definition, tail=other.tail,
                     v_gene=other.v_gene, j_gene=other.j_gene)

    def graft_cdrs_onto_human_germline(self, v_gene=None, j_gene=None,
                                       backmutate_vernier=False, backmutations: List[Union['Position',str]] = []):
        """Graft CDRs from this Chain onto the nearest human germline sequence

        :param v_gene: Use defined V germline allele (e.g. IGHV1-18*01), gene (e.g. IGHV1-18) or family (e.g. IGHV1)
        :param j_gene: Use defined J germline allele (e.g. IGHJ1*01) or gene (e.g. IGHJ1)
        :param backmutate_vernier: Also graft all Kabat Vernier positions from this chain (perform backmutations)
        :param backmutations: List of positions that should additionally be grafted from this chain (str or or :class:`Position`)
        :return: Chain with CDRs grafted from this chain and frameworks from TODO
        """
        germline_chain = self.find_merged_human_germline(v_gene=v_gene, j_gene=j_gene)

        if self.scheme != 'imgt' or self.cdr_definition != 'imgt':
            germline_chain = germline_chain.renumber(self.scheme, self.cdr_definition)

        return self.graft_cdrs_onto(germline_chain, backmutate_vernier=backmutate_vernier, backmutations=backmutations)

    def _parse_position(self, position: Union[int, str, 'Position'], allow_raw=False):
        """Create :class:`Position` key object from string or int.

        Note: The position should only be used for indexing, CDR definition is not preserved!

        :param position: Numeric or string position representation
        :param allow_raw: Also allow unaligned numeric (int) indexing from 0 to length of sequence - 1
        :return: new Position object, should only be used for indexing, CDR definition is not preserved!
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
        return self.get_position_by_raw_index(position)

    def get_position_by_raw_index(self, index):
        """Get Position object at corresponding raw numeric position"""
        return list(self.positions.keys())[index]

    def find_human_germlines(self, limit=10, v_gene=None, j_gene=None, unique=True, rank_by_frameworks=False) -> Tuple[List['Chain'], List['Chain']]:
        """Find most identical V and J germline sequences based on IMGT alignment

        :param limit: Number of best matching germlines to return
        :param v_gene: Filter germlines to specific V gene name
        :param j_gene: Filter germlines to specific J gene name
        :param unique: Skip germlines with duplicate amino acid sequence
        :param rank_by_frameworks: Prioritize framework sequence identity when finding nearest matches
        :return: list of top V chains, list of top J chains
        """
        from abnumber.germlines import get_imgt_v_chains, get_imgt_j_chains

        chain = self if self.scheme == 'imgt' and self.cdr_definition == 'imgt' else self.renumber('imgt')
        v_chains = list(get_imgt_v_chains(chain.chain_type).values())
        j_chains = list(get_imgt_j_chains(chain.chain_type).values())

        if v_gene:
            if v_gene.startswith('IGKV') and self.chain_type == 'L':
                raise NotImplementedError('Cannot graft lambda chain into kappa chain')
            if v_gene.startswith('IGLV') and self.chain_type == 'K':
                raise NotImplementedError('Cannot graft kappa chain into lambda chain')
            v_chains = [chain for chain in v_chains if chain.name.startswith(v_gene)]
            if not v_chains:
                print('Available V genes:', get_imgt_v_chains(chain.chain_type).keys())
                raise ValueError(f'No V genes found for "{chain.chain_type}" chain gene name "{v_gene}"')

        if j_gene:
            j_chains = [chain for chain in j_chains if chain.name.startswith(j_gene)]
            if not j_chains:
                print('Available J genes:', get_imgt_j_chains(chain.chain_type).keys())
                raise ValueError(f'No J genes found for "{chain.chain_type}" chain gene name "{j_gene}"')
            
        if unique:
            v_chains = _get_unique_chains(v_chains)
            j_chains = _get_unique_chains(j_chains)

        v_alignments = [chain.align(germline) for germline in v_chains]
        v_ranks = np.array([-alignment.num_identical(ignore_cdrs=rank_by_frameworks) - (alignment.num_similar() * 0.01) for alignment in v_alignments]).argsort(kind='stable')[:limit]
        top_v_chains = [v_chains[r] for r in v_ranks]

        j_alignments = [chain.align(germline) for germline in j_chains]
        j_ranks = np.array([-alignment.num_identical(ignore_cdrs=rank_by_frameworks) - (alignment.num_similar() * 0.01) for alignment in j_alignments]).argsort(kind='stable')[:limit]
        top_j_chains = [j_chains[r] for r in j_ranks]

        return top_v_chains, top_j_chains

    def find_merged_human_germline(self, top=0, v_gene=None, j_gene=None, rank_by_frameworks=False) -> 'Chain':
        """Find n-th most identical V and J germline sequence based on IMGT alignment and merge them into one Chain

        :param top: Return top N most identical germline (0-indexed)
        :param v_gene: Filter germlines to specific V gene name
        :param j_gene: Filter germlines to specific J gene name
        :param rank_by_frameworks: Prioritize framework sequence identity when finding nearest matches
        :return: merged germline sequence Chain object
        """
        v_chains, j_chains = self.find_human_germlines(limit=top+1, v_gene=v_gene, j_gene=j_gene, rank_by_frameworks=rank_by_frameworks)
        v_chain = v_chains[top]
        j_chain = j_chains[top]
        return v_chain.merge(j_chain)

    def merge(self, other: 'Chain') -> 'Chain':
        """Merge with another chain to create a new Chain object, positions from this chain taking priority"""

        merged_dict = {
            **{pos: aa for pos, aa in other},
            **{pos: aa for pos, aa in self}
        }

        return Chain(
            sequence=None,
            aa_dict=merged_dict,
            chain_type=self.chain_type,
            scheme='imgt',
            tail='',
            name=f'{self.name or ""} {other.name or ""}'.strip() or None
        )

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

        Region is an uppercase string, one of: ``"FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"``

        :return: Dictionary of Region name -> Dictionary of (:class:`Position` -> Amino acid)
        """
        return OrderedDict(
            FR1=self.fr1_dict,
            CDR1=self.cdr1_dict,
            FR2=self.fr2_dict,
            CDR2=self.cdr2_dict,
            FR3=self.fr3_dict,
            CDR3=self.cdr3_dict,
            FR4=self.fr4_dict
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
    def fr1_seq(self):
        """Unaligned string representation of the Framework 1 region sequence"""
        return ''.join(self.fr1_dict.values())

    @property
    def cdr1_seq(self):
        """Unaligned string representation of the CDR 1 region sequence"""
        return ''.join(self.cdr1_dict.values())

    @property
    def fr2_seq(self):
        """Unaligned string representation of the Framework 2 region sequence"""
        return ''.join(self.fr2_dict.values())

    @property
    def cdr2_seq(self):
        """Unaligned string representation of the CDR 2 region sequence"""
        return ''.join(self.cdr2_dict.values())

    @property
    def fr3_seq(self):
        """Unaligned string representation of the Framework 3 region sequence"""
        return ''.join(self.fr3_dict.values())

    @property
    def cdr3_seq(self):
        """Unaligned string representation of the CDR 3 region sequence"""
        return ''.join(self.cdr3_dict.values())

    @property
    def fr4_seq(self):
        """Unaligned string representation of the Framework 4 region sequence"""
        return ''.join(self.fr4_dict.values())


class RawChainAccessor:
    def __init__(self, chain: Chain):
        self.chain = chain

    def __getitem__(self, item):
        if isinstance(item, slice):
            if item.step is not None and item.step != 1:
                raise IndexError(f'Slicing with step != 1 is not implemented, got: {item}')
            if item.start is not None and not is_integer(item.start):
                raise IndexError(f'Expected int start index for chain.raw, got {type(item.start)}: {item.start}')
            if item.stop is not None and not is_integer(item.stop):
                raise IndexError(f'Expected int end index for chain.raw, got {type(item.stop)}: {item.stop}')
            return self.chain.slice(start=item.start, stop=item.stop, stop_inclusive=False, allow_raw=True)
        if not is_integer(item):
            raise IndexError(f'Expected int indexing for chain.raw, got {type(item)}: {item}')
        pos = self.chain.get_position_by_raw_index(item)
        return self.chain[pos]




