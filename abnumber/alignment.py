from typing import Union

from abnumber.common import is_similar_residue, is_integer
from abnumber.position import Position
import itertools


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
    def __init__(self, positions, residues, scheme, chain_type, names=None):
        assert isinstance(positions, list), 'Expected list of positions and residues. ' \
                                            'Use chain.align(other) to create an alignment.'
        assert len(positions) == len(residues)
        unique_cdr_definitions = set(pos.cdr_definition for pos in positions)
        assert len(unique_cdr_definitions) <= 1, f'Aligned chains should use the same CDR definitions, got: {unique_cdr_definitions}'
        num_seqs = len(residues[0])
        self.names = names if names is not None else ([None] * num_seqs)
        assert len(self.names) == num_seqs, 'Number of names should match number of sequences'
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

    def __contains__(self, item):
        position = self._parse_position(item)
        assert position.scheme == self.scheme, f'Expected {self.scheme} scheme, got {position.scheme}'
        return position in self.positions

    def __getitem__(self, item):
        if isinstance(item, slice):
            if item.step is not None and item.step != 1:
                raise IndexError(f'Slicing with step != 1 is not implemented, got: {item}')
            return self.slice(start=item.start, stop=item.stop)
        pos = self._parse_position(item)
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

        start = self._parse_position(start, allow_raw=allow_raw) if start is not None else None
        stop = self._parse_position(stop, allow_raw=allow_raw) if stop is not None else None

        new_positions = []
        new_residues = []
        for pos, residues in zip(self.positions, self.residues):
            if start is not None and pos < start:
                continue
            if stop is not None and (pos > stop or (not stop_inclusive and pos >= stop)):
                break
            new_positions.append(pos)
            new_residues.append(residues)

        return Alignment(positions=new_positions, residues=new_residues, scheme=self.scheme, chain_type=self.chain_type, names=self.names)

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
        return self.positions[position]

    def format(self, mark_identity=True, mark_cdrs=True, names=False):
        """Format alignment to string

        :param mark_identity: Add BLAST style middle line showing identity (``|``), similar residue (``+``) or different residue (``.``)
        :param mark_cdrs: Add line highlighting CDR regions using ``^``
        :param names: Add chain.name to the beginning of each sequence
        :return: formatted string
        """

        def _identity_symbol(a, b):
            return '|' if a == b else ('+' if is_similar_residue(a, b) else '.')

        lines = []
        longest_name = max(len(name or '') for name in self.names) if hasattr(self, 'names') and self.names else 0
        for i in range(len(self.residues[0])):
            name = self.names[i] if hasattr(self, 'names') and self.names else None
            if mark_identity and i != 0:
                lines.append((' '*(longest_name+1) if names else '') + ''.join(_identity_symbol(aas[i], aas[i-1]) for pos, aas in self))
            lines.append(((name or '').rjust(longest_name) + ' ' if names else '') + ''.join(aas[i] for pos, aas in self))
        if mark_cdrs:
            if self.positions[0].cdr_definition == 'kabat':
                lines.append((' '*(longest_name+1) if names else '') + ''.join('^' if pos.is_in_cdr() else ("°" if pos.is_in_vernier() else ' ') for pos in self.positions))
            else:
                lines.append((' '*(longest_name+1) if names else '') + ''.join('^' if pos.is_in_cdr() else ' ' for pos in self.positions))
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

    def num_identical(self, ignore_cdrs = False):
        """Get number of positions with identical residues

        :param ignore_cdrs: Ignore CDR regions when counting identical residues
        """
        if ignore_cdrs:
            return sum(len(set(aas)) == 1 for pos, aas in zip(self.positions, self.residues) if not pos.is_in_cdr())
        return sum(len(set(aas)) == 1 for aas in self.residues)

    def num_similar(self):
        """Get number of positions with similar residues based on BLOSUM62"""
        return sum(len(set(aas)) == 1 or all(is_similar_residue(a, b) for a, b in itertools.combinations(set(aas), 2))
                   for aas in self.residues)

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
            if item.start is not None and not is_integer(item.start):
                raise IndexError(f'Expected int start index for alignment.raw, got {type(item.start)}: {item.start}')
            if item.stop is not None and not is_integer(item.stop):
                raise IndexError(f'Expected int end index for alignment.raw, got {type(item.stop)}: {item.stop}')
            return self.alignment.slice(start=item.start, stop=item.stop, stop_inclusive=False, allow_raw=True)
        if not is_integer(item):
            raise IndexError(f'Expected int indexing for alignment.raw, got {type(item)}: {item}')
        pos = self.alignment.positions[item]
        return self.alignment[pos]
