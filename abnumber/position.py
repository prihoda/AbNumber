import copy
from typing import List, Union

from abnumber.common import _validate_chain_type, SCHEME_POSITION_TO_REGION, SCHEME_VERNIER, POS_REGEX, SCHEME_HALLMARK


class Position:
    """Numbered position using a given numbering scheme

    Used as a key to store Position -> Amino acid information.

    Position objects are sortable according to the schema simply using ``sorted()``.
    """
    def __init__(self, chain_type: str, number: int, letter: str, scheme: str):
        _validate_chain_type(chain_type)
        self.chain_type: str = chain_type
        self.number: int = int(number)
        self.letter: str = letter.strip()
        self.scheme: str = scheme
        self.cdr_definition: str = self.scheme
        self.cdr_definition_position: int = self.number

    def copy(self):
        return copy.copy(self)

    def _key(self):
        # Note: We are not including chain_type, but just Heavy/Light flag, to keep Kappa and Lambda chain positions equal
        return self.chain_type_prefix(), self.number, self.letter, self.scheme

    def __repr__(self):
        return f'{self.chain_type_prefix()}{self.number}{self.letter} ({self.scheme})'

    def __str__(self):
        return self.format()

    def set_cdr_definition(self, cdr_definition: str, cdr_definition_position: int):
        assert cdr_definition is not None, 'cdr_definition is required'
        assert cdr_definition_position is not None, 'cdr_definition_position is required'
        self.cdr_definition = cdr_definition
        self.cdr_definition_position = cdr_definition_position

    def format(self, chain_type=True, region=False, rjust=False, ljust=False, fillchar=' '):
        """Format Position to string

        :param chain_type: Add chain type prefix (H/L)
        :param region: Add region prefix (FR1, CDR1, ...)
        :param rjust: Align text to the right
        :param ljust: Align text to the left
        :param fillchar: Characer to use for alignment padding
        :return: formatted string
        """
        formatted = f'{self.number}{self.letter}'
        if chain_type:
            formatted = f'{self.chain_type_prefix()}{formatted}'
        if region:
            formatted = f'{self.get_region()} {formatted}'
        just = 4 + 1* int(chain_type) + 5 * int(region)
        if rjust:
            formatted = formatted.rjust(just, fillchar)
        if ljust:
            formatted = formatted.ljust(just, fillchar)
        return formatted

    def __hash__(self):
        return self._key().__hash__()

    def __eq__(self, other):
        return isinstance(other, Position) and self._key() == other._key()

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
        letter_ord = 0
        remaining_letters = self.letter
        while remaining_letters:
            letter_ord *= 256
            letter_ord += ord(remaining_letters[0])
            remaining_letters = remaining_letters[1:]
        if self.scheme == 'imgt':
            if self.number in [33, 61, 112]:
                # position 112 is sorted in reverse
                letter_ord = -letter_ord
        elif self.scheme in ['chothia', 'kabat', 'aho']:
            # all letters are sorted alphabetically for these schemes
            pass
        else:
            raise NotImplementedError(f'Cannot compare positions of scheme: {self.scheme}')
        return self.is_heavy_chain(), self.number, letter_ord

    def get_region(self):
        """Get string name of this position's region

        :return: uppercase string, one of: ``"FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"``
        """
        if self.cdr_definition in SCHEME_POSITION_TO_REGION:
            regions = SCHEME_POSITION_TO_REGION[self.cdr_definition]
        else:
            regions = SCHEME_POSITION_TO_REGION[f'{self.cdr_definition}_{self.chain_type}']
        return regions[self.cdr_definition_position]

    def is_in_cdr(self):
        """Check if given position is found in the CDR regions"""
        return self.get_region().lower().startswith('cdr')

    def is_in_vernier(self):
        if self.scheme == 'kabat':
            return self.number in SCHEME_VERNIER.get(f'{self.scheme}_{self.chain_type}', [])
        elif self.cdr_definition == 'kabat':
            return self.cdr_definition_position in SCHEME_VERNIER.get(f'{self.cdr_definition}_{self.chain_type}', [])
        else:
            raise NotImplementedError('Vernier zone identification is currently supported '
                                      f'only with Kabat numbering or CDR definitions, got: {self.scheme}+{self.cdr_definition}')

    def is_vhh_hallmark(self):
        if self.scheme == 'kabat':
            return self.number in SCHEME_HALLMARK.get(f'{self.scheme}_{self.chain_type}', [])
        elif self.cdr_definition == 'kabat':
            return self.cdr_definition_position in SCHEME_HALLMARK.get(f'{self.cdr_definition}_{self.chain_type}', [])
        else:
            raise NotImplementedError('Hallmark zone identification is currently supported '
                                      f'only with Kabat numbering or CDR definitions, got: {self.scheme}+{self.cdr_definition}')

    @classmethod
    def from_string(cls, position, chain_type, scheme):
        """Create Position object from string, e.g. "H5"

        Note that Positions parsed from string do not support separate CDR definitions.
        """
        match = POS_REGEX.match(position.upper())
        _validate_chain_type(chain_type)
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


def sort_positions(positions: List[str], chain_type: str, scheme: str) -> List:
    """Sort position strings to correct order based on given scheme"""
    has_prefix = [p.startswith('H') or p.startswith('L') for p in positions]
    assert all(has_prefix) or not any(has_prefix), 'Inconsistent position prefix'
    has_prefix = all(has_prefix)

    position_objects = [Position.from_string(p, chain_type=chain_type, scheme=scheme) for p in positions]

    return [p.format(chain_type=has_prefix) for p in sorted(position_objects)]
