from abnumber import *
from abnumber import Chain as OldChain


class Chain(OldChain):
    def __init__(self, sequence, scheme, cdr_definition=None, name=None, assign_germline=False, allowed_species=None, use_anarcii=True, anarcii_args=None, **kwargs):
        super().__init__(sequence, scheme, cdr_definition=cdr_definition, name=name, assign_germline=assign_germline, allowed_species=allowed_species, use_anarcii=use_anarcii, anarcii_args=anarcii_args, **kwargs)

    @classmethod
    def batch(cls, seq_dict: dict, scheme: str, cdr_definition=None, assign_germline=False, allowed_species=None, multiple_domains=False, use_anarcii=True, anarcii_args=None):
        return super().batch(seq_dict, scheme, cdr_definition=cdr_definition, assign_germline=assign_germline, allowed_species=allowed_species, multiple_domains=multiple_domains, use_anarcii=use_anarcii, anarcii_args=anarcii_args)

    def renumber(self, scheme=None, cdr_definition=None, allowed_species=None, use_anarcii=True, anarcii_args=None):
        return super().renumber(scheme=scheme, cdr_definition=cdr_definition, allowed_species=allowed_species, use_anarcii=use_anarcii, anarcii_args=anarcii_args)
