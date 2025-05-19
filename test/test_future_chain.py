import pytest
from abnumber.future import Chain, ChainParseError, Position
import numpy as np

def test_light_chain_from_str():
    var = 'ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIK'
    tail = 'RTV'
    chain = Chain(var + tail, scheme='imgt')
    assert chain.seq == var
    assert chain.tail == tail

    assert chain == Chain(var + tail, scheme='imgt')
    assert chain != Chain(var + 'AAA', scheme='imgt')

    expected_format = '''
ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIK
                          ^^^^^^                 ^^^                                    ^^^^^^^^^          
    '''

    assert chain.format().strip() == expected_format.strip()

    assert chain.cdr1_seq == 'QGISSA'
    assert chain.cdr2_seq == 'DAS'
    assert chain.cdr3_seq == 'QQFNSYPLT'


def test_heavy_chain_from_str():
    var = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS'
    tail = 'AKTTAPSVYPLA'
    chain = Chain(var + tail, scheme='imgt')
    assert chain.seq == var
    assert chain.tail == tail

    assert chain == Chain(var + tail, scheme='imgt')
    assert chain != Chain(var + 'AAA', scheme='imgt')

    expected_format = '''
QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
                         ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^           
'''

    assert chain.format().strip() == expected_format.strip()

    assert chain.cdr1_seq == 'GYTFTRYT'
    assert chain.cdr2_seq == 'INPSRGYT'
    assert chain.cdr3_seq == 'ARYYDDHYCLDY'


def test_heavy_chain_slice():
    seq = 'EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSS'
    chain = Chain(seq, scheme='imgt')

    assert str(chain['4']) == 'L'
    assert str(chain[:'5']) == 'EVQLQ'
    assert str(chain['H111':]) == 'DDERGHYCLDYWGQGTTLTVSS'
    assert str(chain['3':'2']) == ''
    assert str(chain['1000':'2000']) == ''

    assert str(chain.raw[3]) == 'L'
    assert str(chain.raw[:5]) == 'EVQLQ'
    numpy_int = np.array([5], dtype=np.int64)[0]
    assert str(chain.raw[:numpy_int]) == 'EVQLQ'
    assert str(chain.raw[numpy_int]) == 'Q'

    expected_slice = ', '.join([f'{p}={aa}' for p, aa in chain['111':'112']])
    assert expected_slice == 'H111=D, H111A=D, H111B=E, H112B=R, H112A=G, H112=H'


def test_heavy_alignment_slice():
    seq = 'EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSS'
    chain = Chain(seq, scheme='imgt')
    seq2 = 'QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYDDYLDRWGQGTTLTVSSAKTTAP'
    chain2 = Chain(seq2, scheme='imgt')
    alignment = chain.align(chain2)

    expected_format = '''
EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPS-RGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSS
+|||.||||||.||||+|||||||||||.||||||||||||||||+||||||||.|.||||||||||||||||||||||||||.+||||||||||||||||..........||.|||||||||||
QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYD-------DYLDRWGQGTTLTVSS
                         ^^^^^^^^                 ^^^^^^^^^                                      ^^^^^^^^^^^^^^^^^
'''
    assert expected_format.strip() == alignment.format().strip()

    assert alignment['111B'] == ('E', '-')
    assert alignment.raw[105] == ('E', '-')


def test_same_scheme_different_cdr_definition_positions_are_equal():
    seq = 'EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSS'
    chain1 = Chain(seq, scheme='imgt', cdr_definition='imgt')
    chain2 = Chain(seq, scheme='imgt', cdr_definition='kabat')

    position1 = list(chain1.positions)[0]
    position2 = list(chain2.positions)[0]
    assert position1 == position2, 'Positions with same numbering scheme should be equal no matter the CDR definition'


def test_position_hashing():
    seq = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS'
    chain_imgt_def = Chain(seq, scheme='imgt', cdr_definition='imgt')
    chain_north_def = Chain(seq, scheme='imgt', cdr_definition='north')

    pos_imgt_def = list(chain_imgt_def.positions)[0]
    pos_north_def = list(chain_north_def.positions)[0]

    assert hash(pos_imgt_def) == hash(pos_north_def), 'Positions hash should only be based on numbering scheme, not CDR definition'
    assert pos_imgt_def == pos_north_def, 'Positions equality should only be based on numbering scheme, not CDR definition'
    assert not(pos_imgt_def > pos_north_def), 'Position sorting should only be based on numbering scheme, not CDR definition'
    assert not(pos_imgt_def < pos_north_def), 'Position sorting should only be based on numbering scheme, not CDR definition'
    assert pos_imgt_def.cdr_definition != pos_north_def.cdr_definition

    kappa_pos = Position(chain_type='K', number=1, letter='', scheme='imgt')
    lambda_pos = Position(chain_type='L', number=1, letter='', scheme='imgt')
    assert kappa_pos == lambda_pos
    assert Position(chain_type='K', number=2, letter='', scheme='imgt') > lambda_pos


@pytest.mark.parametrize("scheme", ['chothia', 'kabat', 'imgt'])
def test_invalid_chain_raises_error(scheme):
    with pytest.raises(ChainParseError):
        Chain('AAA', scheme=scheme)


@pytest.mark.parametrize("scheme", ['chothia', 'kabat', 'imgt'])
def test_chain_with_invalid_chars_raises_error(scheme):
    with pytest.raises(ChainParseError):
        Chain('QVQ LVQ', scheme=scheme)


def test_multiple_chains_raises_error():
    with pytest.raises(ChainParseError):
        Chain('QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSSQVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS', scheme='imgt')


def test_aho_without_cdr_definition_raises_error():
    with pytest.raises(ValueError):
        Chain('QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS', scheme='aho')


def test_alignment_with_different_schemes_raises_error():
    seq = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS'
    chain1 = Chain(seq, scheme='imgt')
    chain2 = Chain(seq, scheme='chothia')
    with pytest.raises(AssertionError):
        chain1.align(chain2)


def test_alignment_with_different_cdr_definitions_raises_error():
    seq = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS'
    chain1 = Chain(seq, scheme='imgt', cdr_definition='imgt')
    chain2 = Chain(seq, scheme='imgt', cdr_definition='chothia')
    with pytest.raises(AssertionError):
        chain1.align(chain2)


@pytest.mark.parametrize("scheme", ['chothia', 'kabat', 'imgt', 'aho'])
@pytest.mark.parametrize("cdr_definition", ['chothia', 'kabat', 'imgt', 'north'])
@pytest.mark.parametrize("seq_cdrs", [
    ('QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS', {
        'chothia': ['GYTFTRY', 'NPSRGY', 'YYDDHYCLDY'],
        'kabat': ['RYTMH', 'YINPSRGYTNYNQKFKD', 'YYDDHYCLDY'],
        'imgt': ['GYTFTRYT', 'INPSRGYT', 'ARYYDDHYCLDY'],
        'north': ['KASGYTFTRYTMH', 'YINPSRGYTN', 'ARYYDDHYCLDY']
    }),
    ('ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIK', {
        'chothia': ['RASQGISSALA', 'DASNLES', 'QQFNSYPLT'],
        'kabat': ['RASQGISSALA', 'DASNLES', 'QQFNSYPLT'],
        'imgt': ['QGISSA', 'DAS', 'QQFNSYPLT'],
        'north': ['RASQGISSALA', 'YDASNLES', 'QQFNSYPLT']
    })
])
def test_scheme_and_cdr_definition(scheme, cdr_definition, seq_cdrs):
    seq, expected_cdrs_by_definition = seq_cdrs
    expected_cdrs = expected_cdrs_by_definition[cdr_definition]

    chain = Chain(
        seq,
        scheme=scheme,
        cdr_definition=cdr_definition
    )
    print('\n' + chain.format())

    assert chain.seq == seq

    assert chain.cdr1_seq == expected_cdrs[0]
    assert chain.cdr2_seq == expected_cdrs[1]
    assert chain.cdr3_seq == expected_cdrs[2]

    merged_cdrs = ''
    for pos, aa in chain:
        if pos.is_in_cdr():
            merged_cdrs += aa

    assert merged_cdrs == ''.join(expected_cdrs), 'CDRs based on Position.is_in_cdr should agree with Chain definitions'

    assert chain[:'5'].seq == chain.seq[:5]
    assert chain.raw[:5].seq == chain.seq[:5]
    assert chain.raw[-5:].seq == chain.seq[-5:]
    assert chain[:'200'].seq == chain.seq

    sliced_chain = chain['6':]
    assert sliced_chain.seq == chain.seq[5:]
    assert sliced_chain.cdr1_seq == expected_cdrs[0]
    assert sliced_chain.cdr2_seq == expected_cdrs[1]
    assert sliced_chain.cdr3_seq == expected_cdrs[2]

@pytest.mark.parametrize("scheme", ['chothia', 'kabat', 'imgt', 'aho'])
@pytest.mark.parametrize("cdr_definition", ['chothia', 'kabat', 'imgt', 'north'])
def test_chain_renumber(scheme, cdr_definition):
    seq = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS'

    assert Chain(seq, scheme=scheme, cdr_definition=cdr_definition) == Chain(seq, scheme='imgt').renumber(scheme, cdr_definition)
    assert Chain(seq, scheme=scheme, cdr_definition=cdr_definition) == Chain(seq, scheme='chothia').renumber(scheme, cdr_definition)
    assert Chain(seq, scheme=scheme, cdr_definition=cdr_definition) == Chain(seq, scheme='kabat').renumber(scheme, cdr_definition)

    assert Chain(seq, scheme=scheme, cdr_definition=cdr_definition).cdr3_seq == Chain(seq, scheme=scheme, cdr_definition='imgt').renumber(cdr_definition=cdr_definition).cdr3_seq
    assert Chain(seq, scheme=scheme, cdr_definition=cdr_definition).cdr3_seq == Chain(seq, scheme=scheme, cdr_definition='chothia').renumber(cdr_definition=cdr_definition).cdr3_seq
    assert Chain(seq, scheme=scheme, cdr_definition=cdr_definition).cdr3_seq == Chain(seq, scheme=scheme, cdr_definition='kabat').renumber(cdr_definition=cdr_definition).cdr3_seq


def test_multiple_domains():
    vh = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS'
    vl = 'ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIK'
    with pytest.raises(ChainParseError):
        chain = Chain(vh + vl, 'imgt')


def test_germline_assignment():
    light_seq = 'ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIK'
    light_chain = Chain(light_seq, scheme='imgt', assign_germline=True)
    assert light_chain.v_gene == 'IGKV1-13*02'
    assert light_chain.j_gene == 'IGKJ4*01'
    light_chain = Chain(light_seq, scheme='imgt', assign_germline=True, allowed_species='mouse')
    assert light_chain.v_gene == 'IGKV11-125*01'
    assert light_chain.j_gene == 'IGKJ1*01'

    heavy_seq = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS'
    heavy_chain = Chain(heavy_seq, scheme='imgt', assign_germline=True)
    assert heavy_chain.v_gene == 'IGHV1-4*01'
    assert heavy_chain.j_gene == 'IGHJ2*01'
    heavy_chain = Chain(heavy_seq, scheme='imgt', assign_germline=True, allowed_species='alpaca')
    assert heavy_chain.v_gene == 'IGHV3S1*01'
    assert heavy_chain.j_gene == 'IGHJ4*01'

    assert heavy_chain.raw[:10].v_gene == heavy_chain.v_gene
    assert heavy_chain.raw[:10].j_gene == heavy_chain.j_gene

    chain_without_germline = Chain(heavy_seq, scheme='imgt', assign_germline=False)
    assert chain_without_germline.v_gene is None
    assert chain_without_germline.j_gene is None


def test_nearest_j_region():
    heavy_seq = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    heavy_chain = Chain(heavy_seq, scheme='imgt')

    nearest_v, nearest_j = heavy_chain.find_human_germlines(1)

    assert nearest_j[0].name == 'IGHJ6*01'


def test_batch():
    chains, errors = Chain.batch({
        'A': 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS',
        'B': 'EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSS',
        'C': 'FOO',
        'D': 'EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSSEVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSS'
    }, scheme='imgt')
    assert len(chains) == 2
    assert chains['A'].raw[0] == 'Q'
    assert chains['B'].raw[0] == 'E'
    assert 'C' not in chains
    assert errors['C'] == 'Variable chain sequence not recognized: "FOO"'
    assert 'D' not in chains
    assert errors['D'] == 'Found 2 antibody domains: "EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSSEVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSS"'


def test_batch_multiple_domains():
    chains, errors = Chain.batch({
        'A': 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS',
        'B': 'EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSS',
        'C': 'FOO',
        'D': 'EVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSSGGGGSQVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYSEDDERGHYCLDYWGQGTTLTVSSS'
    }, scheme='imgt', multiple_domains=True)
    assert len(chains) == 3
    assert len(chains['A']) == 1
    assert chains['A'][0].raw[0] == 'Q'
    assert len(chains['B']) == 1
    assert chains['B'][0].raw[0] == 'E'
    assert 'C' not in chains
    assert errors['C'] == 'Variable chain sequence not recognized: "FOO"'
    assert len(chains['D']) == 2
    assert chains['D'][0].raw[0] == 'E'
    assert chains['D'][0].tail == 'GGGGS'
    assert chains['D'][1].raw[0] == 'Q'
    assert chains['D'][1].tail == 'S'


def test_multiple_domains():
    vh = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    vl = 'ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIK'
    chains = Chain.multiple_domains('MELVIS' + vh + 'GGGS' + vl + 'CCC', scheme='imgt')
    assert len(chains) == 2
    assert chains[0].seq == vh
    assert chains[0].tail == 'GGGS'
    assert chains[1].seq == vl
    assert chains[1].tail == 'CCC'

def test_nearest_graft_imgt():
    vh = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    #              IMGT CDRs       ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^
    hu = 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTRYTMHWVRQAPGQGLEWMGIINPSRGYTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    chain = Chain(vh, 'imgt')
    assert chain.graft_cdrs_onto_human_germline().seq == hu


def test_nearest_graft_kabat():
    vh = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    #      °       KABAT CDRs       °°°°^^^^^           °°°^^^^^^^^^^^^^^^^^ ° ° ° °    °                 °°^^^^^^^^^^°
    hu = 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTRYTMHWVRQAPGQGLEWMGYINPSRGYTNYNQKFKDRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    chain = Chain(vh, 'kabat')
    assert chain.graft_cdrs_onto_human_germline().seq == hu


def test_imgt_minimal_framework_will_not_contain_x():
    vh = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    #              IMGT CDRs       ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^
    fw = 'EVQLLESGGGLVQPGGSLRLSCAXXXXXXXXXXXXXXRQAPGKGLEWVXXXXXXXXXXXXXDSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYXXXXXXXXXXXXXXXXXGTLVTVSS'
    hu = 'EVQLLESGGGLVQPGGSLRLSCAASGYTFTRYTMHWVRQAPGKGLEWVGYINPSRGYTNYNDSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARYYDDHYCLDYWGQGTLVTVSS'
    #                            **************           *************                                 *****************
    chain = Chain(vh, 'imgt')
    other = Chain(fw, 'imgt')
    assert chain.graft_cdrs_onto(other).seq == hu



def test_nearest_graft_kabat_insertion_is_not_grafted():
    vh = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMPQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    hu = 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTRYTMHWVRQAPGQGLEWMGYINPSRGYTNYNQKFKDRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    chain = Chain(vh, 'kabat')
    assert chain.graft_cdrs_onto_human_germline().seq == hu


def test_contains():
    vh1 = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMPQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    chain1 = Chain(vh1, 'kabat')
    assert '1' in chain1
    vh2 = 'VQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMPQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    chain2 = Chain(vh2, 'kabat')
    assert '1' not in chain2
    assert '1' in chain1.align(chain2)


def test_setitem():
    vh = 'VQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMPQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    chain = Chain(vh, 'kabat')
    assert chain.seq == vh
    assert '1' not in chain
    chain['1'] = 'X'
    assert chain.seq == 'X' + vh
    assert '1' in chain
    chain['3'] = 'E'
    assert chain.seq == 'XVELQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMPQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS'


def test_setitem_north():
    vh = 'VQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMPQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS'
    chain = Chain(vh, 'aho', 'north')
    assert chain.seq == vh
    assert '1' not in chain
    position = Position('H', 1, '', 'aho')
    position.set_cdr_definition('north', 1)
    chain[position] = 'X'
    assert chain.seq == 'X' + vh
    assert '1' in chain
    position = Position('H', 3, '', 'aho')
    position.set_cdr_definition('north', 3)
    chain[position] = 'E'
    assert chain.seq == 'XVELQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMPQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTVTVSS'

