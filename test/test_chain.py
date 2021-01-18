import pytest
from abnumber import Chain, ChainParseError
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


@pytest.mark.parametrize("scheme", ['chothia', 'kabat', 'imgt'])
def test_invalid_chain_raises_error(scheme):
    with pytest.raises(ChainParseError):
        Chain('AAA', scheme=scheme)


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

    assert Chain(seq, scheme=scheme, cdr_definition=cdr_definition) == Chain(seq, scheme='imgt').renumber(scheme)
    assert Chain(seq, scheme=scheme, cdr_definition=cdr_definition) == Chain(seq, scheme='chothia').renumber(scheme)
    assert Chain(seq, scheme=scheme, cdr_definition=cdr_definition) == Chain(seq, scheme='kabat').renumber(scheme)

    assert Chain(seq, scheme=scheme, cdr_definition=cdr_definition).cdr3_seq == Chain(seq, scheme=scheme, cdr_definition='imgt').renumber(cdr_definition=cdr_definition).cdr3_seq
    assert Chain(seq, scheme=scheme, cdr_definition=cdr_definition).cdr3_seq == Chain(seq, scheme=scheme, cdr_definition='chothia').renumber(cdr_definition=cdr_definition).cdr3_seq
    assert Chain(seq, scheme=scheme, cdr_definition=cdr_definition).cdr3_seq == Chain(seq, scheme=scheme, cdr_definition='kabat').renumber(cdr_definition=cdr_definition).cdr3_seq


def test_multiple_domains():
    vh = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS'
    vl = 'ELVMTQSPSSLSASVGDRVNIACRASQGISSALAWYQQKPGKAPRLLIYDASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFAIYYCQQFNSYPLTFGGGTKVEIK'
    with pytest.raises(ChainParseError):
        chain = Chain(vh + vl, 'imgt')